"""One supervised job with explicit limits, reclaimable-cache accounting and tree stop."""
from pathlib import Path
import datetime
import json
import os
import signal
import subprocess
import time
import psutil

ROOT = Path(__file__).resolve().parents[1]


def put(path, value):
    temporary = path.with_suffix(path.suffix + '.pending')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def cgroup_memory(directory=Path('/sys/fs/cgroup')):
    usage = int((directory / 'memory.current').read_text())
    stats = dict(line.split() for line in (directory / 'memory.stat').read_text().splitlines())
    events = dict(line.split() for line in (directory / 'memory.events').read_text().splitlines())
    inactive = int(stats.get('inactive_file', 0))
    # Active file cache is also reclaimable. Some cgroup kernels do not expose
    # aggregate 'kernel', so conservatively sum its principal components.
    kernel = int(stats.get('kernel', sum(int(stats.get(k, 0)) for k in
                                       ('kernel_stack', 'pagetables', 'percpu', 'sock', 'slab'))))
    nonreclaimable = sum(int(stats.get(k, 0)) for k in
                        ('anon', 'shmem', 'unevictable', 'file_dirty', 'file_writeback')) + kernel
    return dict(cgroup_GiB=usage / 2**30, inactive_file_GiB=inactive / 2**30,
                working_set_GiB=max(0, usage - inactive) / 2**30,
                nonreclaimable_upper_GiB=nonreclaimable / 2**30,
                active_file_GiB=int(stats.get('active_file', 0)) / 2**30,
                anon_GiB=int(stats.get('anon', 0)) / 2**30,
                oom_events=int(events.get('oom', 0)), oom_kill_events=int(events.get('oom_kill', 0)))


def exceeded(task, sample, initial):
    if sample.get('nonreclaimable_upper_GiB', 0) >= task.get('max_nonreclaimable_GiB', float('inf')):
        return 'CGROUP_NONRECLAIMABLE_LIMIT'
    if 'max_working_set_GiB' in task and sample['working_set_GiB'] >= task['max_working_set_GiB']:
        return 'CGROUP_WORKING_SET_LIMIT'
    if 'max_cgroup_GiB' in task and sample['cgroup_GiB'] >= task['max_cgroup_GiB']:
        return 'CGROUP_MEMORY_LIMIT'
    if sample.get('tree_RSS_GiB', 0) >= task.get('max_tree_RSS_GiB', float('inf')):
        return 'PROCESS_TREE_RSS_LIMIT'
    if sample['oom_events'] > initial['oom_events'] or sample['oom_kill_events'] > initial['oom_kill_events']:
        return 'CGROUP_OOM_EVENT'
    return None


def stop_tree(proc, observed, grace_seconds=15):
    # Process identities also cover descendants that create a separate MPI session.
    try:
        for child in psutil.Process(proc.pid).children(recursive=True):
            observed[(child.pid, child.create_time())] = child
    except psutil.Error:
        pass
    alive = []
    for child in observed.values():
        try:
            if child.is_running():
                child.terminate()
                alive.append(child)
        except psutil.Error:
            pass
    if proc.poll() is None:
        proc.terminate()
    # Some container kernels reject psutil's pidfd_open with EINVAL. Process
    # identity/status polling avoids that API and still detects PID reuse.
    deadline = time.monotonic() + grace_seconds
    survivors = alive
    while survivors:
        remaining = []
        for child in survivors:
            try:
                if child.is_running() and child.status() != psutil.STATUS_ZOMBIE:
                    remaining.append(child)
            except psutil.Error:
                pass
        survivors = remaining
        if not survivors or time.monotonic() >= deadline:
            break
        time.sleep(.1)
    for child in survivors:
        try:
            child.kill()
        except psutil.Error:
            pass
    try:
        proc.wait(timeout=3)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()


def main():
    task = json.loads((ROOT / 'TASK.json').read_text())
    out = ROOT / 'monitoring' / task['name']
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    env.update(task.get('environment', {}))
    initial = cgroup_memory()
    start = time.monotonic()
    samples, requested, observed = [], [], {}
    reason = exceeded(task, initial, initial)
    proc = None
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda *_: requested.append('USER_OR_SUPERVISOR_STOP'))
    with (out / 'job.log').open('x') as log:
        if reason is None:
            proc = subprocess.Popen(task['command'], cwd=task.get('cwd', str(ROOT)), env=env,
                                    stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            while True:
                rss = 0
                try:
                    parent = psutil.Process(proc.pid)
                    for child in [parent] + parent.children(recursive=True):
                        try:
                            observed[(child.pid, child.create_time())] = child
                            rss += child.memory_info().rss
                        except psutil.Error:
                            pass
                except psutil.Error:
                    pass
                elapsed = time.monotonic() - start
                row = dict(elapsed_s=elapsed, tree_RSS_GiB=rss / 2**30, **cgroup_memory())
                samples.append(row)
                put(out / 'LIVE.json', dict(status='RUNNING', pid=proc.pid, task=task['name'], **row))
                if requested:
                    reason = requested[0]
                elif elapsed > task['max_seconds']:
                    reason = 'TIME_LIMIT'
                else:
                    reason = exceeded(task, row, initial)
                if reason:
                    stop_tree(proc, observed)
                    break
                if proc.poll() is not None:
                    break
                time.sleep(1)
        else:
            log.write('PRELAUNCH_RESOURCE_GATE: ' + reason + '\n')
    result = dict(status='EXITED', exit_code=proc.wait() if proc else None, stop_reason=reason,
                  pid=proc.pid if proc else None, task=task, elapsed_s=time.monotonic() - start,
                  peak_tree_RSS_GiB=max((row['tree_RSS_GiB'] for row in samples), default=0),
                  peak_cgroup_GiB=max((row['cgroup_GiB'] for row in samples), default=initial['cgroup_GiB']),
                  peak_working_set_GiB=max((row['working_set_GiB'] for row in samples), default=initial['working_set_GiB']),
                  initial_memory=initial, samples=samples,
                  memory_accounting='Both working_set=current-inactive_file and a conservative anon+shmem+kernel+unevictable+dirty/writeback upper estimate are recorded; task selects its limit. Tree RSS and OOM guards remain. No cache flushing or cgroup modification.',
                  completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(), automatic_next_job=False)
    put(out / 'EXECUTION.json', result)
    put(out / 'LIVE.json', {key: value for key, value in result.items() if key != 'samples'})


if __name__ == '__main__':
    main()
