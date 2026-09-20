"""Small local process runner, with complete logs and measured elapsed time."""
import os
import signal
import subprocess
import time
from pathlib import Path
from .provenance import now, write_json

ROOT = Path(__file__).resolve().parents[2]


def read_resource_log(path):
    """Keep wall-clock and CPU measurements distinct from Python monotonic time."""
    values = {}
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if line.startswith('Elapsed (wall clock)'):
            clock = line.split(': ', 1)[1]
            seconds = 0.
            for part in clock.split(':'):
                seconds = 60*seconds + float(part)
            values.update(wall_elapsed_text=clock, wall_elapsed_s=seconds)
        elif line.startswith('User time (seconds):'):
            values['user_cpu_s'] = float(line.rsplit(':', 1)[1])
        elif line.startswith('System time (seconds):'):
            values['system_cpu_s'] = float(line.rsplit(':', 1)[1])
        elif line.startswith('Maximum resident set size (kbytes):'):
            values['peak_rss_kib'] = int(line.rsplit(':', 1)[1])
    return values


def run_logged(command, label, cwd=None, env=None, timeout=None):
    logs = ROOT / 'logs/sv1'
    logs.mkdir(parents=True, exist_ok=True)
    stdout_path = logs / (label + '.log')
    resource_path = logs / (label + '_resources.txt')
    measured = ['/usr/bin/time', '-v', '-o', str(resource_path), *map(str, command)]
    start = time.monotonic()
    print('Starting', label, flush=True)
    timed_out = False
    with stdout_path.open('w') as stream:
        try:
            process = subprocess.Popen(measured, cwd=cwd or ROOT, env=env, stdout=stream,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            code = None
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    elapsed = time.monotonic() - start
    rss = None
    if resource_path.exists():
        for line in resource_path.read_text().splitlines():
            if 'Maximum resident set size (kbytes):' in line:
                rss = int(line.rsplit(':', 1)[1])
    result = {'timestamp': now(), 'command': list(map(str, command)), 'cwd': str(cwd or ROOT),
              'exit_code': code, 'timed_out': timed_out, 'elapsed_s': elapsed,
              'peak_rss_kib': rss, 'memory_scope': 'GNU time maximum process RSS, not an MPI sum',
              'log': str(stdout_path.relative_to(ROOT)), 'resource_log': str(resource_path.relative_to(ROOT))}
    write_json(logs / (label + '.json'), result)
    print('Finished', label, 'exit', code, 'elapsed', round(elapsed, 3), flush=True)
    return result


def sv_environment():
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
               QT_QPA_PLATFORM='offscreen', OMP_NUM_THREADS='1',
               XDG_CONFIG_HOME=str(ROOT / 'outputs/sv1/app_config'),
               XDG_CACHE_HOME=str(ROOT / 'outputs/sv1/app_cache'))
    env.pop('PYTHONPATH', None)
    return env
