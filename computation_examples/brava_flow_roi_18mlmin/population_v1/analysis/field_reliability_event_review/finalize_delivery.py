"""Delivery-only integrity and receipt packaging, after all supervised jobs exit."""
from common import *
import datetime
import os
import psutil


def main():
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    own = psutil.Process()
    excluded = {os.getpid(), *(p.pid for p in own.parents())}
    active = []
    for proc in psutil.process_iter(['pid', 'cmdline', 'status']):
        try:
            args = proc.info['cmdline'] or []
            if proc.pid not in excluded and proc.info['status'] != psutil.STATUS_ZOMBIE:
                if any(str(HERE) in arg for arg in args):
                    active.append(dict(pid=proc.pid, command=args))
        except psutil.Error:
            pass
    assert not active, active
    integrity = {}
    for title, manifest in [('this_review', HERE/'INPUT_SHA256.json'),
                            ('prior_original_evidence', OLD/'INPUT_SHA256.json')]:
        expected = read(manifest)
        actual = {name: sha(BASE/name) for name in expected}
        differences = {name: {'expected': digest, 'actual': actual[name]}
                       for name, digest in expected.items() if actual[name] != digest}
        assert not differences, differences
        integrity[title] = dict(count=len(expected), all_unchanged=True,
                               manifest_sha256=sha(manifest), mismatches=differences)
    frozen_sha = sha(HERE/'EVENT_FIELD_PROTOCOL.json')
    assert frozen_sha == '9e1aa92438e86dadaddc79c5d837d73fe558af447cf4d986d56bbdc2d54a5aea'
    put('FINAL_INPUT_INTEGRITY.json', dict(checked_UTC=now, groups=integrity,
                                        frozen_protocol_sha256=frozen_sha))
    jobs = []
    for path in sorted((HERE/'monitoring').glob('*/EXECUTION.json')):
        d = read(path)
        assert d['status'] == 'EXITED'
        live = read(path.parent/'LIVE.json')
        assert live['status'] == 'EXITED'
        log = path.parent/'job.log'
        samples = d['samples']
        outcome = 'COMPLETED_DIAGNOSTIC_STAGE'
        if path.parent.name.startswith('prepare_coordinates') and d['exit_code'] != 0:
            outcome = 'STOPPED_BY_AGENT_VTK_BOUNDS_STALL; resolved using conservative candidates'
        elif path.parent.name == 'sample_fields_attempt2':
            outcome = 'INCOMPLETE_CV_DESPITE_EXIT_ZERO; replaced by corrected M0 attempt3'
        elif path.parent.name == 'spatial_features':
            outcome = 'FAILED_KEY_TYPE; corrected in attempt2'
        jobs.append(dict(name=path.parent.name, status=d['status'], exit_code=d['exit_code'],
            stop_reason=d['stop_reason'], interpreted_outcome=outcome,
            elapsed_s=d['elapsed_s'], completed_UTC=d['completed_UTC'], task=d['task'],
            peak_tree_RSS_GiB=d['peak_tree_RSS_GiB'],
            peak_cgroup_GiB=d['peak_cgroup_GiB'], peak_working_set_GiB=d['peak_working_set_GiB'],
            peak_nonreclaimable_upper_GiB=max((s['nonreclaimable_upper_GiB'] for s in samples), default=0),
            oom_event_increase=max((s['oom_events']-d['initial_memory']['oom_events'] for s in samples), default=0),
            execution_file=str(path.relative_to(HERE)), execution_sha256=sha(path),
            log_file=str(log.relative_to(HERE)), log_sha256=sha(log),
            log_tail=log.read_text()[-2000:], automatic_next_job=False))
    total = sum(j['elapsed_s'] for j in jobs)
    assert total < 5400
    assert read(HERE/'TEST_RESULTS.json')['exit_code'] == 0
    qc = read(HERE/'VISUAL_QC.json')
    for item in qc['figures']:
        assert sha(HERE/item['file']) == item['sha256'], item['file']
    results = read(HERE/'FIELD_RELIABILITY_RESULTS.json')
    put('EXECUTION_RECEIPT.json', dict(status='COMPLETED_ALL_OWN_JOBS_EXITED', completed_UTC=now,
        new_CFD_count=0, reused_CFD_count=3, reused_CFD=NAMES, reference_field_count=2,
        reused_payloads={n: read(HERE/'FIELD_COMPARABILITY.json')['fields'][n]['raw_payload'] for n in NAMES},
        new_mesh_count=0, new_particle_integrations=0, new_population_count=0,
        supervised_diagnostic_budget_s=5400, supervised_elapsed_s_including_failed_attempts=total,
        timing_scope='Sum of recorded supervised jobs; excludes interactive inspection, coding, transfers and final metadata hashing. No new CFD runtime.',
        jobs=jobs, active_project_processes=active, automatic_next_job=False,
        GPU_peak_tensor_allocated_GiB=max(read(HERE/f'{n}_SAMPLING.json')['GPU_peak_allocated_GiB'] for n in NAMES),
        GPU_peak_scope='PyTorch allocator peak in each field-sampling stage, not whole-device or other-project memory.',
        final_tests=read(HERE/'TEST_RESULTS.json'), input_integrity='FINAL_INPUT_INTEGRITY.json',
        outcome=results['status'], physical_validation=results['physical_validation'],
        log_storage='Original monitoring logs retained on server; ignored by Git, with per-job hashes and excerpts here.'))
    ignored = {'monitoring', '__pycache__', 'runtime', '.pytest_cache'}
    files = [p for p in HERE.rglob('*') if p.is_file() and not (set(p.relative_to(HERE).parts) & ignored)
             and p.name not in {'TASK.json', 'OUTPUT_SHA256.json'}]
    for p in files:
        if p.suffix == '.json':
            def reject(value):
                raise ValueError(f'{p}: non-finite JSON literal {value}')
            json.loads(p.read_text(), parse_constant=reject)
    put('OUTPUT_SHA256.json', {str(p.relative_to(HERE)): sha(p) for p in sorted(files)})
    print(json.dumps(dict(files=len(files), jobs=len(jobs), supervised_seconds=total,
                         new_CFD_count=0, active_project_processes=active, input_integrity=integrity)))


if __name__ == '__main__':
    main()
