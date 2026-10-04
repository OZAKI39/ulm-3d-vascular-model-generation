"""Run unchanged project tests and record their collection and failure evidence."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

AUDIT = Path(__file__).resolve().parent
PROJECT = AUDIT.parents[1]
PYTHON = Path('/home/lzy/projects/computation_examples/_shared/python_environment/bin/python')


def run_tests(stage):
    import pytest

    collected, outcomes, failures = [], {}, {}

    class Recorder:
        def pytest_collection_finish(self, session):
            collected.extend(item.nodeid for item in session.items)

        def pytest_runtest_logreport(self, report):
            if report.when == 'call' or report.failed or report.skipped:
                outcomes[report.nodeid] = report.outcome
            if report.failed:
                failures[report.nodeid + '::' + report.when] = report.longreprtext

    args = ['tests', '-q', '-p', 'no:cacheprovider',
            '--basetemp', str(Path(os.environ['CLOT_CLEANUP_VALIDATION_TMP']) / (stage.lower() + '_pytest'))]
    start = time.monotonic()
    result = pytest.main(args, plugins=[Recorder()])
    counts = {key: list(outcomes.values()).count(key) for key in ('passed', 'failed', 'skipped')}
    (AUDIT / (stage + '_TESTS.json')).write_text(json.dumps({
        'command': [sys.executable, '-m', 'pytest'] + args,
        'exit_code': int(result), 'elapsed_seconds': time.monotonic() - start,
        'counts': counts, 'collected_count': len(collected),
        'collected_nodeids': collected, 'failures': failures,
        'environment': {k: os.environ.get(k) for k in [
            'PYTHONDONTWRITEBYTECODE', 'NUMBA_CACHE_DIR', 'MPLCONFIGDIR',
            'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']},
    }, ensure_ascii=False, indent=2) + '\n')
    return int(result)


if __name__ == '__main__':
    stage = sys.argv[1]
    if stage == '_pytest':
        raise SystemExit(run_tests(sys.argv[2]))
    assert stage in ('BASELINE', 'FINAL')
    state_path = AUDIT / 'VALIDATION_TEMP.json'
    if stage == 'BASELINE':
        temp = tempfile.mkdtemp(prefix='clot_repository_cleanup_v1_')
        state_path.write_text(json.dumps({'path': temp, 'owned_by_this_task': True}, indent=2) + '\n')
    else:
        temp = json.loads(state_path.read_text())['path']
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1', NUMBA_CACHE_DIR=temp + '/numba',
               MPLCONFIGDIR=temp + '/mpl', CLOT_CLEANUP_VALIDATION_TMP=temp,
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMBA_NUM_THREADS='1')
    env['PYTHONPATH'] = str(PROJECT) + os.pathsep + env.get('PYTHONPATH', '')
    with (AUDIT / (stage + '.log')).open('w') as log:
        proc = subprocess.run([str(PYTHON), str(Path(__file__).resolve()), '_pytest', stage],
                              cwd=PROJECT, env=env, stdout=log, stderr=subprocess.STDOUT)
    path = AUDIT / (stage + '_TESTS.json')
    if path.exists():
        data = json.loads(path.read_text())
        print(json.dumps({'stage': stage, 'exit_code': proc.returncode,
                          'counts': data['counts'], 'collected_count': data['collected_count'],
                          'elapsed_seconds': data['elapsed_seconds']}), flush=True)
    else:
        print('Test invocation failed before generating records; inspect ' + str(AUDIT / (stage + '.log')))
        raise SystemExit(proc.returncode or 2)
