"""Use the existing project supervisor; cumulative non-CFD budget is 90 min."""
from common import *
import importlib.util

if __name__ == '__main__':
    stage = sys.argv[1]
    elapsed = sum(read(p)['elapsed_s'] for p in (HERE/'monitoring').glob('*/EXECUTION.json'))
    remaining = 5400 - elapsed
    if remaining <= 0:
        raise SystemExit('NON_CFD_BUDGET_EXHAUSTED')
    name, attempt = stage, 1
    while (HERE/'monitoring'/name).exists():
        attempt += 1
        name = f'{stage}_attempt{attempt}'
    spec = importlib.util.spec_from_file_location('managed', BASE/'production_workflow/scripts/managed_job.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.ROOT = HERE
    put('TASK.json', dict(name=name, command=[sys.executable, str(HERE/(stage+'.py')), *sys.argv[2:]],
        cwd=str(HERE), max_seconds=min(remaining, 1800), max_tree_RSS_GiB=6,
        max_nonreclaimable_GiB=34, environment={'PYTHONUNBUFFERED':'1','PYTHONDONTWRITEBYTECODE':'1'}))
    mod.main()
    result = read(HERE/'monitoring'/name/'EXECUTION.json')
    (HERE/'TASK.json').unlink()
    print(name, result['exit_code'], result['stop_reason'], result['elapsed_s'])
    if result['exit_code'] != 0 or result['stop_reason']:
        raise SystemExit(1)
