"""Bound each native harness/build using the existing project supervisor."""
from common import *
import importlib.util

if __name__=='__main__':
    stage=sys.argv[1]
    if stage not in {'preflight','baseline','build','validate','boundary','finalize'}:
        raise SystemExit('Only explicitly listed non-CFD stages are permitted')
    name=stage;attempt=1
    while (HERE/'monitoring'/name).exists():
        attempt+=1;name=f'{stage}_attempt{attempt}'
    spec=importlib.util.spec_from_file_location('managed',BASE/'production_workflow/scripts/managed_job.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.ROOT=HERE
    put('TASK.json',dict(name=name,command=[sys.executable,str(HERE/(stage+'.py')),*sys.argv[2:]],cwd=str(HERE),
        max_seconds=1800,max_tree_RSS_GiB=6,max_nonreclaimable_GiB=34,
        environment={'PYTHONUNBUFFERED':'1','PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}))
    try:mod.main()
    finally:(HERE/'TASK.json').unlink(missing_ok=True)
    r=read(HERE/'monitoring'/name/'EXECUTION.json')
    print(name,r['exit_code'],r['stop_reason'],r['elapsed_s'])
    if r['exit_code']!=0 or r['stop_reason']:raise SystemExit(1)
