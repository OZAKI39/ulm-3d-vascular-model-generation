"""Reuse existing project supervisor; no unattended jobs or implicit next stage."""
from common import *
import importlib.util

if __name__=='__main__':
    stage=sys.argv[1]
    elapsed=sum(read(p)['elapsed_s'] for p in (HERE/'monitoring').glob('*/EXECUTION.json'))
    remaining=5400-elapsed
    if remaining<=0: raise SystemExit('DIAGNOSTIC_BUDGET_EXHAUSTED')
    spec=importlib.util.spec_from_file_location('managed',BASE/'production_workflow/scripts/managed_job.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.ROOT=HERE
    name=stage; attempt=1
    while (HERE/'monitoring'/name).exists():
        attempt+=1;name=f'{stage}_attempt{attempt}'
    put('TASK.json',dict(name=name,command=[sys.executable,str(HERE/(stage+'.py')),*sys.argv[2:]],
       cwd=str(HERE),max_seconds=min(remaining,1800),max_tree_RSS_GiB=10,max_nonreclaimable_GiB=34,
       environment={'PYTHONUNBUFFERED':'1','PYTHONDONTWRITEBYTECODE':'1'}))
    module.main()
