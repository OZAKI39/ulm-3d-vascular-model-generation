"""Run one remaining frozen short control after inspecting the prior control.

This is not a retry entry point. A0 remains owned by the original queue;
A1/A2/A3 must have no prior paid attempt and use the unchanged frozen protocol.
"""
import argparse
import json
from pathlib import Path
import shutil
import time

from py_scripts.fluid_physics.common import now,sha256_file,write_json
from py_scripts.single_rbc_benchmark.workflow import child_env
from .workflow import load_config,paths,require_authorization,short_spec,isolated_build,counted_run


def execute(task_id):
    c=load_config();b,runs,_=paths(c);auth=require_authorization(c)
    read=lambda p:json.loads(Path(p).read_text())
    plan=read(b/'frozen_benchmark_plan.json');queue=plan['first_executable_queue']
    index=next(i for i,t in enumerate(queue) if t['id']==task_id)
    if index==0:raise ValueError('A0_NOT_A_RETRY_ENTRY')
    ledger=read(runs/'gpu/budget_ledger.json')
    assert all(a['status']!='RUNNING' for a in ledger['attempts'])
    assert not any(a['task_id']==task_id for a in ledger['attempts']), 'PAID_ATTEMPT_ALREADY_EXISTS_NO_RETRY'
    reviews=[]
    for previous in queue[:index]:
        p=b/'runtime_reviews'/previous['id']/'review.json';review=read(p)
        assert review['task']==previous['id'] and review['execution']==read(runs/'gpu'/previous['id']/'execution.json')
        reviews.append(dict(task=previous['id'],review=str(p),sha256=sha256_file(p),status=review['status']))
    task=queue[index]
    used=sum(a['charged_s'] for a in ledger['attempts'])
    assert used+task['allocation_s']<=360, 'FROZEN_FIRST_QUEUE_CAP_EXCEEDED'
    library=isolated_build(c) if task['library']=='isolated' else None
    spec=short_spec(c,task,library);worker=Path(__file__).parent/'mirheo_worker.py'
    identity=dict(spec=spec,worker_sha256=sha256_file(worker))
    decision=dict(recorded_at=now(),task=task,authorization=str(b/'authorization.json'),authorization_sha256=sha256_file(b/'authorization.json'),previous_control_reviews=reviews,
                  reason='Continue a different, explicitly pre-approved control after reviewing the preserved earlier failure/quality evidence. Do not retry the earlier task or start long B/C experiments.',
                  changed_factor=('independent per-PV bouncer instances in the old library' if index==1 else 'isolated ordered native implementation with shared bouncer' if index==2 else 'native trace enabled on the same ordered implementation'),
                  spec=spec,worker_sha256=sha256_file(worker),launcher_sha256=sha256_file(__file__),first_queue_charged_s_before=used,first_queue_limit_s=360,automatic_retry=False)
    write_json(b/(task_id+'_reviewed_launch_plan.json'),decision)
    def create(d):
        start=time.perf_counter();write_json(d/'spec.json',spec)
        write_json(d/'case_creation.json',dict(elapsed_s=time.perf_counter()-start,new_inputs=True))
        shutil.copyfile(worker,d/'worker_source.py')
        return child_env(gpu=True)+['PYTHONPATH='+str(b.parents[2]),'SINGLE_RBC_AUTHORIZED='+spec['plan_sha256'],
                'RBC_REPAIR_NATIVE_TRACE='+('1' if task.get('trace') else '0'),'/usr/bin/mpirun.openmpi','--bind-to','none','-np','2',str(b.parents[2]/'.venv/bin/python'),'-B','-m','py_scripts.single_rbc_repair.mirheo_worker','--spec',str(d/'spec.json')]
    directory,execution,cached=counted_run(c,task_id,'gpu',create,task['allocation_s'],identity)
    result=dict(recorded_at=now(),task=task_id,directory=None if directory is None else str(directory),execution=execution,cached=cached,attempt_launched=directory is not None,new_repeat=directory is not None and not cached,qualified_speedup=None,next_action='CPU review before the next distinct control')
    write_json(b/(task_id+'_reviewed_execution.json'),result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--task',required=True,choices=['A1_independent_installed','A2_shared_ordered_isolated','A3_ordered_instrumented'])
    result=execute(p.parse_args().task)
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if result['execution']['status']!='COMPLETED':raise SystemExit(2)


if __name__=='__main__':main()
