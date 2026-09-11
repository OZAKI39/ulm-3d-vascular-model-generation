"""One evidence-gated continuous observation control, charged to allocation A."""
import json
from pathlib import Path
import shutil
import time

from py_scripts.fluid_physics.common import now, sha256_file, write_json
from py_scripts.single_rbc_benchmark.workflow import child_env
from .workflow import load_config, paths, require_authorization, short_spec, isolated_build, counted_run


def main():
    c=load_config();b,runs,_=paths(c);auth=require_authorization(c)
    task='A4_continuous_native_observation';allocation=90
    read=lambda p:json.loads(Path(p).read_text())
    ledger=read(runs/'gpu/budget_ledger.json')
    assert not any(x['status']=='RUNNING' for x in ledger['attempts'])
    assert not any(x['task_id']==task for x in ledger['attempts']), 'NO_AUTOMATIC_RETRY'
    assert all(x['task_id'].startswith('A') for x in ledger['attempts'])
    used=sum(x['charged_s'] for x in ledger['attempts'])
    assert used+allocation<=auth['gpu_allocation_s']['A_short_controls_and_failure_range']
    prior=[]
    for name in ('A0_shared_installed','A1_independent_installed','A2_shared_ordered_isolated','A3_ordered_instrumented'):
        p=b/'runtime_reviews'/name/'review.json';review=read(p)
        assert review['execution']==read(runs/'gpu'/name/'execution.json')
        prior.append(dict(task=name,review=str(p),sha256=sha256_file(p),status=review['status']))
    library=isolated_build(c)
    spec=short_spec(c,dict(bouncer_policy='shared'),library)
    spec['continuous_observation']=True
    code=Path(__file__).parent
    artifacts=[code/'continuous_worker.py',code/'continuous_protocol.py',Path(__file__),
               b/'run_boundary_audit.json',b/'cpu_continuous_protocol_before_fix.log',b/'cpu_continuous_protocol_after_fix.log']
    hashes={str(p):sha256_file(p) for p in artifacts}
    identity=dict(spec=spec,artifact_sha256=hashes)
    decision=dict(recorded_at=now(),task=task,category='A',allocation_s=allocation,used_A_s_before=used,
                  authorization_sha256=sha256_file(b/'authorization.json'),prior_reviews=prior,
                  changed_factor='Replace repeated run(0)/run(100) observation boundaries with one run(5000), observed by native HDF5 plugins.',
                  unchanged_physics='Same deformable membrane, two fluid species, forces, ordered shared bouncer, dt, geometry and zero-shear endpoint as A3.',
                  purpose='Remove implicit repeated initial classification and force resets; inspect saved physical quality before any long run.',
                  limitations='Native dump synchronization and stochastic trajectories can differ. This is not an exact deterministic race proof or a Gamma=4 result.',
                  automatic_retry=False,artifact_sha256=hashes,spec=spec)
    write_json(b/(task+'_launch_plan.json'),decision)
    def create(d):
        start=time.perf_counter();write_json(d/'spec.json',spec)
        for p,h in hashes.items():assert sha256_file(p)==h
        shutil.copyfile(code/'continuous_worker.py',d/'worker_source.py')
        shutil.copyfile(code/'continuous_protocol.py',d/'continuous_protocol_source.py')
        write_json(d/'case_creation.json',dict(elapsed_s=time.perf_counter()-start,new_inputs=True))
        return child_env(gpu=True)+['PYTHONPATH='+str(b.parents[2]),'SINGLE_RBC_AUTHORIZED='+spec['plan_sha256'],
            'RBC_REPAIR_NATIVE_TRACE=1','/usr/bin/mpirun.openmpi','--bind-to','none','-np','2',str(b.parents[2]/'.venv/bin/python'),'-B',
            '-m','py_scripts.single_rbc_repair.continuous_worker','--spec',str(d/'spec.json')]
    directory,execution,cached=counted_run(c,task,'gpu',create,allocation,identity)
    result=dict(recorded_at=now(),task=task,directory=None if directory is None else str(directory),execution=execution,cached=cached,qualified_speedup=None)
    write_json(b/(task+'_execution.json'),result)
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if execution['status']!='COMPLETED':raise SystemExit(2)


if __name__=='__main__':main()
