"""One new dt control after the recorded pre-bounce WLC blow-up, within A."""
import copy
import json
from pathlib import Path
import shutil
import time

from py_scripts.fluid_physics.common import now,sha256_file,write_json
from py_scripts.single_rbc_benchmark.workflow import child_env
from .workflow import load_config,paths,require_authorization,short_spec,isolated_build,counted_run


def main():
    c=load_config();b,runs,_=paths(c);auth=require_authorization(c)
    task='A6_continuous_preparation_half_dt';allocation=300
    read=lambda p:json.loads(Path(p).read_text())
    ledger=read(runs/'gpu/budget_ledger.json')
    assert not any(x['task_id']==task or x['status']=='RUNNING' for x in ledger['attempts']), 'NO_AUTOMATIC_RETRY'
    assert all(x['task_id'].startswith('A') for x in ledger['attempts'])
    used=sum(x['charged_s'] for x in ledger['attempts'])
    assert used+allocation<=auth['gpu_allocation_s']['A_short_controls_and_failure_range']
    failure=b/'A5_failure_step_geometry.json';review=b/'runtime_reviews/A5_continuous_preparation_30/review/review.json'
    assert read(failure)['states']['old']['edges_beyond_WLC_limit']>0
    assert read(review)['execution']==read(runs/'gpu/A5_continuous_preparation_30/execution.json')
    spec=short_spec(c,dict(bouncer_policy='shared'),isolated_build(c))
    spec['config']=copy.deepcopy(c);spec['config']['protocol']['relaxation_time']=30.
    spec['config']['dpd']['dt']=.0005
    spec['config']['dpd']['wall_relax_steps']=2000
    spec.update(continuous_observation=True,dt=.0005,prep_steps=60000,sample_steps=1000)
    code=Path(__file__).parent
    artifacts=[code/'continuous_worker.py',code/'continuous_protocol.py',Path(__file__),failure,review,b/'material_matching.json']
    hashes={str(p):sha256_file(p) for p in artifacts}
    plan=dict(recorded_at=now(),task=task,category='A',allocation_s=allocation,used_A_s_before=used,
        authorization_sha256=sha256_file(b/'authorization.json'),artifact_sha256=hashes,spec=spec,
        changed_factor='Halve dt on the ordered, continuous preparation protocol. Double steps for the same t*=30 and the same physical frozen-wall preparation time; material parameters and mass stay unchanged.',
        comparison='A5 is retained and charged. This is a new integration-parameter control, not a retry of A5 or of the archived unpatched half-dt run.',
        purpose='Test whether smaller steps prevent the observed WLC domain violation and pre-bounce geometric blow-up; evaluate original shape criteria even if the solver exits normally.',
        estimate='A5 cost 26.72 seconds through about 6407 steps including setup. 60000 steps at the observed scale plus setup/output fits the 300-second reservation; no endpoint timing claim.',
        decision_rule='Native failure or failed preparation blocks B/C. No candidate or particle deletion, capacity increase, correction, membrane pinning, or automatic repeat.',
        automatic_retry=False,qualified_speedup=None)
    write_json(b/(task+'_launch_plan.json'),plan)
    def create(d):
        start=time.perf_counter()
        for p,h in hashes.items():assert sha256_file(p)==h
        write_json(d/'spec.json',spec)
        shutil.copyfile(code/'continuous_worker.py',d/'worker_source.py')
        shutil.copyfile(code/'continuous_protocol.py',d/'continuous_protocol_source.py')
        write_json(d/'case_creation.json',dict(elapsed_s=time.perf_counter()-start,new_inputs=True))
        return child_env(gpu=True)+['PYTHONPATH='+str(b.parents[2]),'SINGLE_RBC_AUTHORIZED='+spec['plan_sha256'],'RBC_REPAIR_NATIVE_TRACE=1',
            '/usr/bin/mpirun.openmpi','--bind-to','none','-np','2',str(b.parents[2]/'.venv/bin/python'),'-B',
            '-m','py_scripts.single_rbc_repair.continuous_worker','--spec',str(d/'spec.json')]
    directory,execution,cached=counted_run(c,task,'gpu',create,allocation,dict(spec=spec,artifact_sha256=hashes))
    result=dict(recorded_at=now(),task=task,directory=None if directory is None else str(directory),execution=execution,cached=cached,qualified_speedup=None)
    write_json(b/(task+'_execution.json'),result)
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if execution['status']!='COMPLETED':raise SystemExit(2)


if __name__=='__main__':main()
