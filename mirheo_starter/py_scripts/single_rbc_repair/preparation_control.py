"""Bounded extension of continuous zero-shear preparation, within allocation A."""
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
    task='A5_continuous_preparation_30';allocation=180
    read=lambda p:json.loads(Path(p).read_text())
    ledger=read(runs/'gpu/budget_ledger.json')
    assert not any(x['task_id']==task or x['status']=='RUNNING' for x in ledger['attempts']), 'NO_AUTOMATIC_RETRY'
    assert all(x['task_id'].startswith('A') for x in ledger['attempts'])
    used=sum(x['charged_s'] for x in ledger['attempts'])
    assert used+allocation<=auth['gpu_allocation_s']['A_short_controls_and_failure_range']
    review=b/'runtime_reviews/A4_continuous_native_observation/review/review.json'
    membership=b/'A4_full_population_membership.json'
    assert all(read(review)['short_screen'].values())
    assert all(x['confirmed_mismatch_count']==0 and x['inner_particles']==642 for x in read(membership)['frames'])
    spec=short_spec(c,dict(bouncer_policy='shared'),isolated_build(c))
    spec['config']=copy.deepcopy(c);spec['config']['protocol']['relaxation_time']=30.
    spec.update(continuous_observation=True,prep_steps=30000,sample_steps=500)
    code=Path(__file__).parent
    artifacts=[code/'continuous_worker.py',code/'continuous_protocol.py',Path(__file__),review,membership,b/'material_matching.json']
    hashes={str(p):sha256_file(p) for p in artifacts}
    plan=dict(recorded_at=now(),task=task,category='A',allocation_s=allocation,used_A_s_before=used,
        authorization_sha256=sha256_file(b/'authorization.json'),artifact_sha256=hashes,spec=spec,
        changed_factor='Extend continuous zero-shear preparation from t*=5 to the predeclared maximum t*=30; no shear, membrane or solvent parameter changes.',
        observation='Native HDF5 every 500 steps (0.5 time units), native collision trace still every 250 steps; one run call. This output cadence change is disclosed and no deterministic trajectory equivalence is claimed.',
        criteria=read(b/'material_matching.json')['common_preparation_criteria'],
        decision_rule='Review geometry, all stored population membership and preparation windows. A quality failure stops further expensive B/C, without automatic retry or threshold changes.',
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
