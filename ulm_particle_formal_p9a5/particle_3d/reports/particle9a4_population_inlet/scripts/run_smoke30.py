"""Thirty accepted P9-A.4 births, frozen P9-A.1 dynamics, isolated server only."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import importlib.util,json,multiprocessing,os,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.population_inlet_p9a4 import load_new_environment,make_source,REPORT_RELATIVE
from particle_3d.continuous_infusion import sha256,canonical_bytes
R=ROOT/REPORT_RELATIVE


def dump(path,value):
 with Path(path).open('xb') as f:f.write(canonical_bytes(value))


def main():
 if not str(ROOT).startswith('/workspace/particle9a4_population_inlet_'):
  raise ValueError('Isolated authorized P9-A.4 server directory required')
 for name,digest in json.loads((ROOT/'deployment_manifest.json').read_text()).items():
  if sha256(ROOT/name)!=digest:raise ValueError('Deployed source/input mismatch: '+name)
 gate=json.loads((R/'data/smoke_gate.json').read_text())
 if not all(gate[k] for k in ('inlet_audit_pass','portable_regression_pass','new_tests_pass','legacy_b_c_pass')):
  raise ValueError('Inlet and regression gates required before smoke')
 contract=json.loads((ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json').read_text())
 ledger_path=R/'data/inlet100k/smoke30_births.json'
 if sha256(ledger_path)!=gate['smoke30_births_sha256']:raise ValueError('Accepted cohort changed')
 cohort=json.loads(ledger_path.read_text());events=cohort['events']
 if len(events)!=30 or [e['particle_id'] for e in events]!=list(range(1,31)):raise ValueError('Only first 30 accepted births authorized')
 out=R/'outputs/smoke30';out.mkdir(parents=True,exist_ok=False)
 env=load_new_environment(ROOT);source=make_source(env,contract)
 if cohort['identity']!=source.identity:raise ValueError('Population/science source identity differs from audit')
 reference=ROOT/'particle_3d/reports/network_derived_flow_mb_validation_v1'
 spec=importlib.util.spec_from_file_location('p9a4_frozen_mb_runner',reference/'scripts/runner.py')
 runner=importlib.util.module_from_spec(spec);sys.modules[spec.name]=runner;spec.loader.exec_module(runner)
 runner.ENV=env;runner.ROOT=ROOT;runner.INITIAL={}
 for e in events:
  sample=env.field.sample(e['birth_center_m']);assert sample.inside_lumen
  runner.INITIAL[str(e['particle_id'])]=dict(position_m=e['birth_center_m'],q=e['q'],velocity_m_s=sample.velocity_m_s.tolist(),omega_s_inv=(.5*sample.vorticity_s_inv).tolist())
 runner.IDENTITY=dict(model='P9A4_POPULATION_WITH_UNCHANGED_P9A1_DYNAMICS',flow_sha256=env.new_flow_sha256,
  source_contract_sha256=sha256(ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json'),
  source_identity=source.identity,cohort_sha256=sha256(ledger_path),selection='FIRST_30_ACCEPTED_SOURCE_EVENTS_NO_OUTLET_CONDITION',paired_old_new_comparison=False)
 dump(out/'initial_states.json',runner.INITIAL);dump(out/'config.json',dict(identity=runner.IDENTITY,workers=6,dt_s=.00025,horizon_s=1.5,initial_velocity='NEW_FIELD_FREE_INITIAL',threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']}))
 t=time.perf_counter();rows=[]
 with ProcessPoolExecutor(max_workers=6,mp_context=multiprocessing.get_context('fork')) as pool:
  tasks=[pool.submit(runner.mb_job,(e,str(out))) for e in events]
  for future in as_completed(tasks):
   row=future.result();rows.append(row);print(json.dumps(row),flush=True)
 dump(out/'completed.json',dict(results=sorted(rows,key=lambda r:r['bubble_id']),runtime_s=time.perf_counter()-t,count=30))
 # Basin labels are computed only after the immutable accepted cohort exists.
 from particle_3d.particle82_point_native import NativePointTracer
 env.native_point=NativePointTracer(env);point_rows=[]
 with ProcessPoolExecutor(max_workers=6,mp_context=multiprocessing.get_context('fork')) as pool:
  tasks=[pool.submit(runner.point_job,(e,str(out/'point'))) for e in events]
  for future in as_completed(tasks):point_rows.append(future.result())
 dump(out/'point_completed.json',dict(count=30,selection='FIRST_30_ACCEPTED_EVENTS_DESCRIPTIVE_ONLY',rows=sorted(point_rows,key=lambda r:r['bubble_id'])))
 print('SMOKE_AND_DESCRIPTIVE_POINT_AUDIT_COMPLETE',flush=True)
if __name__=='__main__':main()
