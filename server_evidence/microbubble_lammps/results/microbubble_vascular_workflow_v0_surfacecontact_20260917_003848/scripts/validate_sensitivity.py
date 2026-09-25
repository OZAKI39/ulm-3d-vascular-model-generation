from pathlib import Path
import json,csv,hashlib,collections,math
S=Path(__file__).resolve().parents[1]
def read(path):
 with path.open() as f:return list(csv.DictReader(f))
base=S/'runs/SENS_T1_D15_mpi1';tr0=read(base/'TRAJECTORIES.csv');c0=read(base/'WALL_SURFACE_CONTACT_EVENTS.csv');base_key=lambda r:(r['step'],r['stage'],r['particle_id'])
accepted0={base_key(r):r for r in c0 if r['accepted']=='1'};out=[]
for D in sorted((S/'runs').glob('SENS_*_mpi1')):
 cfg=dict(line.split() for line in (D/'case.cfg').read_text().splitlines() if line and not line.startswith('#'));state=json.loads((D/'RUN_STATE.json').read_text());assert state['time_s']==.51 and state['reason']=='MAX_PHYSICAL_TIME'
 tr=read(D/'TRAJECTORIES.csv');c=read(D/'WALL_SURFACE_CONTACT_EVENTS.csv');accepted={base_key(r):r for r in c if r['accepted']=='1'};assert len(tr)==len(tr0) and accepted.keys()==accepted0.keys()
 maxpos=maxvel=0.;mode_changes=cluster_changes=candidate_changes=0
 for a,b in zip(tr,tr0):
  assert (a['step'],a['time_s'],a['particle_id'])==(b['step'],b['time_s'],b['particle_id'])
  maxpos=max(maxpos,math.sqrt(sum((float(a[k+'_m'])-float(b[k+'_m']))**2 for k in 'xyz')))
 for key,a in accepted.items():
  b=accepted0[key];cluster_changes+=int(a['surface_cluster_count']!=b['surface_cluster_count']);mode_changes+=int(a['constraint_mode']!=b['constraint_mode']);candidate_changes+=int(a['candidate_triangle_count']!=b['candidate_triangle_count']);maxvel=max(maxvel,math.sqrt(sum((float(a['used_v'+k])-float(b['used_v'+k]))**2 for k in 'xyz')))
 assert cluster_changes==0 and mode_changes==0 and maxpos<1e-14 and maxvel<1e-12,(D.name,maxpos,maxvel,cluster_changes,mode_changes)
 out.append(dict(case=D.name,tie_factor=float(cfg['tie_factor']),smooth_deg=float(cfg['smooth_dihedral_threshold_deg']),status='PASS',time=state['time_s'],trajectory_rows=len(tr),accepted_stage_rows=len(accepted),candidate_count_changes=candidate_changes,cluster_count_changes=cluster_changes,mode_changes=mode_changes,max_position_difference_m=maxpos,max_used_velocity_difference_m_s=maxvel,trajectory_sha256=hashlib.sha256((D/'TRAJECTORIES.csv').read_bytes()).hexdigest(),constraint_mode_counts=dict(collections.Counter(a['constraint_mode'] for a in accepted.values()))))
assert len(out)==9
(S/'validation/REAL_LONG_TOLERANCE_SENSITIVITY.json').write_text(json.dumps(dict(status='PASS',scope='Nine complete t=0 to 0.51 s real LONG prefix integrations; synthetic sensitivity in SYNTHETIC_SURFACE_QUALIFICATION.json',reference='SENS_T1_D15_mpi1',cases=out),indent=2)+'\n');print('REAL_LONG_TOLERANCE_SENSITIVITY_PASS',[(r['case'],r['candidate_count_changes'],r['max_position_difference_m']) for r in out])
