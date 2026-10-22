"""Independent completed-run audit and physical-time visualization data.
Does not modify solver states, sources, flow, or terminal classifications.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,sys,collections
import numpy as np
from scipy.spatial import cKDTree
ap=argparse.ArgumentParser();ap.add_argument('--run-name',required=True);ap.add_argument('--run-collection',choices=['runs','validation_runs'],default='runs');ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];run=root/args.run_collection/args.run_name
sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.provenance import write_json,digest,canonical
from particle_3d.wall_gap import wall_gap
from particle_3d.particle_shapes import Sphere
manifest=json.loads((run/'MANIFEST.json').read_text());summary=json.loads((run/'SUMMARY.json').read_text())
assert summary['status']=='COMPLETED_WITH_EXPLICIT_TERMINALS'
out=root/'analysis'/args.run_name;out.mkdir(parents=True,exist_ok=True)
(out/'tracks').mkdir(exist_ok=True)
def ledger(name):
 with gzip.open(run/(name+'.jsonl.gz'),'rt') as f:return [json.loads(l) for l in f]
source=ledger('source');terminals=ledger('terminals');byterminal={r['particle_id']:r for r in terminals}
accepted=[r for r in source if r['status']=='ACCEPTED_BIRTH'];bybirth={r['particle_id']:r for r in accepted}
assert len(bybirth)==len(accepted)==summary['accepted_birth_count']
assert len(byterminal)==len(terminals)==len(accepted)
assert len(source)==summary['source_event_count']
assert all(r['diameter_draw_count']==r['position_draw_count']==r['admission_check_count']==1 for r in source)
assert all(r['particle_id'] is None for r in source if r['status']!='ACCEPTED_BIRTH')
assert all(r['birth_time_s']<=summary['infusion_stop_time_s'] for r in source)
assert all(source[k]['birth_time_s']<source[k+1]['birth_time_s'] for k in range(len(source)-1))
ledger_hashes={}
for name in ['source','terminals','interaction_audit']:
 h=hashlib.sha256()
 with gzip.open(run/(name+'.jsonl.gz'),'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 ledger_hashes[name]=h.hexdigest();assert h.hexdigest()==manifest['ledger_scientific_hashes'][name]
parts={k:[] for k in ['t','id','x','v','gap','tetra']};states=[];nearfield_min=float('inf');pairchecks=0;next_pair_time=0.;finite=True;max_age_error=0.;max_kinematic_error=0.
for chunk in manifest['chunks']:
 path=run/chunk['path'];assert digest(path)==chunk['file_sha256']
 with np.load(path) as d:
  h=hashlib.sha256()
  for k in sorted(d.files):
   a=np.ascontiguousarray(d[k]);h.update(canonical(dict(name=k,shape=a.shape,dtype=str(a.dtype))));h.update(a.tobytes())
   if np.issubdtype(a.dtype,np.floating):finite=finite and bool(np.isfinite(a).all())
  assert h.hexdigest()==chunk['array_sha256']
  times=np.repeat(d['state_time_s'],np.diff(d['state_offsets']))
  max_age_error=max(max_age_error,float(np.max(abs(d['age']-(times-d['birth_time'])))))
  for key,value in [('t',times),('id',d['particle_id']),('x',d['x']),('v',d['v']),('gap',d['wall_gap']),('tetra',d['tetra_id'])]:parts[key].append(value)
  for j,t in enumerate(d['state_time_s']):
   a,b=d['state_offsets'][j:j+2];states.append((float(t),int(b-a)))
   if t<next_pair_time or b-a<2:continue
   next_pair_time=t+.05;x=d['x'][a:b];r=d['radius'][a:b]
   pairs=cKDTree(x).query_pairs(2*r.max()+2e-9,output_type='ndarray');pairchecks+=1
   if len(pairs):
    gaps=np.linalg.norm(x[pairs[:,0]]-x[pairs[:,1]],axis=1)-r[pairs[:,0]]-r[pairs[:,1]]
    nearfield_min=min(nearfield_min,float(gaps.min()))
 print('AUDITED_CHUNK',chunk['index'],flush=True)
assert finite and max_age_error<1e-12
p={k:np.concatenate(v) for k,v in parts.items()};del parts
order=np.lexsort((p['t'],p['id']));p={k:v[order] for k,v in p.items()};del order
boundaries=np.r_[0,np.flatnonzero(np.diff(p['id']))+1,len(p['id'])]
env=load_environment(json.loads((root/'data/environment.json').read_text()))
rows=[];exitchecks=[];failures=[];wall_samples=[]
frames=576;frame_time=np.linspace(0,summary['physical_time_s'],frames)
frame_x=np.full((frames,len(accepted),3),np.nan);trackpoints=[];trackoffsets=[0]
for n,(a,b) in enumerate(zip(boundaries[:-1],boundaries[1:])):
 pid=int(p['id'][a]);birth=bybirth[pid];terminal=byterminal[pid]
 t=p['t'][a:b];keep=np.r_[np.diff(t)>0,True];t=t[keep];x=p['x'][a:b][keep];v=p['v'][a:b][keep];gap=p['gap'][a:b][keep];tetra=p['tetra'][a:b][keep]
 assert np.all(np.diff(t)>0) and abs(t[0]-birth['birth_time_s'])<1e-12
 assert abs(t[-1]-terminal['physical_time_s'])<1e-10,(pid,t[-1],terminal['physical_time_s'])
 error=np.linalg.norm(np.diff(x,axis=0)-np.diff(t)[:,None]*v[1:],axis=1)
 if len(error):max_kinematic_error=max(max_kinematic_error,float(error.max()))
 row={k:terminal.get(k) for k in ['particle_id','source_event_id','realization_id','status','terminal_reason','outlet','birth_time_s','exit_time_s','residence_time_s','physical_time_s','trigger']}
 row.update(diameter_um=birth['diameter_um'],radius_m=birth['radius_m'],samples=len(t),minimum_recorded_wall_gap_m=float(gap.min()),terminal_x_m=x[-1,0],terminal_y_m=x[-1,1],terminal_z_m=x[-1,2])
 rows.append(row)
 if terminal['status']=='EXITED':
  hit=env.classifier.first_event(x[-2],x[-1])
  exitchecks.append(dict(particle_id=pid,recorded_outlet=terminal['outlet'],independent_outlet=None if hit is None else hit.role,
      segment_fraction=None if hit is None else hit.segment_fraction,PASS=hit is not None and hit.role==terminal['outlet']))
 if terminal['status']=='SOLVER_FAILURE':
  f=env.field.sample(x[-1]);ev=terminal.get('support_evidence',[]);g=wall_gap(Sphere(x[-1],birth['radius_m']),env.wall)
  record=dict(row,actual_FEM_velocity_m_s=f.velocity_m_s.tolist(),actual_FEM_speed_m_s=float(np.linalg.norm(f.velocity_m_s)),
      support_history_count=len(ev),CPU_GPU_velocity_max_error_m_s=max((float(np.max(abs(np.array(e['velocity'])-e['independently_recomputed_velocity_m_s']))) for e in ev),default=None),
      terminal_tetra_id=int(tetra[-1]),triangle_id=int(g.wall_triangle_id),wall_gap_m=float(g.gap_m),
      contact_ranks=sorted(set(e['solver'].get('contact_redundancy',{}).get('rank_after',0) for e in ev)))
  if tetra[-1]>=0:record['tetra_vertex_speeds_m_s']=np.linalg.norm(env.field.velocity_nodes_m_s[env.field.tetra[tetra[-1]]],axis=1).tolist()
  failures.append(record)
 # Stratified and smallest recorded gaps, finite-triangle independent check.
 for k in sorted(set([0,len(t)-1,int(np.argmin(gap))])):
  wall_samples.append((pid,float(t[k]),x[k],birth['radius_m'],float(gap[k]),k==0))
 # Raw birth stores the admission lower-bound clearance and tetra_id=-1.
 # This derived export fills exact birth metadata explicitly; raw chunks stay intact.
 gap=gap.copy();tetra=tetra.copy();gap[0]=wall_gap(Sphere(x[0],birth['radius_m']),env.wall).gap_m;tetra[0]=env.field.locate(x[0])[0]
 np.savez_compressed(out/'tracks'/f'mb_{pid:06d}.npz',physical_time_s=t,age_s=t-birth['birth_time_s'],x_m=x,velocity_m_s=v,wall_gap_m=gap,tetra_id=tetra)
 live=(frame_time>=t[0])&(frame_time<t[-1])
 for j in range(3):frame_x[live,n,j]=np.interp(frame_time[live],t,x[:,j])
 # Exact vertices, display decimation only; raw accepted endpoints remain above.
 stride=max(1,int(np.ceil((len(t)-1)/1499)));ix=np.unique(np.r_[np.arange(0,len(t),stride),len(t)-1]);trackpoints.append(x[ix]);trackoffsets.append(trackoffsets[-1]+len(ix))
print('EXIT_PARITY',sum(r['PASS'] for r in exitchecks),len(exitchecks),flush=True)
wall_rows=[]
for pid,t,x,r,stored,is_birth in wall_samples:
 g=wall_gap(Sphere(x,r),env.wall);wall_rows.append(dict(particle_id=pid,time_s=t,recorded_gap_m=stored,independent_gap_m=g.gap_m,error_m=abs(stored-g.gap_m) if not is_birth else 0.,kind='BIRTH_ADMISSION_LOWER_BOUND' if is_birth else 'ACCEPTED_EXACT_GAP',lower_bound_violation_m=max(0.,stored-g.gap_m) if is_birth else 0.,triangle_id=int(g.wall_triangle_id)))
wallerror=max(r['error_m'] for r in wall_rows);birth_bound_violation=max(r['lower_bound_violation_m'] for r in wall_rows);minwall=min(r['independent_gap_m'] for r in wall_rows)
np.savez_compressed(out/'visualization_data.npz',frame_time_s=frame_time,frame_x_m=frame_x,
 particle_id=np.array([r['particle_id'] for r in rows]),track_points_m=np.concatenate(trackpoints),track_offsets=np.array(trackoffsets))
def csvwrite(name,items):
 if not items:return
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
csvwrite('trajectory_catalog.csv',rows);csvwrite('independent_outlet_checks.csv',exitchecks);csvwrite('independent_wall_checks.csv',wall_rows)
# Counts at all accepted snapshot times; cross-check full source/outcome accounting.
birth_times=np.sort([r['birth_time_s'] for r in accepted]);term_times=np.sort([r['physical_time_s'] for r in terminals]);exittimes={o:np.sort([r['physical_time_s'] for r in terminals if r['outlet']==o]) for o in summary['outlets']}
timeline=[]
for t in np.unique(np.r_[0,frame_time,birth_times,term_times]):
 born=int(np.searchsorted(birth_times,t,side='right'));done=int(np.searchsorted(term_times,t,side='right'))
 timeline.append(dict(physical_time_s=float(t),accepted_births=born,active_count=born-done,terminated=done,
     **{o:int(np.searchsorted(times,t,side='right')) for o,times in exittimes.items()}))
csvwrite('population_timeline.csv',timeline)
interaction_rows=[]
with gzip.open(run/'interaction_audit.jsonl.gz','rt') as f:
 for line in f:
  r=json.loads(line)
  for pair in r.get('pair_ids',[]):interaction_rows.append(dict(physical_time_s=r['physical_time_s'],particle_i=pair[0],particle_j=pair[1],contact=pair in r.get('pair_contact_ids',[])))
csvwrite('pair_interactions.csv',interaction_rows)
write_json(out/'failure_diagnostics.json',failures)
coverage={k:0 for k in summary['outlets']};firstcoverage=None
for r in sorted(terminals,key=lambda r:r['physical_time_s']):
 if r['status']=='EXITED':coverage[r['outlet']]+=1
 if firstcoverage is None and min(coverage.values())>=20:firstcoverage=r['physical_time_s']
test_override=manifest['identity']['population_contract'].get('fixed_infusion_test_override_s')
if test_override is None:assert firstcoverage==summary['infusion_stop_time_s']
result=dict(status='AUDIT_COMPLETED_WITH_EXPLICIT_LIMITATIONS',run_name=args.run_name,flow_sha256=env.flow_sha256,
 production_identity=manifest['identity'],postprocessing_script_sha256=digest(Path(__file__)),
 source_and_terminal_accounting_pass=True,all_file_and_scientific_hashes_pass=True,finite_arrays=finite,
 maximum_age_identity_error_s=max_age_error,maximum_held_velocity_path_error_m=max_kinematic_error,
 recorded_particle_states=len(p['id']),independent_outlet_checks=len(exitchecks),independent_outlet_pass=all(r['PASS'] for r in exitchecks),
 independent_wall_checks=len(wall_rows),birth_admission_lower_bound_max_violation_m=birth_bound_violation,independent_wall_max_error_m=wallerror,independent_wall_min_gap_m=minwall,
 wall_roundoff_m=env.wall.roundoff_m,pair_sampled_states=pairchecks,
 sampled_near_contact_minimum_raw_gap_m=None if not np.isfinite(nearfield_min) else nearfield_min,
 pair_check_scope='KD tree at >=0.05 s intervals, all pairs within 2*r_max+2 nm; not every-step independent proof',
 coverage_stop_time_s=firstcoverage,source_stop_rule_reproduced=firstcoverage==summary['infusion_stop_time_s'],fixed_infusion_test_override_s=test_override,summary=summary,
 unique_nearfield_pairs=len(set((r['particle_i'],r['particle_j']) for r in interaction_rows)),
 unique_contact_pairs=len(set((r['particle_i'],r['particle_j']) for r in interaction_rows if r['contact'])),
 failure_count=len(failures),scientific_integrity_pass=summary['scientific_integrity_pass'],
 animation_time='ABSOLUTE_PHYSICAL_TIME; held-velocity interpolation; no age alignment',
 diameter_um=dict(min=min(r['diameter_um'] for r in rows),max=max(r['diameter_um'] for r in rows)))
oldrun=root/'runs/realization_0000_final'
if args.run_name=='realization_0000_contactfix' and (oldrun/'MANIFEST.json').exists():
 oldmanifest=json.loads((oldrun/'MANIFEST.json').read_text())
 comparison=dict(scope='SOURCE_LEDGER_AND_PREVIOUSLY_FAILED_ID_ONLY; OLD_RUN_NOT_VALID_PRODUCTION',
  source_scientific_hash_identical=ledger_hashes['source']==oldmanifest['ledger_scientific_hashes']['source'])
 with gzip.open(oldrun/'terminals.jsonl.gz','rt') as f:
  oldterminal=next((json.loads(l) for l in f if json.loads(l)['particle_id']==1131),None)
 if oldterminal is not None and 1131 in byterminal:
  track=np.load(out/'tracks/mb_001131.npz');ot=oldterminal['physical_time_s'];t=track['physical_time_s'];x=track['x_m'];oldx=np.array(oldterminal['support_evidence'][-1]['position_m'])
  later=x[t>ot];comparison.update(particle_id=1131,previous_terminal_time_s=ot,
   revised_terminal={k:v for k,v in byterminal[1131].items() if k!='support_evidence'},later_recorded_states=len(later),
   maximum_later_displacement_from_previous_stop_m=float(np.linalg.norm(later-oldx,axis=1).max()) if len(later) else 0.)
 write_json(out/'contact_fix_comparison.json',comparison);result['contact_fix_comparison']=comparison
write_json(out/'INDEPENDENT_AUDIT.json',result);print(json.dumps(result,indent=2),flush=True)
