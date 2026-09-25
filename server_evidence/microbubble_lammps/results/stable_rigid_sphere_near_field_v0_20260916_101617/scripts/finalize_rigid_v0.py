"""Independent finalizer: recompute from CSV states/steps, never call the C++ library."""
from pathlib import Path
import json,sys,itertools,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'));sys.path.insert(0,str(R/'src'))
from validate_cases import evaluate_case,gap_convergence,table,xyz,qvals
from reference_rigid_sphere_lubrication import pair_matrix,coefficients,FlowGradientReference

def save(name,d): (R/name).parent.mkdir(exist_ok=True,parents=True);(R/name).write_text(json.dumps(d,indent=2)+'\n')
def compare(a,b,mapids=None):
 ta=table(R/'cases'/a/'TRAJECTORIES.csv');tb=table(R/'cases'/b/'TRAJECTORIES.csv');ma=mv=mo=mg=0
 for tag in np.unique(ta['particle_id']):
  aa=ta[ta['particle_id']==tag];bb=tb[tb['particle_id']==(mapids.get(int(tag),int(tag)) if mapids else tag)];assert len(aa)==len(bb);assert np.max(abs(aa['time_s']-bb['time_s']))<1e-12
  ma=max(ma,float(np.max(np.linalg.norm(xyz(aa)-xyz(bb),axis=1))));qa=qvals(aa);qb=qvals(bb);mv=max(mv,float(np.max(abs(qa[:,:3]-qb[:,:3]))));mo=max(mo,float(np.max(abs(qa[:,3:]-qb[:,3:]))));mg=max(mg,float(np.max(abs(aa['nearest_gap_m']-bb['nearest_gap_m']))))
 return {'max_position_m':ma,'max_velocity_m_s':mv,'max_omega_rad_s':mo,'max_gap_m':mg}

def main(full=False):
 cases=json.loads((R/'contracts/CASE_INDEX.json').read_text());checks={};summaries={}
 for c in cases:
  r=evaluate_case(c['name']);summaries[c['name']]=r;checks[c['name']]=r['status']=='PASS';print(c['name'],r['status'],flush=True)
 cg=gap_convergence();checks['Cgap']=cg['status']=='PASS'
 h=compare('G_three_particle','H_permuted',{1:3,2:1,3:2});h['status']='PASS' if h['max_position_m']<=1e-11 and h['max_velocity_m_s']<=1e-10 and h['max_omega_rad_s']<=1e-7 else 'FAIL';save('validation/PERMUTATION_AUDIT.json',h);checks['permutation']=h['status']=='PASS'
 m=compare('I_mpi1','I_mpi4');m['cross_rank_active_pair_observations']=json.loads((R/'cases/I_mpi4/RUN_STATE.json').read_text())['cross_rank_active_pair_observations'];m['status']='PASS' if m['max_position_m']<=1e-10 and m['max_velocity_m_s']<=1e-10 and m['max_omega_rad_s']<=1e-7 and m['cross_rank_active_pair_observations']>0 else 'FAIL';save('validation/MPI_AUDIT.json',m);checks['MPI']=m['status']=='PASS'
 oldlocal=Path('/home/lzy/projects/compre_output/passive_transport_v0/20260915_233516');oldremote=Path('/workspace/microbubble_lammps/results/passive_transport_v0_20260915_233516');old=oldlocal if oldlocal.exists() else oldremote
 # Exact old passive case2 on common integer steps. New omega does not feed back on isolated V.
 a=table(R/'cases/J_passive/TRAJECTORIES.csv');b=table(old/'cases/case2_c025/trajectory_rank0.csv');aa=np.array([a[a['step']==k][0] for k in b['step']],dtype=a.dtype);px=float(np.max(abs(xyz(aa)-xyz(b))));pv=float(np.max(abs(qvals(aa)[:,:3]-np.column_stack([b[k] for k in ['vx_m_s','vy_m_s','vz_m_s']]))));j={'status':'PASS' if px<=1e-12 and pv<=1e-12 else 'FAIL','old_case':'case2_c025','compared_steps':len(b),'old_final_step':int(b['step'][-1]),'new_final_step':int(a['step'][-1]),'max_position_diff_m':px,'max_velocity_diff_m_s':pv,'old_rotation':'absent','new_rotation':'Omega_inf, no center trajectory feedback outside cutoff'};save('validation/PASSIVE_REGRESSION.json',j);j['terminal_position_diff_m']=float(np.max(abs(xyz(a)[-1]-xyz(b)[-1])));j['terminal_time_diff_s']=float(abs(a['time_s'][-1]-b['time_s'][-1]));j['extra_terminal_roundoff_step']=j['new_final_step']-j['old_final_step'];save('validation/PASSIVE_REGRESSION.json',j);checks['passive_regression']=j['status']=='PASS' and j['terminal_position_diff_m']<=1e-12 and j['terminal_time_diff_s']<=1e-12
 # Nonempty pair outside the cutoff, directly matched to old recorded passive trajectories.
 jp=table(R/'cases/J2_passive_separated_pair/TRAJECTORIES.csv');jo=table(R/'provenance/J2_OLD_SELECTED_TRAJECTORIES.csv');px=pv=0
 for tag in [1,2]:
  aa=jp[jp['particle_id']==tag];bb=jo[jo['particle_id']==tag];matched=np.array([aa[aa['step']==k][0] for k in bb['step']],dtype=aa.dtype);px=max(px,float(np.max(abs(xyz(matched)-xyz(bb)))));pv=max(pv,float(np.max(abs(qvals(matched)[:,:3]-np.column_stack([bb[k] for k in ['vx_m_s','vy_m_s','vz_m_s']])))))
 j1=jp[jp['particle_id']==1];j2=jp[jp['particle_id']==2];assert len(j1)==len(j2)
 jr1=j1['diameter_um']*.5e-6;jr2=j2['diameter_um']*.5e-6;jg=np.linalg.norm(xyz(j1)-xyz(j2),axis=1)-jr1-jr2;cutoff=.2*jr1*jr2/(jr1+jr2)
 endpoint_position=endpoint_time=0
 for tag in [1,2]:
  aa=jp[jp['particle_id']==tag];bb=jo[jo['particle_id']==tag];endpoint_position=max(endpoint_position,float(np.max(abs(xyz(aa)[-1]-xyz(bb)[-1]))));endpoint_time=max(endpoint_time,float(abs(aa['time_s'][-1]-bb['time_s'][-1])))
 outside=bool(np.all(jg>=cutoff) and np.all(jp['active_pair_count']==0))
 jpcheck={'status':'PASS' if px<=1e-12 and pv<=1e-12 and outside and endpoint_position<=1e-12 and endpoint_time<=1e-12 else 'FAIL','original_ids':[1,2],'old_reference':'provenance/J2_OLD_SELECTED_TRAJECTORIES.csv','max_position_difference_m':px,'max_velocity_difference_m_s':pv,'all_pairs_outside_cutoff':outside,'independent_minimum_gap_m':float(jg.min()),'maximum_cutoff_m':float(cutoff.max()),'terminal_position_difference_m':endpoint_position,'terminal_time_difference_s':endpoint_time,'N':2};save('validation/NONEMPTY_PAIR_PASSIVE_REGRESSION.json',jpcheck);checks['nonempty_pair_passive_regression']=jpcheck['status']=='PASS'
 # All K gaps are recomputed offline from recorded coordinates. This is not a production neighbor search.
 k=table(R/'cases/K_real_8/TRAJECTORIES.csv');min_gap=np.inf;records=[]
 for step in np.unique(k['step']):
  rr=k[k['step']==step];x=xyz(rr);a=rr['diameter_um']*.5e-6;nearest=np.full(len(rr),np.inf);paircount=np.zeros(len(rr),int)
  for i,l in itertools.combinations(range(len(rr)),2):
   gap=np.linalg.norm(x[i]-x[l])-a[i]-a[l];nearest[i]=min(nearest[i],gap);nearest[l]=min(nearest[l],gap);min_gap=min(min_gap,gap)
   if gap<.2*a[i]*a[l]/(a[i]+a[l]):paircount[i]+=1;paircount[l]+=1
  for i in range(len(rr)):records.append([step,rr['time_s'][i],rr['particle_id'][i],nearest[i],paircount[i]])
 np.savetxt(R/'validation/CASE_K_NEAREST_PAIR_GAPS.csv',records,delimiter=',',header='step,time_s,particle_id,nearest_gap_m,active_pair_count',comments='',fmt='%.17g')
 state=json.loads((R/'cases/K_real_8/RUN_STATE.json').read_text());oldc=json.loads((old/'cases/case5_c025/CASE_CONTRACT.json').read_text());newc=json.loads((R/'cases/K_real_8/CASE_CONTRACT.json').read_text());identity=np.array_equal(np.array(newc['initial_positions_m']),np.array(oldc['initial_positions_m'])) and np.array_equal(2e6*np.array(newc['radii_m']),np.array(oldc['diameters_um'])) and newc['field_sha256']==oldc['field_sha256']
 # Compare floating identities using exact stored particle radii rather than a multiply round-trip.
 oldparticles=(old/'cases/case5_c025/particles.data').read_text().split('Atoms # sphere\n\n')[1].split('\n\nVelocities')[0].strip().splitlines();oldr=np.array([float(line.split()[2])*.5 for line in oldparticles]);identity=bool(np.array_equal(np.array(newc['initial_positions_m']),np.array(oldc['initial_positions_m'])) and np.array_equal(np.array(newc['radii_m']),oldr) and newc['field_sha256']==oldc['field_sha256'])
 kcheck={'status':'PASS' if identity and (state['time_s']>.008514 or state['reason']=='WALL_SAFETY_STOP') and min_gap>=-1e-12 else 'FAIL','same_population_positions_field':identity,'old_stop_time_s':.005675914432167097,'new_time_s':state['time_s'],'reason':state['reason'],'minimum_pair_gap_m':float(min_gap),'max_angular_speed_rad_s':float(np.max(k['angular_speed_rad_s'])),'wall_safety_stop':state['reason']=='WALL_SAFETY_STOP','exact_nearest_gap_source':'validation/CASE_K_NEAREST_PAIR_GAPS.csv','raw_neighbor_gap_scope':'Raw trajectory uses nearest pair within 8um LAMMPS list; 1e100 sentinel denotes no such neighbor; final VTP uses exact offline all-pair nearest gap.'};save('validation/CASE_K_REPLAY_AUDIT.json',kcheck);checks['old8']=kcheck['status']=='PASS'
 evb=table(R/'cases/B_tangential/PAIR_HYDRODYNAMIC_EVENTS.csv');initial_normal=float(abs(evb['normal_relative_speed_m_s'][0]));evg=table(R/'cases/G_three_particle/PAIR_HYDRODYNAMIC_EVENTS.csv');initial_pairs=int(np.count_nonzero(evg['step']==0));dd=table(R/'cases/D_pass_by/TRAJECTORIES.csv');motion={'initial_B_normal_speed_m_s':initial_normal,'G_initial_simultaneous_pairs':initial_pairs,'D_max_angular_speed_rad_s':float(np.max(dd['angular_speed_rad_s']))};motion['status']='PASS' if initial_normal<=1e-12 and initial_pairs>=2 and motion['D_max_angular_speed_rad_s']>0 else 'FAIL';save('validation/REQUIRED_MOTION_COVERAGE.json',motion);checks['required_motion_coverage']=motion['status']=='PASS'
 checks['contact_continuation']=json.loads((R/'cases/F_contact_cg04/RUN_STATE.json').read_text())['constraint_steps']>0
 # Pure prescribed kinematic modes supplement torque-driven transient C cases.
 mode_data=[];a=b=1e-6;hgap=2.5e-8;xij=np.array([-(a+b+hgap),0,0]);M,c=pair_matrix(a,b,xij)
 for mode in ['normal','tangent','C1_equal_transverse','C2_opposite_transverse']:
  q=np.zeros((2,6))
  if mode=='normal':q[0,0]=1e-3;q[1,0]=-1e-3
  if mode=='tangent':q[0,1]=1e-3;q[1,1]=-1e-3
  if mode=='C1_equal_transverse':q[:,5]=1000
  if mode=='C2_opposite_transverse':q[:,5]=[1000,-1000]
  f=(-M@q.reshape(-1)).reshape(2,6);mode_data.append({'mode':mode,'q_SI':q.tolist(),'force_torque_SI':f.tolist(),'positive_dissipation_W':float(q.reshape(-1)@M@q.reshape(-1))})
 save('validation/PURE_KINEMATIC_MODES_REFERENCE.json',{'role':'independent analytical mode diagnostics; C++ full matrix agreement covered by coefficient/matrix audits','modes':mode_data,'C3':'BLOCKED_TWIST_NOT_IMPLEMENTED'})
 if full:
  import vtk
  from vtk.util.numpy_support import numpy_to_vtk,vtk_to_numpy
  reader=vtk.vtkSTLReader();reader.SetFileName(str(R/'provenance/closed_geometry_m.stl'));reader.Update();dist=vtk.vtkImplicitPolyDataDistance();dist.SetInput(reader.GetOutput());minimum=np.inf;count=0;bad=0;maxdiff=0
  def distance(x):return abs(dist.EvaluateFunction(x))
  def segment_safe(a,b,min_allowed,depth=0):
   da=distance(a);db=distance(b);length=np.linalg.norm(b-a)
   if min(da,db)<min_allowed-1e-14:return False
   if min(da,db)-length/2>=min_allowed:return True
   if depth>=20:return False
   mid=.5*(a+b);return segment_safe(a,mid,min_allowed,depth+1) and segment_safe(mid,b,min_allowed,depth+1)
  for name in ['J_passive','K_real_8']:
   rr=table(R/'cases'/name/'TRAJECTORIES.csv');ss=table(R/'cases'/name/'INTEGRATION_STAGES.csv');contract=json.loads((R/'cases'/name/'CASE_CONTRACT.json').read_text());radii=dict(zip(contract['ids'],contract['radii_m']));margin=5.996975424319603e-7
   for row in rr:
    x=np.array([row[v] for v in ['x_m','y_m','z_m']]);gap=distance(x)-radii[int(row['particle_id'])];minimum=min(minimum,gap);maxdiff=max(maxdiff,abs(gap-row['surface_wall_gap_m']));count+=1;bad+=gap<margin-1e-14
   for tag in radii:
    t=rr[rr['particle_id']==tag];sx=ss[ss['particle_id']==tag];xs=xyz(sx);threshold=radii[tag]+margin
    for z in xs:count+=1;bad+=distance(z)<threshold-1e-14
    for aa,bb in zip(xyz(t)[:-1],xyz(t)[1:]):bad+=not segment_safe(aa,bb,threshold)
   points=vtk.vtkPoints();cloudpoints=np.vstack([xyz(rr),xyz(ss)]);points.SetData(numpy_to_vtk(cloudpoints,deep=True));cloud=vtk.vtkPolyData();cloud.SetPoints(points);enclosed=vtk.vtkSelectEnclosedPoints();enclosed.SetInputData(cloud);enclosed.SetSurfaceData(reader.GetOutput());enclosed.SetTolerance(1e-9);enclosed.Update();bad+=int(np.count_nonzero(vtk_to_numpy(enclosed.GetOutput().GetPointData().GetArray('SelectedPoints'))==0))
  wall={'status':'PASS' if bad==0 and maxdiff<=1e-12 else 'FAIL','independent_engine':'VTK exact STL point/triangle distance and recursive segment certificate','points':count,'violations':int(bad),'minimum_surface_wall_gap_m':float(minimum),'max_cpp_vtk_distance_error_m':float(maxdiff)};save('validation/LOCAL_INDEPENDENT_WALL_AUDIT.json',wall);checks['independent_wall']=wall['status']=='PASS'
 result={'status':'PASS_WITH_TWIST_PENDING' if all(checks.values()) else 'FAIL','checks':checks,'cases':summaries,'twist':'PENDING_BLOCKER','full_rotational_lubrication':'NOT_PASS_TWIST_PENDING','C3':'BLOCKED_TWIST_NOT_IMPLEMENTED','selected_C_gap':cg['selected_C_gap'],'MPI':m,'permutation':h,'passive_regression':j,'Case_K':kcheck,'all_validation_recomputed_from_raw_CSV':True,'Cplusplus_loaded':False,'local_independent_wall_required_for_final':not full,'human_visual_review':'PENDING'}
 save('LOCAL_INDEPENDENT_FINALIZER.json' if full else 'NUMERICAL_FINALIZER.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['cases']},indent=2));return 0 if all(checks.values()) else 2
if __name__=='__main__':raise SystemExit(main('--full-local' in sys.argv))
