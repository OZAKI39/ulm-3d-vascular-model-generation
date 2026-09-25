from pathlib import Path
import json,csv,sys
import numpy as np,yaml,h5py,vtk
from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray,vtk_to_numpy
S=Path(__file__).resolve().parents[1];cfg=yaml.safe_load((S/'configs/workflow_v0.yaml').read_text());G=dict(np.load(S/'geometry/GEOMETRY_ARRAYS.npz'));manifest=json.loads((S/'geometry/BOUNDARY_MANIFEST.json').read_text());inlet=next(p for p in manifest['ports'] if p['role']=='INLET_PORT');origin=np.array(inlet['center_m']);normal=np.array(inlet['outward_unit_normal']);axis=np.eye(3)[np.argmin(abs(normal))];e1=np.cross(normal,axis);e1/=np.linalg.norm(e1);e2=np.cross(normal,e1)
sys.path.insert(0,str(S/'inputs'));from sonovue_sampler import SonoVueDistribution
law=SonoVueDistribution(S/cfg['bubbles']['sonovue_distribution']);seed=cfg['bubbles']['random_seed']
def poly(faces):
 m=vtk.vtkPolyData();p=vtk.vtkPoints();p.SetData(numpy_to_vtk(G['points'],deep=True));m.SetPoints(p);c=vtk.vtkCellArray();c.SetCells(len(faces),numpy_to_vtkIdTypeArray(np.c_[np.full(len(faces),3),faces].astype(np.int64).ravel(),deep=True));m.SetPolys(c);return m
closed=poly(G['faces']);walls=poly(G['faces'][np.isin(G['classes'],[0,1,2])]);locator=vtk.vtkStaticCellLocator();locator.SetDataSet(walls);locator.BuildLocator();inside=vtk.vtkImplicitPolyDataDistance();inside.SetInput(closed)
with h5py.File(S/cfg['paths']['frozen_flow'],'r') as f:
 dims=f['dims'][:].astype(int);floworigin=f['origin_m'][:];dx=float(f['dx_m'][()]);indices=f['linear_index'][:];fluid=f['fluid_linear_index'][:];vel=f['Velocity_m_s'][:];attrs={k:(v.decode() if isinstance(v,bytes) else str(v)) for k,v in f.attrs.items()}
assert attrs['coordinate_units']=='m' and attrs['velocity_units']=='m/s'
def flowquery(p):
 q=(p-floworigin)/dx
 if np.any(q<0) or np.any(q>=dims-1):return 'OUTSIDE',np.zeros(3)
 b=np.floor(q).astype(int);frac=q-b;result=np.zeros(3)
 for z in [0,1]:
  for y in [0,1]:
   for x in [0,1]:
    ijk=b+[x,y,z];index=int(ijk[0]+dims[0]*(ijk[1]+dims[1]*ijk[2]));ff=np.searchsorted(fluid,index);ii=np.searchsorted(indices,index)
    if ff>=len(fluid) or fluid[ff]!=index:return 'SOLID',np.zeros(3)
    if ii>=len(indices) or indices[ii]!=index:return 'MISSING',np.zeros(3)
    w=np.prod(np.where(np.array([x,y,z]),frac,1-frac));result+=w*vel[ii]
 return 'VALID',result
def wallquery(p):
 q=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.);locator.FindClosestPoint(p,q,cell,sub,d2);return np.sqrt(float(d2)),q,int(cell)
def section(offset):
 point=origin-offset*normal;plane=vtk.vtkPlane();plane.SetOrigin(point);plane.SetNormal(normal);cut=vtk.vtkCutter();cut.SetCutFunction(plane);cut.SetInputData(closed);cut.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION);cut.Update();strip=vtk.vtkStripper();strip.SetInputConnection(cut.GetOutputPort());strip.JoinContiguousSegmentsOn();strip.Update();m=strip.GetOutput();p=vtk_to_numpy(m.GetPoints().GetData());cells=vtk_to_numpy(m.GetLines().GetData());loops=[];j=0
 while j<len(cells):
  count=int(cells[j]);ids=cells[j+1:j+count+1];loop=p[ids];j+=count+1
  if len(loop)>=4 and np.linalg.norm(loop[0]-loop[-1])<1e-10:loops.append(loop[:-1])
 assert loops,'STOP_INJECTION_SECTION_NOT_CLOSED';loop=min(loops,key=lambda p:np.linalg.norm(p.mean(axis=0)-point));xy=np.c_[(loop-point)@e1,(loop-point)@e2];area=.5*abs(np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1)));widths=np.ptp(xy,axis=0);return point,loop,area,widths
capids=np.flatnonzero(G['classes']==3);captri=G['points'][G['faces'][capids]];probs=G['area'][capids]/G['area'][capids].sum();summary=[];allpositions=[];offset_scan=[]
for offset in np.array([.5,1,2,3,4,5])*1e-6:
 point,loop,area,widths=section(offset);offset_scan.append({'offset_m':float(offset),'section_area_m2':area,'widths_m':widths.tolist(),'upper_bound_inscribed_radius_m':float(min(widths)/2)})
for count in cfg['bubbles']['counts']:
 name={1:'W1',8:'W2',16:'W3'}[count];case=S/'cases'/name;case.mkdir(exist_ok=True);u=np.random.default_rng(seed).random(count);diam=law.inverse_cdf(u)*1e-6;radii=diam/2;margin=dx;offset=max(4*dx,float(radii.max())+margin);plane,loop,section_area,widths=section(offset);diskarea=float(np.pi*np.sum(radii*radii));radiusbound=float(min(widths)/2)
 rows=[{'particle_id':i+1,'diameter_m':float(diam[i]),'radius_m':float(radii[i]),'sample_uniform_0_1':float(u[i]),'random_seed':seed} for i in range(count)]
 with (case/'POPULATION_DRAWN.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 feasible_area=bool(diskarea<=section_area);feasible_size=bool(np.all(radii+margin<=radiusbound));accepted=[];reason='';attempts=[]
 # A failed necessary area/width condition proves that position retries cannot succeed.
 if not feasible_area:reason='COPLANAR_DISK_AREA_EXCEEDS_FLUID_SECTION'
 elif not feasible_size:reason='BUBBLE_RADIUS_EXCEEDS_SECTION_WIDTH_BOUND_WITH_CLEARANCE'
 else:
  for i,a in enumerate(radii):
   rng=np.random.default_rng(seed+1000*count+i);done=False
   for attempt in range(cfg['injection']['maximum_position_attempts_per_bubble']):
    tid=rng.choice(len(captri),p=probs);r1,r2=rng.random(2);s=np.sqrt(r1);weights=np.array([1-s,s*(1-r2),s*r2]);point=weights@captri[tid]-offset*normal;distance,nearest,wallcell=wallquery(point);status,velocity=flowquery(point)
    if inside.EvaluateFunction(point)>=0 or distance-a<margin or status!='VALID':continue
    if any(np.linalg.norm(point-v['position'])<a+v['radius']+margin for v in accepted):continue
    accepted.append({'id':i+1,'position':point,'radius':a,'gap':distance-a,'velocity':velocity,'sample_uniform':u[i],'attempts':attempt+1});done=True;attempts.append(attempt+1);break
   if not done:reason='POSITION_REJECTION_BUDGET_EXHAUSTED';attempts.append(cfg['injection']['maximum_position_attempts_per_bubble']);break
 result={'case':name,'count':count,'seed':seed,'status':'PASS' if len(accepted)==count else 'INITIALIZATION_FAILED','offset_m':offset,'offset_rule':'max(4*dx,maximum_radius+dx)','initial_clearance_m':margin,'initial_clearance_rule':'one frozen field voxel','plane_center_m':plane.tolist(),'plane_normal':normal.tolist(),'section_area_m2':section_area,'bubble_disks_area_m2':diskarea,'disk_area_over_available':diskarea/section_area,'section_widths_m':widths.tolist(),'inscribed_radius_upper_bound_m':radiusbound,'maximum_bubble_radius_m':float(radii.max()),'necessary_coplanar_area_condition':feasible_area,'necessary_individual_size_with_margin_condition':feasible_size,'reason':reason,'position_retry_budget':cfg['injection']['maximum_position_attempts_per_bubble'],'actual_attempts_per_placed_bubble':attempts,'geometric_infeasibility_short_circuits_futile_retries':not feasible_area or not feasible_size,'radii_never_resampled_or_shrunk':True,'accepted_count':len(accepted)};summary.append(result)
 np.savez_compressed(case/'INJECTION_SECTION.npz',points=loop,plane_center=plane,normal=normal,axis1=e1,axis2=e2)
 if len(accepted)==count:
  out=[]
  for v in accepted:
   i=v['id']-1;out.append({**rows[i],'initial_x':v['position'][0],'initial_y':v['position'][1],'initial_z':v['position'][2],'initial_gap_m':v['gap'],'flow_vx_m_s':v['velocity'][0],'flow_vy_m_s':v['velocity'][1],'flow_vz_m_s':v['velocity'][2]});allpositions.append(v['position'].tolist())
  with (case/'BUBBLES_INITIAL.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 (case/'INITIALIZATION_STATUS.json').write_text(json.dumps(result,indent=2)+'\n')
report={'status':'PASS' if all(r['status']=='PASS' for r in summary) else 'FAIL','cases':summary,'offset_sensitivity':offset_scan,'packing_proof':'All sphere centres on one injection plane imply disjoint radius-a disks contained within its fluid polygon. Therefore sum(pi*a_i^2)<=polygon area is necessary. A disk also requires width>=2*(a+clearance) in every tangent direction. Failed necessary inequalities cannot be repaired by more position retries.','flow_sampling':'Preflight independent Python calculation equivalent to frozen C++ eight-corner trilinear sampler; Native donor execution is NOT_RUN after initialization STOP; this check is not a runtime validation.','FLOW_PHYSICS_STATUS':'ENGINEERING_TRANSIENT_FIELD_ONLY'};(S/'validation/INITIALIZATION_PREFLIGHT.json').write_text(json.dumps(report,indent=2)+'\n')
if report['status']=='FAIL':
 stop={'STOP_REASON':'STOP_INLET_BATCH_INITIALIZATION_INFEASIBLE','EXPECTED':'W1=1, W2=8 and W3=16 fixed SonoVue radii, simultaneously on one inward-offset inlet plane, without wall/pair overlap and with successful frozen-flow query.','ACTUAL':summary,'IMPACT':'W1 initialization can be checked, but W2/W3 cannot initialize under the frozen single-plane batch configuration. No RK2/MPI workflow cases have been started.','NEXT_OPTIONS':['User authorizes a 3D inlet buffer/slab batch placement policy; individually oversized bubbles still require an explicit admission policy or a wider inlet.','User chooses a wider real inlet/geometry compatible with the intended size range; establish matching frozen flow rather than scaling existing geometry.','User explicitly limits the first stage to W1 or authorizes a disclosed geometry-conditioned size population. No implicit filtering or seed selection is applied.'],'USER_RULE':'Request section57: stop when inlet injection plane cannot be sampled reliably; do not hack completion.'};(S/'validation/STOP.json').write_text(json.dumps(stop,indent=2)+'\n')
print(json.dumps(report,indent=2))
