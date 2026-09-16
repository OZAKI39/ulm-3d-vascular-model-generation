from pathlib import Path
import sys,json,hashlib,math,os,itertools
import numpy as np,vtk
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from reference_flow_sampler import ReferenceField
from reference_drag import tau

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def csv(path,a,header):np.savetxt(path,a,delimiter=',',header=header,comments='',fmt='%.17g')
selected=json.loads((P/'provenance/POPULATION_IDENTITY.json').read_text());dt=json.loads((P/'contracts/PARTICLE_TIMESTEP_CONTRACT.json').read_text())['dt_particle_s'];syn=json.loads((P/'contracts/SYNTHETIC_FIELDS_CONTRACT.json').read_text());real=ReferenceField(P/'fields/FROZEN_FLOW_FIELD_V0.h5')
queries=P/'validation/queries';queries.mkdir(exist_ok=True)
rng=np.random.default_rng(2026091801);positions=np.array(syn['origin_m'])+rng.uniform(.01,19.99,(10000,3))*syn['dx_m']
for name in ['uniform','linear']:csv(queries/(name+'.csv'),np.column_stack((np.arange(len(positions)),positions)),'id,x,y,z')
interior=real.interior_lower_ids();deep=real.interior_lower_ids(pad=1);assert len(deep)>=1000 and len(interior)>=10000
rng=np.random.default_rng(2026091802);nodeids=rng.choice(deep,1000,replace=False);csv(queries/'real_nodes.csv',np.column_stack((np.arange(1000),real.xyz(nodeids))),'id,x,y,z');np.save(queries/'real_node_ids.npy',nodeids)
cellids=rng.choice(interior,10000,replace=False);points=real.xyz(cellids)+rng.uniform(.05,.95,(10000,3))*real.dx;csv(queries/'real_interior.csv',np.column_stack((np.arange(10000),points)),'id,x,y,z')
invalid=np.array([[100e-6,0,0],[0,-100e-6,0],[0,0,100e-6],[np.nan,0,0]])
csv(queries/'invalid_domain.csv',np.column_stack((np.arange(4),invalid)),'id,x,y,z')
csv(queries/'invalid_missing.csv',np.array([[0,1e-6,1e-6,1e-6]]),'id,x,y,z')
# The minimum corner of the real grid is documented external solid; validate reference status.
solid=real.origin+real.dx*.1;assert real.sample(solid)[1][0]==2;csv(queries/'invalid_solid.csv',np.r_[0,solid][None,:],'id,x,y,z')
rng=np.random.default_rng(2026091803);d=rng.uniform(.75e-6,5.25e-6,10000);uf=rng.uniform(-.01,.01,(10000,3));vp=rng.uniform(-.01,.01,(10000,3));vp[:100]=uf[:100]
csv(queries/'force.csv',np.column_stack((np.arange(10000),d,uf,vp)),'id,d,ufx,ufy,ufz,vx,vy,vz')
# Diameter already selected. Deterministic safe-point search cannot redraw or resize it.
reader=vtk.vtkSTLReader();reader.SetFileName(str(P/'provenance/closed_geometry_m.stl'));reader.Update();mesh=reader.GetOutput();distance=vtk.vtkImplicitPolyDataDistance();distance.SetInput(mesh)
rng=np.random.default_rng(2026091900);candidates=rng.choice(interior,min(4096,len(interior)),replace=False);centers=real.xyz(candidates)+.5*real.dx
distances=np.array([abs(distance.EvaluateFunction(x)) for x in centers]);order=np.argsort(-distances,kind='stable');chosen=None
for k in order:
 x=centers[k];poly=vtk.vtkPolyData();pts=vtk.vtkPoints();pts.InsertNextPoint(*x);poly.SetPoints(pts);inside=vtk.vtkSelectEnclosedPoints();inside.SetInputData(poly);inside.SetSurfaceData(mesh);inside.SetTolerance(1e-7);inside.Update()
 if inside.IsInside(0):chosen=int(k);break
assert chosen is not None,'BLOCKED_NO_SAFE_INTERIOR_PATH'
spawn=centers[chosen];wall=float(distances[chosen]);diam=selected['C']['diameters_um'][0]*1e-6;available=wall-diam/2-3*real.dx
if available<=0:raise RuntimeError('PARTIAL_BLOCKED_GEOMETRY_SAFETY')
max_u=float(np.max(np.linalg.norm(real.u,axis=1)));horizon=min(5*float(tau(diam)),available/(4*max_u));csteps=int(math.floor(horizon/dt));assert csteps>0
# All force components are convex interpolation of frozen node vectors. |v| <= max_u
# for this zero-v0, positive-damping stable test; bound checked against output every step.
safety={'status':'PASS','selection_rule':'closest sampled diameter to target d50 fixed before fit; select largest exact STL distance among 4096 deterministic all-fluid cell centers that are enclosed','candidate_seed':2026091900,'candidate_count':len(candidates),'spawn_m':spawn.tolist(),'selected_lower_cell_index':int(candidates[chosen]),'diameter_m':diam,'sphere_radius_m':diam/2,'center_wall_distance_m':wall,'required_center_wall_distance_m':diam/2+3*real.dx,'available_center_displacement_m':available,'all_8_corners_fluid':True,'closed_surface_enclosed':True,'max_field_speed_m_s':max_u,'duration_rule':'min(5*tau_p, available_displacement/(4*max_field_speed)); floor to integer dt','steps':csteps,'physical_time_s':csteps*dt,'guaranteed_distance_rule':'distance_to_closed_surface(x) >= distance_at_spawn - norm(x-spawn), checked every step; no wall force','geometry_sha256':sha(P/'provenance/closed_geometry_m.stl'),'wall_distance_reference':'VTK implicit polydata exact triangle distance; enclosure via vtkSelectEnclosedPoints'}
put(P/'contracts/CASE_C_SAFE_INTERIOR_CONTRACT.json',safety)
np.savez_compressed(P/'validation/SAFE_INTERIOR_SEARCH.npz',candidate_lower_ids=candidates,positions_m=centers,center_wall_distances_m=distances)
allcases=[]
def add(name,pop,sel,kind='uniform',ranks=1,steps=None,engine='cpu'):
 ds=np.asarray(selected[pop]['diameters_um'])[sel]*1e-6;ds=np.atleast_1d(ds);ids=np.asarray(selected[pop]['source_bubble_ids'])[sel];ids=np.atleast_1d(ids);n=len(ds)
 sp={'name':name,'population':pop,'source_bubble_ids':ids.tolist(),'N':n,'diameters_m':ds.tolist(),'mpi_ranks':ranks,'engine':engine,'field':'fields/'+({'uniform':'UNIFORM_FIELD.h5','linear':'LINEAR_FIELD.h5','real':'FROZEN_FLOW_FIELD_V0.h5'}[kind]),'field_kind':kind,'dt_s':dt,'steps':steps if steps is not None else int(math.ceil(5*float(tau(ds[0]))/dt)),'rho_technical_kg_m3':1000.,'mu_pa_s':.001,'initial_velocity_m_s':[0.,0.,0.],'trajectory_stride':25 if n>1 else 1,'pair_interaction':'ZERO','role':'COUPLING_ENGINEERING_TEST_ONLY'}
 sp['field_sha256']=sha(P/sp['field']);sp['canonical_real_field_sha256']=sha(P/'fields/FROZEN_FLOW_FIELD_V0.h5')
 if kind=='real':xyz=spawn[None,:];boxlo=real.origin;boxhi=real.origin+(real.dims-1)*real.dx;sp['safety']=safety
 elif n>1:xyz=(np.array(np.unravel_index(np.arange(n),(10,10,10))).T+.5)*8e-6-40e-6;boxlo=np.array([-40e-6]*3);boxhi=-boxlo
 else:xyz=np.zeros((1,3));boxlo=np.array([-40e-6]*3);boxhi=-boxlo
 sp['initial_positions_m']=xyz.tolist() if n==1 else {'grid_shape':[10,10,10],'spacing_m':8e-6,'box_min_m':[-40e-6]*3};sp['box_min_m']=boxlo.tolist();sp['box_max_m']=boxhi.tolist()
 p=P/'cases'/name;p.mkdir(parents=True);put(p/'CASE_CONTRACT.json',sp)
 with (p/'particles.data').open('w') as f:
  f.write(f'Frozen SonoVue continuous sampler diameters; SI coupling engineering test\n\n{n} atoms\n1 atom types\n\n')
  for a,k in enumerate('xyz'):f.write(f'{boxlo[a]:.17g} {boxhi[a]:.17g} {k}lo {k}hi\n')
  f.write('\nAtoms # sphere\n\n');np.savetxt(f,np.column_stack((np.arange(1,n+1),np.ones(n),ds,np.full(n,1000.),xyz)),fmt=['%d','%d']+['%.17g']*5)
  f.write('\nVelocities\n\n');np.savetxt(f,np.column_stack((np.arange(1,n+1),np.zeros((n,6)))),fmt=['%d']+['%.17g']*6)
 relative=os.path.relpath(P/sp['field'],p);safe=sp.get('safety');args=[wall,*spawn,3*real.dx] if safe else [-1,0,0,0,0]
 cmds=['units si','atom_style sphere','boundary p p p','newton off',f'processors {ranks} 1 1','read_data particles.data','comm_modify vel yes','pair_style zero 6e-6','pair_coeff * *','neighbor 1e-6 bin','neigh_modify every 1 delay 0 check no','fix integrator all nve/sphere',f'fix coupling all frozen/flow/drag {relative} 0.001 {sp["trajectory_stride"]} trajectory '+' '.join(format(x,'.17g') for x in args),f'timestep {dt:.17g}','thermo '+str(max(1,sp['steps']//10)),'thermo_style custom step atoms time ke','thermo_modify lost error norm no format float %.17g','run 0','write_dump all custom initial.dump id diameter mass xu yu zu vx vy vz fx fy fz modify sort id format float %.17g',f'run {sp["steps"]}','write_dump all custom final.dump id diameter mass xu yu zu vx vy vz fx fy fz modify sort id format float %.17g','print "COUPLING_CASE_COMPLETED"']
 (p/'in.lammps').write_text('\n'.join(cmds)+'\n');allcases.append(sp)
add('case_a_uniform','A',0);add('case_b_linear','B',0,'linear');add('case_c_real','C',0,'real',steps=csteps)
for i,q in enumerate(['d10','d50','d90']):add('case_d_'+q,'D',i)
add('case_e_mpi1','E',slice(None),ranks=1,steps=1000);add('case_e_mpi4','E',slice(None),ranks=4,steps=1000);add('kokkos_compatibility','E',slice(None),steps=1000,engine='gpu')
gates={'coordinate_m':1e-12,'uniform_sampler_m_s':1e-13,'linear_sampler_m_s':1e-12,'cpp_python_sampler_m_s':1e-12,'real_node_identity_m_s':1e-12,'force_relative':1e-12,'force_zero_absolute_N':1e-25,'analytic_velocity_normalized':1e-3,'analytic_position_normalized':1e-3,'linear_velocity_normalized':2e-3,'linear_position_normalized':2e-3,'mpi_position_m':1e-10,'mpi_velocity_m_s':1e-10,'diameter_roundtrip_um':1e-10,'final_slip_over_initial_max':.01,'nonfinite':0,'lost_atoms':0,'invalid_query':0,'run_wall_seconds_max':300}
put(P/'contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json',{'status':'FROZEN_BEFORE_FORMAL_CASE_A','architecture':'FROZEN_FLOW_ONE_WAY','COUPLING_DIRECTION':'FLUID_TO_PARTICLE_ONLY','Palabos_runtime':False,'interpolation':'TRILINEAR, all 8 corners must be fluid and have velocity','invalid_query_policy':'HARD_FAIL; no clamping, nearest fallback or extrapolation; interpolation upper edge excluded if it cannot provide 8 distinct in-bounds corners','node_identity_test_selection':'1000 random fluid nodes with all adjacent 27 nodes valid, avoiding floating-point floor ambiguity near an integer boundary; no coordinate snapping','field_format':'FROZEN_FLOW_FIELD_V0','force_model':'TECHNICAL_STOKES_DRAG_V0','force_equation':'3*pi*mu*d*(u_f-v_p)','mu_pa_s':.001,'rho_technical_kg_m3':1000.,'particle_dt_s':dt,'cases':allcases,'gates':gates,'field_sha256':sha(P/'fields/FROZEN_FLOW_FIELD_V0.h5'),'source_particle_selection_contract_sha256':sha(P/'contracts/PARTICLE_SELECTION_CONTRACT.json'),'coupling_in_memory':True,'adapter_type':'external custom Fix registered via public Modify::fix_map in a thin embedding executable, linked to existing unmodified CPU/GPU static libraries; no LAMMPS core edit or baseline rebuild','force_evaluation_phase':'POST_FORCE uses current LAMMPS half-step velocity; unmodified nve/sphere updates final velocity. Trajectory records both actual applied force with its evaluation velocity and diagnostic Stokes force at completed velocity.','position_error_scale':'norm(initial fluid velocity)*tau_p; never absolute world coordinate norm','velocity_error_scale':'norm(initial fluid velocity)','real_field_scientific_status':'ENGINEERING_TRANSIENT_FIELD_ONLY','TWO_WAY_COUPLING':'OFF','WALL_PHYSICS':'OFF','RBC':'OFF','COUPLING_GPU_ACCELERATION':'PENDING','scientific_transport_model':'PENDING','automatic_solver_retry':False,'runtime_environment':{'OMP_NUM_THREADS':'1','HWLOC_COMPONENTS':'-gl','CUDA_VISIBLE_DEVICES':'0'},'gpu_flags':['-k','on','g','1','-sf','kk','-pk','kokkos','neigh','half','newton','off','gpu/aware','off']})
print(json.dumps({'CASE_PREPARATION':'PASS','cases':[{k:s[k] for k in ['name','N','steps','dt_s']} for s in allcases],'safe_interior':safety},indent=2))
