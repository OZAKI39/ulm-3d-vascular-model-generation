from pathlib import Path
import json,hashlib
import numpy as np,h5py
R=Path(__file__).resolve().parents[1];C=json.loads((R/'provenance/reference_audit/WALL_REFERENCE_CONVENTION.json').read_text());a=C['radii_m']['d50'];cases=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def field(name,u,A=None):
 dims=np.array([61,41,41],np.int64);dx=2e-6;origin=np.array([-60e-6,-40e-6,-20e-6]);ids=np.arange(np.prod(dims),dtype=np.uint64);ij=np.column_stack([ids%dims[0],ids//dims[0]%dims[1],ids//(dims[0]*dims[1])]);x=origin+ij*dx;velocity=np.tile(np.asarray(u,dtype=float),(len(x),1))
 if A is not None:velocity+=x@np.array(A).T
 p=R/'fields'/name
 with h5py.File(p,'w') as f:
  f.attrs['format_version']='FROZEN_FLOW_FIELD_V0';f.attrs['coordinate_units']='m';f.attrs['velocity_units']='m/s';f.attrs['metadata_json']=json.dumps({'role':'TECHNICAL_SYNTHETIC_TEST_ONLY','constant_velocity':u,'gradient_s_inverse':A})
  for name,data in [('dims',dims),('origin_m',origin),('dx_m',dx),('linear_index',ids),('fluid_linear_index',ids),('Velocity_m_s',velocity)]:f.create_dataset(name,data=data,compression='gzip' if np.ndim(data) else None)
field('WALL_ZERO_FIELD.h5',[0,0,0]);field('WALL_MPI_FIELD.h5',[.002,0,0]);field('WALL_COMBINED_FIELD.h5',[0,0,0],[[-1000,0,0],[0,500,0],[0,0,500]])
def make(name,radii,positions,fieldname,T,dt=5e-4,cwall=.4,hydro=True,hard=True,control=False,drive=(0,0,0),torque=(0,0,0),mpi=1,pattern=0,pair=False,wall='FLAT',gpu=False,box=None,ids=None,source='TECHNICAL_SYNTHETIC_TEST_ONLY'):
 d=R/'cases'/name;d.mkdir();radii=np.asarray(radii);positions=np.asarray(positions);ids=list(range(1,len(radii)+1)) if ids is None else ids
 lo,hi=box if box else ([-50e-6,-30e-6,-10e-6],[50e-6,30e-6,50e-6])
 text='Overdamped sphere storage; density placeholder never enters equations\n\n'+str(len(radii))+' atoms\n1 atom types\n\n'
 for j,axis in enumerate('xyz'):text+=f'{lo[j]:.17g} {hi[j]:.17g} {axis}lo {axis}hi\n'
 text+='\nAtoms # sphere\n\n'
 for tag,rr,x in zip(ids,radii,positions):text+=f'{tag} 1 {2*rr:.17g} 1 '+' '.join(f'{v:.17g}' for v in x)+'\n'
 text+='\nVelocities\n\n'+''.join(f'{tag} 0 0 0 0 0 0\n' for tag in ids);(d/'particles.data').write_text(text)
 cfg={'field':'../../fields/'+fieldname,'data':'particles.data','wall':wall,'table':'../../tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5','dt_max':dt,'c_gap':.4,'c_wall':cwall,'max_time':T,'stride':1,'steps_cap':200000,'wall_resistance':int(hydro),'hard_wall':int(hard),'pair_enabled':int(pair),'control_overlap':int(control),'drive_pattern':pattern,'drive_x':drive[0],'drive_y':drive[1],'drive_z':drive[2],'torque_x':torque[0],'torque_y':torque[1],'torque_z':torque[2]}
 (d/'case.cfg').write_text(''.join(f'{k} {v}\n' for k,v in cfg.items()))
 contract={'name':name,'status':'FROZEN_BEFORE_RUNTIME','N':len(radii),'radii_m':radii.tolist(),'initial_positions_m':positions.tolist(),'ids':ids,'box_min_m':lo,'box_max_m':hi,'config':cfg,'mpi_ranks':mpi,'kokkos':gpu,'source':source,'field_sha256':sha(R/'fields'/fieldname),'table_sha256':sha(R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5'),'particle_input_sha256':sha(d/'particles.data'),'C_adv':.25,'drive_units':'du m/s = external force / isolated bulk drag; dw rad/s = external torque / isolated bulk rotational drag','wall_background_closure':'excess*(q-q_inf); no ambient strain wall forcing'}
 (d/'CASE_CONTRACT.json').write_text(json.dumps(contract,indent=2)+'\n');cases.append(contract)
for cw in [.4,.2,.1]:make('F_hard_cw'+str(int(cw*10)),[a],[[0,0,1.2*a]],'WALL_ZERO_FIELD.h5',.003,cwall=cw,drive=(.001,0,-.005),torque=(0,20,0))
make('F_no_wall',[a],[[0,0,1.2*a]],'WALL_ZERO_FIELD.h5',1.3*.2*a/.005,dt=1e-6,hydro=False,hard=False,control=True,drive=(.001,0,-.005),torque=(0,20,0))
make('F_resistance_only',[a],[[0,0,1.2*a]],'WALL_ZERO_FIELD.h5',.003,hard=False,control=True,drive=(.001,0,-.005),torque=(0,20,0))
make('F_handoff',[a],[[0,0,20.9*a]],'WALL_ZERO_FIELD.h5',.001,dt=1e-5,drive=(.001,0,.002))
r=np.full(128,C['radii_m']['d10']);x=np.array([[(i%16-7.5)*3e-6,(i//16-3.5)*3e-6,r[i]*(1.1 if i%2==0 else 23.)] for i in range(128)])
for mpi in [1,4]:make('MPI'+str(mpi)+'_128',r,x,'WALL_MPI_FIELD.h5',.002,dt=1e-4,mpi=mpi,pattern=1,drive=(.0001,0,-.0002),torque=(0,10,0),pair=True)
make('KOKKOS_SMOKE',r[:8],x[:8],'WALL_MPI_FIELD.h5',.0002,dt=1e-4,gpu=True,pattern=1,drive=(.0001,0,-.0002),torque=(0,10,0),pair=True)
sep=2*a+.02*a/2
make('COMBINED_PAIR_WALL',[a,a],[[-sep/2,0,1.1*a],[sep/2,0,1.1*a]],'WALL_COMBINED_FIELD.h5',.0002,dt=1e-5,drive=(0,0,-.005),pair=True)
old=json.loads((R/'provenance/PREVIOUS_REAL_8_CASE_CONTRACT.json').read_text());old_data=json.loads((Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['pair_baseline'])/'provenance/old_case5/CASE_CONTRACT.json').read_text());box=(old_data['box_min_m'],old_data['box_max_m'])
make('J_real_8',old['radii_m'],old['initial_positions_m'],'FROZEN_FLOW_FIELD_V0.h5',.02,dt=old['dt_max'],pair=True,wall='../../provenance/closed_geometry_m.stl',box=box,ids=old['ids'],source='EXACT_FROZEN_PREVIOUS_REAL_8_POPULATION; geometry gate may block before timestep')
first=json.loads((R/'contracts/CASE_I_INITIAL_SELECTION.json').read_text());extra=json.loads((R/'contracts/CASE_I_SUPPLEMENTAL_SELECTION.json').read_text());chosen=first['chosen'] or extra['chosen']
if chosen:
 make('I_real_single',[chosen['radius_m']],[chosen['position_m']],'FROZEN_FLOW_FIELD_V0.h5',.02,dt=old['dt_max'],wall='../../provenance/closed_geometry_m.stl',box=box,source='D50_INITIAL_GAP_0P05_TO_0P2_STATIC_SELECTED')
else:
 d=R/'cases/I_real_single';d.mkdir();state={'status':'BLOCKED_NO_VALID_NEAR_WALL_START','reason':'No frozen d50 local-plane-valid and flow-stencil-valid candidate in the recorded static searches','accepted_steps':0,'time_s':0.,'target_time_s':.02,'solver_launched':False,'radius_changed':False,'gate_relaxed':False};(d/'RUN_STATE.json').write_text(json.dumps(state,indent=2)+'\n')
(R/'contracts/CASE_INDEX.json').write_text(json.dumps(cases,indent=2)+'\n')
contract={'contract_name':'MICROBUBBLE_WALL_HYDRODYNAMICS_V0','status':'FROZEN_BEFORE_RUNTIME','model':'RMBW_LOOKUP_LOCAL_PLANE_RESISTANCE_V0',
 'table_sha256':sha(R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5'),'table_contract_sha256':sha(R/'contracts/RMBW_WALL_TABLE_CONTRACT.json'),
 'mu_Pa_s':.001,'rho_kg_m3':1000,'particle_inertia':False,'RBC':False,'buoyancy':False,'lift':False,'adhesion':False,'bubble_bubble_twist':'PENDING',
 'flow_sha256':sha(R/'fields/FROZEN_FLOW_FIELD_V0.h5'),'geometry_sha256':sha(R/'provenance/closed_geometry_m.stl'),'real_flow_status':'ENGINEERING_TRANSIENT_FIELD_ONLY',
 'matrix_convention_sha256':sha(R/'provenance/reference_audit/WALL_REFERENCE_CONVENTION.json'),'bulk_count':1,'wall_assembly':'R_global=R_bulk+R_pair_excess+R_wall_excess; b_wall=R_wall_excess*q_inf',
 'wall_background_closure':'RESISTANCE_ONLY; ambient strain-induced wall force NOT_INCLUDED_IN_V0','interpolation':'piecewise linear log epsilon; full6x6; memory resident',
 'range':[.001,20.],'near_field_high_confidence':[.001,.2],'qualified_reference':[.001,5.],'limited_reference':[5.,20.],'minimum_gap_regularization':'.001*a; numerical only','above_range':'zero excess; DEVELOPMENT BULK HANDOFF',
 'local_plane':json.loads((R/'contracts/GEOMETRY_AUDIT_PLAN.json').read_text()),'tangent_basis':'least parallel global axis, Gram-Schmidt, right-handed; simultaneous tangent sign chosen to maximize dot with previous t1',
 'wall_constraint':'base_gap+dt*base_normal dot V>=0; unified deterministic active set; wall by particleID/triangleID then pair minID/maxID',
 'swept_acceptance':'continuous segment sphere/triangle gap verified using BVH 1-Lipschitz distance bound with adaptive interval subdivision; depth cap30 rejects, never accepts unproven interval',
 'C_adv':.25,'C_gap':.4,'C_wall_candidates':[.4,.2,.1],'C_wall_reference':.1,'selection_rule':'largest passing all Cwall gates',
 'gates':{'interpolation_action_relative':1e-3,'reciprocity_relative':1e-10,'scaled_min_eigenvalue':0.,'geometry_distance_m':1e-12,'rotation_identity_relative':1e-12,'accepted_overlap_m':-1e-12,'Cwall_position_m':1e-8,'Cwall_gap_m':1e-8,'Cwall_omega_relative_to_peak':.02,'Cwall_integrated_rotation_rad':1e-3,'MPI_position_m':1e-10,'MPI_velocity_m_s':1e-10,'MPI_omega_rad_s':1e-7,'linear_residual':1e-10,'RK2_reference_position_m':1e-9,'RK2_reference_velocity_m_s':1e-9,'RK2_reference_omega_rad_s':1e-6},
 'seeds':{'table':2026091601,'geometry':2026091602,'rotation':2026091603},'case_index_sha256':sha(R/'contracts/CASE_INDEX.json'),'case_J_duration_s':.02,'case_I_duration_s':.02,'initial_conditions':'explicit case contracts; CaseI may remain BLOCKED without a valid start',
 'reference_limitations':json.loads((R/'contracts/RMBW_WALL_TABLE_CONTRACT.json').read_text())['limitations'],'maximum_possible_success':'PASS_WITH_REFERENCE_LIMITATIONS','human_visual_review':'PENDING'}
(R/'contracts/MICROBUBBLE_WALL_HYDRODYNAMICS_V0_CONTRACT.json').write_text(json.dumps(contract,indent=2)+'\n');print('CASES_FROZEN',len(cases),'I_valid_start',bool(chosen))
