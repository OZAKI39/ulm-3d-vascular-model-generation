from pathlib import Path
import sys,json,hashlib,shutil
import numpy as np,h5py
R=Path(__file__).resolve().parents[1];ROOT=R/'cases';ROOT.mkdir(exist_ok=True)
# Every synthetic flow is specified before its trajectory is run.
def field(name,A):
 p=R/'fields'/name
 dims=np.array([31,31,31],np.int64);dx=2e-6;origin=np.full(3,-30e-6);ids=np.arange(np.prod(dims),dtype=np.uint64);ij=np.column_stack([ids%dims[0],ids//dims[0]%dims[1],ids//(dims[0]*dims[1])]);x=origin+dx*ij;U=x@np.array(A).T
 with h5py.File(p,'w') as f:
  f.attrs['format_version']='FROZEN_FLOW_FIELD_V0';f.attrs['coordinate_units']='m';f.attrs['velocity_units']='m/s';f.attrs['metadata_json']=json.dumps({'role':'synthetic validation only','analytic_gradient_s_inverse':np.array(A).tolist()})
  for n,d in [('dims',dims),('origin_m',origin),('dx_m',dx),('linear_index',ids),('fluid_linear_index',ids),('Velocity_m_s',U)]:f.create_dataset(n,data=d)
 return name
zero=field('ZERO_FIELD.h5',np.zeros((3,3)));normal=field('NORMAL_FIELD.h5',np.diag([-100,50,50]));shear=field('SHEAR_FIELD.h5',[[0,0,0],[1000,0,0],[0,0,0]]);strong=field('CONVERGING_FIELD.h5',np.diag([-100000,50000,50000]));three=field('THREE_FIELD.h5',[[-200,1000,0],[400,100,0],[0,0,100]])
old=json.loads((R/'provenance/CASE5_GEOMETRY_CONTRACT.json').read_text())
# Read exact original case contract from original evidence, locally or on Vast.
if (R/'provenance/old_case5/CASE_CONTRACT.json').exists():oldcase=R/'provenance/old_case5'
else:oldcase=Path('/workspace/microbubble_lammps/results/passive_transport_v0_20260915_233516/cases/case5_c025')
old=json.loads((oldcase/'CASE_CONTRACT.json').read_text())
basea=np.array([.65e-6,1.05e-6]);rr=basea.sum()+.05*basea.prod()/basea.sum();basepos=np.array([[-rr/2,0,0],[rr/2,0,0]])
index=[]
def case(name,a,x,field,dt,T,mpi=1,mode=0,drive=0,cgap=.4,ids=None,wall=False,gpu=False,source='synthetic'):
 d=ROOT/name;d.mkdir(exist_ok=True);a=np.array(a);x=np.array(x);ids=np.arange(1,len(a)+1) if ids is None else np.array(ids)
 lo=old['box_min_m'] if wall else [-32e-6]*3;hi=old['box_max_m'] if wall else [32e-6]*3
 data='Rigid sphere storage; placeholder density 1 kg/m3 never enters overdamped solve\n\n'+str(len(a))+' atoms\n1 atom types\n\n'
 for axis,k in zip('xyz',range(3)):data+=f'{lo[k]:.17g} {hi[k]:.17g} {axis}lo {axis}hi\n'
 data+='\nAtoms # sphere\n\n'
 for tag,ai,xx in zip(ids,a,x):data+=f'{tag} 1 {2*ai:.17g} 1 '+ ' '.join(format(v,'.17g') for v in xx)+'\n'
 data+='\nVelocities\n\n'
 for tag in ids:data+=f'{tag} 0 0 0 0 0 0\n'
 (d/'particles.data').write_text(data)
 cfg={'name':name,'N':len(a),'radii_m':a.tolist(),'initial_positions_m':x.tolist(),'ids':ids.tolist(),'field':'fields/'+field,'field_sha256':hashlib.sha256((R/'fields'/field).read_bytes()).hexdigest(),'dt_max':dt,'C_adv':.25,'C_gap':cgap,'max_time':T,'mpi_ranks':mpi,'drive_mode':mode,'drive_scale':drive,'drive_interpretation':'0=frozen field;1=constant external z torque/isolated rotational drag;2=opposite external z torques;3=frozen deterministic external forces and torques per ID','wall':wall,'gpu':gpu,'source':source,'steps_cap':200000,'stride':1,'formal_population':'NO; validation only'}
 (d/'CASE_CONTRACT.json').write_text(json.dumps(cfg,indent=2)+'\n');index.append(cfg)
 values={'field':'../../fields/'+field,'data':'particles.data','dt_max':dt,'c_gap':cgap,'max_time':T,'margin':old['margin_m'] if wall else 0,'stride':1,'drive_mode':mode,'drive_scale':drive,'steps_cap':200000,'wall':'../../provenance/closed_geometry_m.stl' if wall else 'NONE'}
 (d/'case.cfg').write_text(''.join(f'{k} {v}\n' for k,v in values.items()))
for cg in [.4,.2,.1]:
 suffix='cg'+str(int(cg*10)).zfill(2)
 case('A_normal_'+suffix,basea,basepos,normal,5e-5,.02,cgap=cg)
 case('F_contact_'+suffix,basea,basepos,strong,5e-7,2e-4,cgap=cg)
case('B_tangential',basea,basepos,shear,1e-5,.002)
case('C1_equal_transverse',[1e-6,1e-6],[[-1.0125e-6,0,0],[1.0125e-6,0,0]],zero,1e-5,.001,mode=1,drive=1000)
case('C2_opposite_transverse',[1e-6,1e-6],[[-1.0125e-6,0,0],[1.0125e-6,0,0]],zero,1e-5,.001,mode=2,drive=1000)
case('D_pass_by',basea,basepos,shear,1e-5,.005)
case('E_arbitrary_pair',basea,basepos,zero,1e-5,.001,mode=3,drive=.0002)
a=np.array([.65e-6,1.05e-6,.85e-6]);s01=a[0]+a[1]+.03*a[0]*a[1]/(a[0]+a[1]);s02=a[0]+a[2]+.03*a[0]*a[2]/(a[0]+a[2]);s12=a[1]+a[2]+.03*a[1]*a[2]/(a[1]+a[2]);theta=np.arccos((s01*s01+s02*s02-s12*s12)/(2*s01*s02));x=np.array([[0,0,0],[s01,0,0],[s02*np.cos(theta),s02*np.sin(theta),0]]);x-=x.mean(axis=0)
case('G_three_particle',a,x,three,1e-5,.002)
case('H_permuted',a,x,three,1e-5,.002,ids=[3,1,2])
case('I_mpi1',basea,basepos,shear,1e-5,.002,mpi=1)
case('I_mpi4',basea,basepos,shear,1e-5,.002,mpi=4)
# Isolated original real-field case 2 reproduces old passive trajectory for all 3010 steps.
case('J_passive',np.array([1.9429000435380954e-6/2]),[[.0001433584056461481,.00012636571998308154,.0001245268122626217]],'FROZEN_FLOW_FIELD_V0.h5',5.791749420578671e-5,.17433165755941799,wall=True,source='original passive case2')
case('K_real_8',np.array(old['diameters_um'])*.5e-6,old['initial_positions_m'],'FROZEN_FLOW_FIELD_V0.h5',old['dt_s'],.02,wall=True,source='exact original passive case5; unchanged population and field')
case('kokkos_smoke',basea,basepos,shear,1e-5,.0005,gpu=True)
(R/'contracts/CASE_INDEX.json').write_text(json.dumps(index,indent=2)+'\n')
(R/'contracts/CASE_C3_TWIST_STATUS.json').write_text(json.dumps({'status':'BLOCKED_TWIST_NOT_IMPLEMENTED','reason':'Non-singular twist collective/self contributions and matching to retained drag not independently validated','no_twist_case_run':True},indent=2)+'\n')
print('PREPARED',len(index),'cases; C3 explicitly blocked')
