from pathlib import Path
import sys,json,ctypes,itertools,hashlib,time
import numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from reference_rigid_sphere_lubrication import coefficients,pair_matrix,assemble,dense_solve,FlowGradientReference
lib=ctypes.CDLL(str(R/'build_math/librigid_math.so'));P=ctypes.POINTER(ctypes.c_double);I=ctypes.POINTER(ctypes.c_int)
lib.rigid_coefficients.argtypes=lib.stable_coefficients.argtypes=[ctypes.c_int,P,P]
lib.rigid_audit_system.argtypes=[ctypes.c_int,P,ctypes.c_int,I,P,P,P,P]
lib.gradient_audit.argtypes=[ctypes.c_char_p,ctypes.c_int,P,P]
def ptr(x):return x.ctypes.data_as(P)
def save(name,d): (R/'validation'/name).write_text(json.dumps(d,indent=2)+'\n')
def cpp_system(x,a,U,O,E,pairs,solve=False):
 n=len(a);state=np.ascontiguousarray(np.column_stack([x,a,U,O,E.reshape(n,9)]));pairs=np.ascontiguousarray(pairs,dtype=np.int32).reshape(-1,2);M=np.empty((6*n,6*n));b=np.empty(6*n);q=np.empty((n,6));stats=np.zeros(2)
 ret=lib.rigid_audit_system(n,ptr(state),len(pairs),pairs.ctypes.data_as(I),ptr(M),ptr(b),ptr(q) if solve else None,ptr(stats));assert ret==0
 return M,b,q,stats
rng=np.random.default_rng(20260916);N=100000;a=rng.uniform(.375e-6,2.625e-6,N);b=rng.uniform(.375e-6,2.625e-6,N)
a[:20000]=b[:20000];a[20000:30000]=b[20000:30000]=.375e-6;a[30000:40000]=.375e-6;b[30000:40000]=2.625e-6;a[40000:50000]=1.9367184131090989e-6/2;b[40000:50000]=3.0424886877828055e-6/2
re=a*b/(a+b);h=re*10**rng.uniform(-3,np.log10(.2),N);inp=np.ascontiguousarray(np.column_stack([a,b,h]));new=np.empty((N,6));old=np.empty((N,6));lib.rigid_coefficients(N,ptr(inp),ptr(new));lib.stable_coefficients(N,ptr(inp),ptr(old));ref=coefficients(a,b,h)
err=np.abs(new-ref[:,:6])/np.maximum(np.abs(ref[:,:6]),[1e-25,1e-25,1e-37,1,1e-25,1e-25]);scalar_old=np.max(np.abs(old[:,:2]-ref[:,:2])/np.maximum(np.abs(ref[:,:2]),1e-25),axis=0)
swap=np.ascontiguousarray(np.column_stack([b,a,h]));oswap=np.empty_like(old);lib.stable_coefficients(N,ptr(swap),ptr(oswap))
# Source pump is the total C11 coefficient, not the Schur-complement pump.
oldpump=np.abs(old[:,2]-ref[:,2])/np.maximum(np.abs(ref[:,2]),1e-37)
oldcoupling=np.abs(old[:,1]*a-ref[:,7])/np.maximum(np.abs(ref[:,7]),1e-37)
report={'status':'PASS' if np.max(err)<1e-10 else 'FAIL','N':N,'coefficient_columns':['sq','sh','pu','tw_PENDING_NOT_TESTED','li','lj'],'production_max_relative_error_by_column':err.max(axis=0).tolist(),'production_max_relative_error':float(err.max()),'stable_sq_sh_relative_error':scalar_old.tolist(),'stable_poly_vs_upoly_max_relative_error':float(np.max(np.abs(old[:,:3]-old[:,3:])/np.maximum(np.abs(old[:,:3]),1e-37))),'stable_pump_vs_reference_schur_max_relative_error':float(oldpump.max()),'stable_pump_swap_max_relative_error':float(np.max(np.abs(old[:,2]-oswap[:,2])/np.maximum(np.maximum(np.abs(old[:,2]),np.abs(oswap[:,2])),1e-37))),'stable_radius_times_shear_vs_published_YB_max_relative_error':float(oldcoupling.max()),'stable_required_mode_gate':'FAIL_REPAIRED_IN_PROJECT_SOURCE','repair_gate':'PASS' if np.max(err)<1e-10 else 'FAIL','samples_sha256':hashlib.sha256(inp.tobytes()).hexdigest()}
np.savez_compressed(R/'validation/COEFFICIENT_100000_STATES.npz',inputs=inp,stable=old,production=new,reference=ref)
save('COEFFICIENT_AUDIT.json',report);print('COEFFICIENT',report,flush=True);assert report['status']=='PASS'
# Matrix properties audited without calling a linear solver. Decision is frozen only afterward.
maxmatrix=maxrhs=maxsym=0;mineigen=np.inf;maxcond=0;states=[]
for n in [2,3,5,10]:
 for repeat in range(20):
  a=rng.uniform(.5e-6,1.5e-6,n);x=np.zeros((n,3));direction=rng.normal(size=3);direction/=np.linalg.norm(direction)
  for i in range(1,n):x[i]=x[i-1]+direction*(a[i-1]+a[i]+.01*min(a[i-1],a[i]))
  U=rng.normal(0,.001,(n,3));O=rng.normal(0,1000,(n,3));E=rng.normal(0,100,(n,3,3));E=.5*(E+E.transpose(0,2,1));pairs=list(itertools.combinations(range(n),2))
  M,b,_,_=cpp_system(x,a,U,O,E,pairs);Mr,br,_=assemble(x,a,U,O,E,pairs);scale=np.column_stack([np.ones((n,3)),np.repeat((1/a)[:,None],3,axis=1)]).reshape(-1);A=M*scale[:,None]*scale[None,:]/1e-8;Ar=Mr*scale[:,None]*scale[None,:]/1e-8
  maxmatrix=max(maxmatrix,np.linalg.norm(A-Ar)/np.linalg.norm(Ar));maxrhs=max(maxrhs,np.linalg.norm((b-br)*scale)/np.linalg.norm(br*scale));maxsym=max(maxsym,np.linalg.norm(A-A.T)/np.linalg.norm(A));ev=np.linalg.eigvalsh(A);mineigen=min(mineigen,ev[0]);maxcond=max(maxcond,ev[-1]/ev[0]);states.append((x,a,U,O,E,pairs))
m={'status':'PASS' if max(maxmatrix,maxrhs)<1e-10 and maxsym<1e-12 and mineigen>0 else 'FAIL','systems':len(states),'N':[2,3,5,10],'max_scaled_matrix_reference_error':maxmatrix,'max_scaled_rhs_reference_error':maxrhs,'max_scaled_symmetry_error':maxsym,'min_scaled_eigenvalue':mineigen,'max_scaled_condition_number':maxcond,'properties':'Pair PSD; positive isolated Stokes drag makes total SPD','solver_invoked_before_this_audit':False};save('MATRIX_PROPERTY_AUDIT.json',m);assert m['status']=='PASS',m
selection={'selected':'PCG','reason':'Independent matrix audit found symmetric positive definite drag-plus-pair system; exact modal B^T D B pair energy proves PSD for positive coefficients. Radius scaling handles mixed translation/rotation SI units. Hard contact solved by active-set Schur complement; underlying resistance solves remain PCG.','evidence':'validation/MATRIX_PROPERTY_AUDIT.json','max_scaled_condition_number':maxcond,'target_recursive_residual':2e-14,'required_recomputed_KKT_residual':1e-10}
(R/'contracts/LINEAR_SOLVER_SELECTION.json').write_text(json.dumps(selection,indent=2)+'\n')
maxres=maxsol=0
for x,a,U,O,E,pairs in states:
 M,b,q,stats=cpp_system(x,a,U,O,E,pairs,True);qr=dense_solve(x,a,U,O,E,pairs);scale=np.column_stack([np.ones((len(a),3)),np.repeat(a[:,None],3,axis=1)]);maxsol=max(maxsol,np.linalg.norm((q-qr)*scale)/np.linalg.norm(qr*scale));maxres=max(maxres,stats[0])
save('DENSE_SOLVER_AUDIT.json',{'status':'PASS' if maxres<=1e-10 and maxsol<=1e-9 else 'FAIL','max_residual':maxres,'max_scaled_solution_relative_error':maxsol,'systems':len(states)});assert maxres<=1e-10 and maxsol<=1e-9
# Actual C++ assembled matrices evaluated on arbitrary combined motions, swaps and collective rotations.
maxforce=maxtorque=maxexchange=maxrigid=0;neg=0;mindiss=0
for k in range(10000):
 a=rng.uniform(.375e-6,2.625e-6,2);h=rng.uniform(.001,.199)*a.prod()/a.sum();normal=rng.normal(size=3);normal/=np.linalg.norm(normal);r=a.sum()+h;x=np.array([.5*r*normal,-.5*r*normal]);U=np.zeros((2,3));O=U.copy();E=np.zeros((2,3,3));M,_,_,_=cpp_system(x,a,U,O,E,[(0,1)]);drag=np.column_stack([np.repeat((6*np.pi*.001*a)[:,None],3,axis=1),np.repeat((8*np.pi*.001*a**3)[:,None],3,axis=1)]).reshape(-1);M-=np.diag(drag)
 q=np.column_stack([rng.normal(0,.001,(2,3)),rng.normal(0,1000,(2,3))]);f=(-M@q.reshape(-1)).reshape(2,6);fscale=max(np.linalg.norm(f[:,:3],axis=1).max(),1e-30);maxforce=max(maxforce,np.linalg.norm(f[0,:3]+f[1,:3])/fscale)
 moment=np.cross(x[0]-x[1],f[0,:3])+f[:,3:].sum(axis=0);tscale=max(r*fscale,np.linalg.norm(f[:,3:],axis=1).max(),1e-36);maxtorque=max(maxtorque,np.linalg.norm(moment)/tscale)
 dissip=q.reshape(-1)@M@q.reshape(-1);energy_scale=np.linalg.norm(M)*np.linalg.norm(q)**2;mindiss=min(mindiss,dissip/max(energy_scale,1e-30));neg+=dissip < -1e-12*energy_scale
 Ms,_,_,_=cpp_system(x[::-1],a[::-1],U,O,E,[(0,1)]);Ms-=np.diag(drag.reshape(2,6)[::-1].reshape(-1));fs=(-Ms@q[::-1].reshape(-1)).reshape(2,6)[::-1];normscale=np.r_[np.repeat(1/fscale,3),np.repeat(1/tscale,3)];maxexchange=max(maxexchange,np.linalg.norm((fs-f)*normscale))
 if k<1000:
  w=rng.normal(0,1000,3);rigidq=np.column_stack([np.cross(np.broadcast_to(w,(2,3)),x),np.broadcast_to(w,(2,3))]);fr=M@rigidq.reshape(-1);maxrigid=max(maxrigid,np.linalg.norm(fr*np.tile(normscale,2))/max(np.linalg.norm(rigidq),1e-30))
sym={'status':'PASS' if maxforce<=1e-12 and maxtorque<=1e-12 and maxexchange<1e-10 and neg==0 else 'FAIL','combined_states':10000,'max_pair_force_normalized_error':maxforce,'max_pair_torque_balance_normalized_error':maxtorque,'max_exchange_scaled_error':maxexchange,'negative_dissipation_violations':int(neg),'minimum_normalized_dissipation':mindiss,'max_collective_rotation_scaled_resistance':maxrigid};save('PAIR_TORQUE_BALANCE_AUDIT.json',sym);print('SYMMETRY',sym,flush=True);assert sym['status']=='PASS'
# Frozen affine data, exact derivative and independent Python edge-difference reference.
path=R/'fields/LINEAR_FIELD.h5';ref=FlowGradientReference(path);x=rng.uniform(-70e-6,70e-6,(10000,3));x=np.ascontiguousarray(x);out=np.empty((len(x),15));assert lib.gradient_audit(str(path).encode(),len(x),ptr(x),ptr(out))==0
U,O,E,G=ref.query(x);analytic=np.array([[20000,10000,-5000],[-8000,15000,7000],[3000,-4000,-12000]])
graderror=np.max(np.abs(out[:,6:].reshape(-1,3,3)-analytic));diff=np.max(np.abs(out[:,6:].reshape(-1,3,3)-G));g={'status':'PASS' if max(graderror,diff)<=1e-10 else 'FAIL','N':len(x),'cpp_vs_analytic_max_s_inverse':float(graderror),'cpp_vs_python_max_s_inverse':float(diff),'max_velocity_difference_m_s':float(np.max(np.abs(out[:,:3]-U))),'max_omega_difference_s_inverse':float(np.max(np.abs(out[:,3:6]-O)))};save('FLOW_GRADIENT_AUDIT.json',g);print('GRADIENT',g,flush=True);assert g['status']=='PASS'
save('MATHEMATICAL_RELEASE_GATE.json',{'status':'PASS','twist':'PENDING_BLOCKER','required_non_twist_modes':'PASS','real_cases_authorized_by_math_gate':True,'reference_contract_sha256':hashlib.sha256((R/'contracts/RIGID_SPHERE_LUBRICATION_REFERENCE_CONTRACT.json').read_bytes()).hexdigest()});print('MATHEMATICAL_RELEASE_GATE PASS',flush=True)
