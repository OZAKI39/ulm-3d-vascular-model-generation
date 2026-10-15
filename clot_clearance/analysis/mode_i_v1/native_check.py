from common import *
import os,sys,time,resource
sys.path.insert(0,str(R))
from analysis.dual_track_v1.native_a import Native
from analysis.dual_track_v1.core_a import Energy,contact
start=time.perf_counter();assert sha(B/'liblammps.so')==EXPECTED and not (O/'NATIVE_STATIC_CHECK.json').exists()
rt=np.load(O/'reference_topology.npz');X=rt['X'];b=rt['bonds'];n=len(X);alive=np.ascontiguousarray(rt['alive'][0],dtype=np.uint8)
g={'X':X,'bonds':b,'bond_types':np.ones(len(b),int),'anchors':np.zeros(n,bool)};e=Energy(g)
native=Native(O/'native_static_identity',g,alive);L=native.L
paths=sorted({l.split()[-1] for l in Path('/proc/self/maps').read_text().splitlines() if 'liblammps.so' in l});assert len(paths)==1 and sha(paths[0])==EXPECTED and Path(paths[0]).resolve()==(B/'liblammps.so').resolve()
m0=L.numpy.extract_atom('d_bpm_peri_reference_m')[:n].copy();assert np.max(abs(m0-e.m)/e.m)<1e-12
states=[('initial_A0',X),('unaccepted_elastic_endpoint',np.load(O/'elastic_precheck.npz')['x'])];records=[]
for name,x in states:
 f,info,th=native.evaluate(x);u,fi=e.eval(x,alive,True);uc,fc,cp=contact(x,g,alive);mi=L.numpy.extract_atom('d_bpm_peri_reference_m')[:n].copy();thi=e.theta(x,alive);diff=np.linalg.norm(f-fi-fc,axis=1)
 rec={'name':name,'U_LPS_J':u,'contact_U_J_independent':uc,'native':info,'max_internal_force_N':float(np.linalg.norm(f,axis=1).max()),'m_relative_max':float(np.max(abs(mi-e.m)/e.m)),'native_vs_independent_force_max_N':float(diff.max()),'native_vs_independent_force_relative':float(np.linalg.norm(f-fi-fc)/max(np.linalg.norm(fi+fc),1e-30)),'theta_max_abs_error':float(abs(th-thi).max()),'atom_IDs_exact':bool(np.array_equal(L.numpy.extract_atom('id')[:n],np.arange(1,n+1))),'velocity_max_m_s':float(abs(L.numpy.extract_atom('v')[:n]).max())}
 if name=='unaccepted_elastic_endpoint':
  ext=np.load(O/'elastic_precheck.npz')['external_force_N'];rec['max_total_residual_force_N']=float(np.linalg.norm(f+ext,axis=1).max());rec['residual_L2_relative']=float(np.linalg.norm(f+ext)/np.linalg.norm(ext))
 np.savez_compressed(O/f'native_{name}.npz',X=X,x=x,native_force_N=f,independent_LPS_force_N=fi,contact_force_N=fc,theta_native=th,theta_independent=thi,reference_m_native=mi,alive=alive,atom_ids=np.arange(1,n+1),velocity=L.numpy.extract_atom('v')[:n].copy())
 records.append(rec)
 assert rec['m_relative_max']<1e-12 and diff.max()<1e-16 and rec['theta_max_abs_error']<1e-12
native.close()
dump(O/'NATIVE_STATIC_CHECK.json',{'status':'PASS_STATIC_FORCE_IDENTITY_ONLY','actual_loaded_library':{p:sha(p) for p in paths},'one_native_instance':True,'intact_reference_run0':1,'static_evaluations_without_integrator':2,'formal_fracture_runs':0,'breaks_enabled':False,'reason':'Read-only coordinate evaluations to verify frozen m, initial reference equilibrium, and the unaccepted independent endpoint; no physical loading trajectory generated.','records':records,'wall_s':time.perf_counter()-start,'CPU_user_s':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'CPU_system_s':resource.getrusage(resource.RUSAGE_SELF).ru_stime,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
print(json.dumps(records,indent=2),flush=True)
