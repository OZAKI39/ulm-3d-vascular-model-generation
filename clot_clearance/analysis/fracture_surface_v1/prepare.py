"""Read-only source/topology gates and preregistration. No ideal plane integration here."""
from common import *
from fractions import Fraction as Q
from datetime import datetime, timezone
import itertools, time, resource, ctypes, shutil

start=time.perf_counter()
assert not (O/'DESIGN_FREEZE.json').exists(), 'Never overwrite a frozen run'
def rel(a,b):
    return float(np.max(np.abs(a-b))/max(float(np.max(np.abs(b))),1e-300))
path=W/'reference_topology.npz'; data=np.load(path)
X,b,r,h,m=[data[k] for k in ['X','bonds','r0_m','h','m']]; n=len(X); nb=len(b)
manifest=R/'verification/y_clot_cfd_and_wide_pvce_dual_track_v1/20261011T154854Z/SERVER_ARTIFACT_MANIFEST.json'
entry=next(x for x in json.loads(manifest.read_text())['files'] if x['path']==str(path.relative_to(R)))
assert sha(path)==entry['sha256']==sha(O/'source_inputs/reference_topology.npz')
dump(O/'HISTORICAL_INPUT_MANIFEST_EXTRACT.json',{'manifest':str(manifest),'manifest_sha256':sha(manifest),'entry':entry})
old=json.loads((W/'DESIGN_FREEZE.json').read_text())
assert [old[k] for k in ['K_Pa','G_Pa','s00','alpha','dx_m','horizon_m']]==[K,G,S00,0.,DX,DELTA]
ids=np.arange(1,n+1); ijk=np.stack(np.meshgrid(np.arange(100),np.arange(8),np.arange(20),indexing='ij'),axis=-1).reshape(-1,3)
axes=[(np.arange(c)+.5)*DX-d/2 for c,d in zip([100,8,20],[.005,PERIOD,.001])]
reference=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3)
assert np.array_equal(X,reference)
assert np.array_equal(b[:,0]*np.int64(n)+b[:,1],np.unique(b[:,0]*np.int64(n)+b[:,1])) and np.all(b[:,0]<b[:,1])
assert np.array_equal(data['bond_reference_id'],np.arange(nb)+1)
saved=np.load(W/'preflight/full/reference.npz'); assert np.array_equal(saved['x'],X)
assert not np.any(saved['native_force']) and not np.any(saved['native_theta'])
keys=[]; offsets=[]
for off in itertools.product(range(-3,4),repeat=3):
    if not 0<sum(v*v for v in off)<=(DELTA/DX)**2:continue
    offsets.append(off); t=ijk+off; valid=(t[:,0]>=0)&(t[:,0]<100)&(t[:,2]>=0)&(t[:,2]<20)
    j=(t[:,0]*8+t[:,1]%8)*20+t[:,2]; i=np.arange(n); valid &= i<j; keys.append(i[valid]*n+j[valid])
expected=np.unique(np.concatenate(keys)); actual=b[:,0].astype(np.int64)*n+b[:,1]
assert np.array_equal(expected,actual)
del keys,expected
image=-np.rint((X[b[:,1],1]-X[b[:,0],1])/PERIOD).astype(np.int8)
d=X[b[:,1]]-X[b[:,0]]; d[:,1]+=image*PERIOD; rcheck=np.linalg.norm(d,axis=1)
hcheck=np.where(abs(rcheck-DELTA)<=DX/2,(DELTA+DX/2-rcheck)/DX,1.)
mc=np.bincount(b.ravel(),weights=np.repeat(rcheck*V*hcheck,2),minlength=n)
assert np.array_equal(image,data['periodic_image_y']) and rel(rcheck,r)<1e-13 and rel(hcheck,h)<1e-13 and rel(mc,m)<1e-13
assert PERIOD>2*DELTA and r.max()<=DELTA*(1+1e-12)
perm=((ijk[:,0]*8+(ijk[:,1]+1)%8)*20+ijk[:,2]); mapped=np.sort(perm[b],axis=1); mk=mapped[:,0].astype(np.int64)*n+mapped[:,1]; ix=np.searchsorted(actual,mk)
assert np.array_equal(actual[ix],mk) and rel(m[perm],m)<1e-13
roi=(X[:,0]>=-.00075)&(X[:,0]<.00075); influence=roi&(np.abs(X[:,2])<=DELTA+DX)
degrees=np.bincount(b.ravel(),minlength=n)
assert np.all(degrees[influence]==len(offsets)) and roi.sum()==4800
topology={'status':'PASS','source':str(path),'source_sha256':sha(path),'historical_manifest_match':True,'nodes':n,'unique_bonds':nb,'cross_Y_seam_bonds':int(np.count_nonzero(image)),'original_atom_ID_definition':'id = ((ix*8+iy)*20+iz)+1; exact match original writer/generator and stored native reference x ordered by ID','atom_ID_range':[1,n],'unique_reference_bond_IDs':True,'independent_integer_stencil_matches_all_original_pairs':True,'pristine_not_A0_A1_A2':True,'ignored_alive_mask_removed_counts':(data['alive']==0).sum(1),'reference_coordinates_exact':True,'native_saved_reference_force_zero':True,'uniform_volume_m3':V,'volume_source':'Original native input set atom * d_vfrac and DX**3 in original writer; no regenerated model used','original_weighted_m_range_m4':[m.min(),m.max()],'m_recomputed_relative_error':rel(mc,m),'r0_recomputed_relative_error':rel(rcheck,r),'taper_recomputed_relative_error':rel(hcheck,h),'periodic_image_exact':True,'periodic_Y_translation_bondset_invariant':True,'periodic_Y_translation_m_relative_error':rel(m[perm],m),'max_r0_m':r.max(),'period_m':PERIOD,'period_greater_than_2delta':True,'interior_directed_neighbors_per_node':len(offsets),'ROI_influence_nodes_all_have_full_stencil':True,'physical_cell_bounds_m':[[-.0025,.0025],[-.0002,.0002],[-.0005,.0005]],'point_center_bounds_m':[X.min(0),X.max(0)]}
dump(O/'PRISTINE_REFERENCE_TOPOLOGY_AUDIT.json',topology)
# Exact polynomial integration gives the paper constants independently of its printed result.
poly={0:Q(1,2),1:Q(-2,3),4:Q(1,6)}; square={}
for i,a in poly.items():
    for j,c in poly.items():square[i+j]=square.get(i+j,Q(0))+a*c
i1=sum(v/Q(k+1) for k,v in poly.items()); i2=sum(v/Q(k+1) for k,v in square.items())
assert i1==Q(1,5) and i2==Q(26,405)
beta=DELTA/5; bp=26*DELTA/405; form1=((9*K-15*G)*bp+15*G*beta)*S00**2
form2=DELTA*S00**2*(26*K/45+55*G/27); kterm=9*K*bp*S00**2; gterm=15*G*(beta-bp)*S00**2
assert abs(form1-form2)<1e-15
dump(O/'THEORY_FORMULA_AUDIT.json',{'status':'PASS','paper':'Zhang & Qiao 2018, DOI 10.1007/s10704-018-0285-8','PDF_pages':[3,4,5,6],'equations':[2,3,4,9,10,11,12,17,18,19,20,21,22,26,27],'three_dimensional':True,'omega':'delta/r (paper); 1/r in locked source','a_prime_Pa':K-5*G/3,'b_times_m_Pa':15*G,'lambda_cap_q':'1/2 - 2*q/3 + q**4/6; 0<=q=z/delta<=1','exact_beta_over_delta':str(i1),'exact_beta_prime_over_delta':str(i2),'beta_m':beta,'beta_prime_m':bp,'formula1':'((9*K-15*G)*beta_prime+15*G*beta)*s00**2','formula2':'delta*s00**2*((26/45)*K+(55/27)*G)','G0_theory_J_m2':form1,'formula2_J_m2':form2,'difference_J_m2':form1-form2,'K_contribution_J_m2':kterm,'G_contribution_J_m2':gterm,'K_fraction':kterm/form1,'G_fraction':gterm/form1,'K_G_split':'Positive volumetric .5*K*theta**2 and deviatoric term; expanded a_prime is NOT the K-only contribution','parameters':{'K_Pa':K,'G_shear_Pa':G,'s00':S00,'delta_m':DELTA,'dx_m':DX,'alpha':0.},'units':'Pa * m * stretch^2 = J/m^2; node energy density is J/m^3','limitations':['ideal simultaneous critical cross-plane point-state extensions','not necessarily globally kinematically realizable','not actual fracture dissipation','not bond-based or 2D formula','no parameter calibration or grid convergence']})
# One historical stored elastic state, re-evaluated with the EXISTING independent .so.
# This is a force/energy algebra check; it is not a new LAMMPS call or equilibrium solve.
q=np.load(W/'preflight/full/periodic_nonaffine.npz'); x=np.ascontiguousarray(q['x']); dd=x[b[:,1]]-x[b[:,0]]; dd[:,1]-=PERIOD*np.rint(dd[:,1]/PERIOD)
extension=np.linalg.norm(dd,axis=1)-r; extension[abs(extension)<2.2204e-16]=0
theta=3*V/m*np.bincount(b.ravel(),weights=np.repeat(h*extension,2),minlength=n)
devsum=np.zeros(n)
for end in [0,1]:
    ii=b[:,end]; ed=extension-theta[ii]*r/3
    devsum+=np.bincount(ii,weights=h/r*ed**2*V,minlength=n)
point_energy=.5*K*theta**2+7.5*G/m*devsum; U=float(np.sum(point_energy)*V)
library=ctypes.CDLL(str(R/'analysis/dual_track_v1/energy_y.so')); fun=library.evaluate_periodic; fun.restype=ctypes.c_double; fun.argtypes=[ctypes.c_int]*2+[ctypes.c_void_p]*6+[ctypes.c_double]*3+[ctypes.c_void_p,ctypes.c_double]
alive=np.ones(nb,np.uint8); force=np.zeros_like(x); args=[np.ascontiguousarray(a) for a in [x,b,alive,r,h,m]]
Ucpp=fun(n,nb,*[a.ctypes.data for a in args],V,K,G,force.ctypes.data,PERIOD)
err=abs(U-Ucpp)/abs(Ucpp); ferr=float(np.linalg.norm(force-q['native_force'])/np.linalg.norm(q['native_force'])); terr=rel(theta,q['native_theta'])
assert err<1e-12 and ferr<1e-9 and terr<1e-9
source={'status':'PASS_FOR_INTACT_REFERENCE_IDEAL_POINT_STATE','source_identity_file':'SOURCE_AND_LIBRARY_IDENTITY.json','source_mapping':{'reference_r0':'bond_bpm_peri.cpp:158-203; original minimum_image at store_data','taper':'bond_bpm_peri.cpp:83-88; h=(delta+dx/2-r0)/dx inside last half-cell; else 1; centers already selected r0<=delta','m':'bond_bpm_peri.cpp:506-546; full virgin sum r0*Vj*h; per-ID reference_m/ready properties persist across sort/migrate/restart','theta':'bond_bpm_peri.cpp:550-602; 3/m sum e*Vj*h; current surviving list, reference r0 and m','LPS_force':'bond_bpm_peri.cpp:355-371; equal-V LPS force density, multiply by V for nodal N','independent_energy':'energy_y.cpp:16-20, .5*(K-5G/3)*sum(V*theta^2)+7.5*G*V^2*sum_unique(h/r0*e^2*(1/mi+1/mj))'},'influence_constant_scaling':'omega_paper=delta*omega_native, m_paper=delta*m_native; lambda,theta, b*omega and energy unchanged if both numerator and m scaled; native m unit m^4, paper m unit m^5','extra_corrections':'No surface/shape tensor or Shojaei weights; equal V, uniform K/G, pristine full family only','point_energy':'W_i=.5*K*theta_i^2+7.5*G/m_i*sum_all_j(omega*Vj*h*(e_ij-theta_i*r0/3)^2)','force_correspondence':'Derivative of independent full-state energy yields V^2*h*[(3K-5G)*(theta_i/mi+theta_j/mj)+15G*e/r0*(1/mi+1/mj)] along current direction','native_PE':'Known incomplete reported deviatoric PE is NOT used; this is force-compatible independent energy','fixed_m_scope':'Eq2-to-Eq9 reduction here uses intact sum omega*r0^2*Vh=m_full. Do not claim the same identity over a deleted family with frozen full m. No damaged-state energy result added in this task.','saved_state_recheck':{'state':str(W/'preflight/full/periodic_nonaffine.npz'),'sha256':sha(W/'preflight/full/periodic_nonaffine.npz'),'point_Eq2_U_J':U,'existing_CXX_U_J':Ucpp,'relative_energy_error':err,'historical_native_force_relative_error':ferr,'historical_native_theta_scaled_error':terr,'new_native_calls':0,'historical_native_library_scope':'Original accepted wide preflight; final double fix source algebra inspected now. No new final native solve.'},'broken_state_lagged_theta':'Not exercised; existing time-layer limitations remain','current_source_bug_found':False}
dump(O/'SOURCE_LPS_CONSISTENCY_AUDIT.json',source)
cgroup={}
for p in ['/sys/fs/cgroup/cpu/cpu.cfs_quota_us','/sys/fs/cgroup/cpu/cpu.cfs_period_us','/sys/fs/cgroup/memory/memory.limit_in_bytes','/sys/fs/cgroup/cpuset/cpuset.cpus']:
    if Path(p).exists():cgroup[p]=Path(p).read_text().strip()
dump(O/'CGROUP_AUDIT.json',{'kind':'v1','values':cgroup,'CPU_quota_cores':int(cgroup['/sys/fs/cgroup/cpu/cpu.cfs_quota_us'])/int(cgroup['/sys/fs/cgroup/cpu/cpu.cfs_period_us']),'task_threads':1,'GPU_used':False})
freeze={'frozen_utc':datetime.now(timezone.utc).isoformat(),'run_id':C['run'],'preconditions':['paper full PDF read in WSL','theory algebra PASS','source force/independent energy correspondence PASS for complete family','pristine hash and independent enumeration PASS'],'reference_sha256':sha(path),'native_library_sha256':EXPECTED,'material':{'K_Pa':K,'G_Pa':G,'s00':S00,'alpha':0.,'dx_m':DX,'delta_m':DELTA,'V_m3':V},'ROI':{'x_half_open_m':[-.00075,.00075],'x_columns_zero_based':[35,64],'column_count':30,'L_ROI_m':30*DX,'Y_period_m':PERIOD,'A_m2':30*DX*PERIOD,'nodes':4800,'selection':'node center X in half-open interval, all Y and Z. Full node families include neighbors outside ROI; plane is infinite. Only node energies are windowed; no clipping of bonds at ROI edge.','minimum_ROI_edge_to_X_wall_m':.00175,'nonzero_energy_nodes_reference_family_untruncated':'Both planes fit within z physical +/-0.5 mm; verified complete 122-neighbor stencil on potential influence band.'},'planes_m':[0.,DX],'classification':'lower side iff reference z < plane; equal assigned upper; exact original float64 coordinates; actual min distance >24 um; unique undirected crossing iff endpoints differ in side','two_sided':'Sum all endpoint point energies on both sides once, each times Vi; divide by ONE projected ROI area. No extra factor 2.','m':'Fixed full original m by atom ID; do not use alive masks','omega':'1/r0; additional consistent delta scaling diagnostic of same main state only','taper':'exact saved original h, independently checked against source; no toggled-taper energy variants','ideal_extension':'cross: s00*r0, noncross: 0, all original bonds retained in deviatoric sum; not global coordinates','sensitivity':'ONE neighboring gap z=dx, same grid phase. Position/finite-height check only; not phase or mesh convergence.','limits':{'formal_main_integrals':1,'formal_control_integrals':1,'new_LAMMPS_instances':0,'equilibrium_solves':0,'CFD_microbubble_AI_runs':0,'material_or_grid_scans':0,'integration_wall_s':60,'process_peak_RSS_MiB':1024,'threads':1},'gates':{'relative_two_form_global_error':1e-12,'relative_pointwise_error_scaled_by_max_W':1e-12,'m_reference_relative_error':1e-13,'lambda_allowance':1e-14,'volume_and_area':'Exact frozen definitions','theory_discrete_agreement':'NO arbitrary percentage gate; report actual difference','geometric_translation_relative_error':1e-12},'stop':'Failure of identity, periodic/taper or energy definition stops certification; never change material/area to force agreement'}
dump(O/'DESIGN_FREEZE.json',freeze)
dump(O/'PREPARATION_RESOURCES.json',{'wall_s':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'formal_plane_integrals':0,'native_instances':0,'existing_independent_energy_evaluations_of_saved_state':1})
print(json.dumps({'status':'GATES_PASSED_AND_DESIGN_FROZEN','G0_theory':form1,'theory_K_fraction':kterm/form1,'reference_bonds':nb,'saved_energy_error':err,'frozen_sha256':sha(O/'DESIGN_FREEZE.json')}))
