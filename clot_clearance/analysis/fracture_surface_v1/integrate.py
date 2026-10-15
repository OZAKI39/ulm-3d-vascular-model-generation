"""Exactly one frozen main plane and one neighboring-gap control; no native calls."""
from common import *
from datetime import datetime, timezone
import time, resource, os

wall=time.perf_counter(); cpu=time.process_time()
assert not (O/'DISCRETE_SURFACE_ENERGY_RESULTS.json').exists(), 'Formal result already exists'
freeze=json.loads((O/'DESIGN_FREEZE.json').read_text()); freeze_sha=sha(O/'DESIGN_FREEZE.json')
assert json.loads((O/'PRISTINE_REFERENCE_TOPOLOGY_AUDIT.json').read_text())['status']=='PASS'
assert json.loads((O/'SOURCE_LPS_CONSISTENCY_AUDIT.json').read_text())['status']=='PASS_FOR_INTACT_REFERENCE_IDEAL_POINT_STATE'
assert sha(W/'reference_topology.npz')==freeze['reference_sha256']
read=time.perf_counter(); raw=np.load(W/'reference_topology.npz')
X,b,r,h,m,image=[raw[k] for k in ['X','bonds','r0_m','h','m','periodic_image_y']]; read_s=time.perf_counter()-read
n=len(X); A=freeze['ROI']['A_m2']; roi=(X[:,0]>=freeze['ROI']['x_half_open_m'][0])&(X[:,0]<freeze['ROI']['x_half_open_m'][1]); nodeV=np.full(n,V)
assert int(roi.sum())==freeze['ROI']['nodes'] and len(np.unique(X[roi,0]))==30
assert A==30*DX*PERIOD
degree=np.bincount(b.ravel(),minlength=n); theory=json.loads((O/'THEORY_FORMULA_AUDIT.json').read_text())
arrays={'atom_id':np.arange(n)+1,'X_m':X,'V_m3':nodeV,'m_full_m4':m,'roi':roi,'planes_m':np.array(freeze['planes_m'])}; audits=[]; rows=[]; identities=[]; layers=[]; times=[]
integral_start=time.perf_counter()
for index,z0 in enumerate(freeze['planes_m']):
    t=time.perf_counter(); side=X[:,2]<z0; cross=side[b[:,0]]!=side[b[:,1]]
    assert np.min(abs(X[:,2]-z0))>24e-6
    if index==0:assert np.array_equal(cross,raw['cross_z0'])
    # Form 1: reference-volume fraction. Endpoint sum includes both sides once.
    mc=np.bincount(b[cross].ravel(),weights=np.repeat(r[cross]*V*h[cross],2),minlength=n)
    lam=mc/m; assert lam.min()>=-1e-14 and lam.max()<=1+1e-14
    WK=4.5*K*S00**2*lam**2; WG=7.5*G*S00**2*(lam-lam**2)
    Wlambda=.5*S00**2*((9*K-15*G)*lam**2+15*G*lam)
    # Form 2: independently form e, theta and every directed deviatoric extension.
    # Zero-e noncross bonds MUST remain in this sum because theta generally is nonzero.
    e=np.zeros(len(b)); e[cross]=S00*r[cross]; thsum=np.zeros(n)
    for end in [0,1]:thsum+=np.bincount(b[:,end],weights=h*e*V,minlength=n)
    theta=3*thsum/m; dev=np.zeros(n); dev_non_cross=np.zeros(n)
    for end in [0,1]:
        ii=b[:,end]; ed=e-theta[ii]*r/3; q=h/r*V*ed**2
        dev+=np.bincount(ii,weights=q,minlength=n)
        dev_non_cross+=np.bincount(ii[~cross],weights=q[~cross],minlength=n)
    W2K=.5*K*theta**2; W2G=7.5*G*dev/m; W2=W2K+W2G
    U1=float(np.sum(Wlambda[roi]*nodeV[roi])); U2=float(np.sum(W2[roi]*nodeV[roi])); value=U1/A
    point_error=float(np.max(abs(Wlambda-W2))/np.max(abs(Wlambda))); global_error=abs(U1-U2)/U1
    split_error=max(float(np.max(abs(WK-W2K))/np.max(WK)),float(np.max(abs(WG-W2G))/np.max(WG)))
    assert global_error<1e-12 and point_error<1e-12 and split_error<1e-12
    active=roi&(mc>0); assert np.all(degree[active]==122)
    # Every active node and its entire reference sphere stay within physical X/Z bounds.
    clearance=np.minimum(.0025-abs(X[active,0]),.0005-abs(X[active,2])); assert clearance.min()>DELTA
    lower=float(np.sum(Wlambda[roi&side]*V)/A); upper=float(np.sum(Wlambda[roi&~side]*V)/A)
    # Scaling omega and m consistently is a same-state algebra diagnostic, not a third geometry.
    scale_error=float(np.max(abs((mc*DELTA)/(m*DELTA)-lam)))
    assert scale_error<1e-14
    incident=cross&(roi[b[:,0]]|roi[b[:,1]]); both=cross&roi[b[:,0]]&roi[b[:,1]]
    bp=.5*float(np.sum(lam[roi]**2*V))/A; beta=.5*float(np.sum(lam[roi]*V))/A
    kvalue=float(np.sum(WK[roi]*V)/A); gvalue=float(np.sum(WG[roi]*V)/A)
    row={'name':['main','neighbor_gap_position_control'][index],'plane_z_m':z0,'A_m2':A,'U_ideal_ROI_J':U1,'G0_discrete_J_m2':value,'K_contribution_J_m2':kvalue,'G_contribution_J_m2':gvalue,'K_fraction':kvalue/value,'G_fraction':gvalue/value,'G0_over_theory':value/theory['G0_theory_J_m2'],'difference_from_theory_J_m2':value-theory['G0_theory_J_m2'],'relative_difference_from_theory':value/theory['G0_theory_J_m2']-1,'beta_discrete_m':beta,'beta_prime_discrete_m':bp,'beta_difference_m':beta-theory['beta_m'],'beta_prime_difference_m':bp-theory['beta_prime_m'],'K_term_difference_from_theory_J_m2':kvalue-theory['K_contribution_J_m2'],'G_term_difference_from_theory_J_m2':gvalue-theory['G_contribution_J_m2'],'lower_side_G_J_m2':lower,'upper_side_G_J_m2':upper,'noncross_bond_deviatoric_contribution_J_m2':float(np.sum((7.5*G/m*dev_non_cross)[roi]*V)/A)}
    rows.append(row)
    audits.append({'plane_z_m':z0,'unique_reference_bonds':len(b),'cross_plane_unique_bonds_whole_specimen':int(cross.sum()),'cross_plane_unique_bonds_incident_to_ROI':int(incident.sum()),'cross_plane_unique_bonds_both_endpoints_in_ROI':int(both.sum()),'cross_plane_directed_endpoint_contributions_in_ROI':int(roi[b[cross]].sum()),'cross_plane_Y_seam_unique_bonds':int(np.count_nonzero(cross&(image!=0))),'cross_plane_Y_seam_unique_bonds_incident_to_ROI':int(np.count_nonzero(incident&(image!=0))),'ROI_nodes_total':int(roi.sum()),'ROI_nodes_nonzero_energy':int(active.sum()),'ROI_lower_active':int((active&side).sum()),'ROI_upper_active':int((active&~side).sum()),'all_specimen_active_nodes':int((mc>0).sum()),'on_plane_nodes':int((X[:,2]==z0).sum()),'lambda_all_range':[lam.min(),lam.max()],'lambda_ROI_active_range':[lam[active].min(),lam[active].max()],'m_full_ROI_active_range_m4':[m[active].min(),m[active].max()],'m_cross_ROI_active_range_m4':[mc[active].min(),mc[active].max()],'all_active_ROI_degrees':np.unique(degree[active]),'minimum_active_ROI_center_boundary_distance_m':clearance.min(),'A_m2':A,'reference_sha256':freeze['reference_sha256'],'unique_periodic_images':True,'extra_factor_two':False})
    identities.append({'plane_z_m':z0,'form1_U_J':U1,'form2_U_J':U2,'absolute_U_difference_J':U1-U2,'relative_global_error':global_error,'pointwise_error_normalized_max_W':point_error,'K_G_separate_pointwise_error':split_error,'theta_error_vs_3lambda_s00':float(np.max(abs(theta-3*lam*S00))),'constant_omega_rescaling_lambda_abs_error':scale_error,'all_non_cross_deviatoric_terms_included':True,'reference_family_unchanged':True,'native_calls':0,'kinematically_realizable_displacement_asserted':False})
    for zz in np.unique(X[:,2]):
        mask=roi&(X[:,2]==zz); qdist=abs(zz-z0)/DELTA; cont=.5-2*qdist/3+qdist**4/6 if qdist<=1 else 0.
        layers.append({'plane_index':index,'z_m':float(zz),'distance_m':float(abs(zz-z0)),'node_count':int(mask.sum()),'lambda_mean':float(lam[mask].mean()),'lambda_min':float(lam[mask].min()),'lambda_max':float(lam[mask].max()),'continuum_cap_lambda':cont,'W_K_J_m3':float(WK[mask].mean()),'W_G_J_m3':float(WG[mask].mean()),'W_total_J_m3':float(Wlambda[mask].mean())})
    for key,arr in {'m_cross_m4':mc,'lambda':lam,'theta':theta,'W_K_J_m3':WK,'W_G_J_m3':WG,'W_ideal_J_m3':Wlambda,'W_point_state_J_m3':W2,'DeltaU_J':Wlambda*V,'active_ROI':active,'cross_plane_bonds':cross}.items():arrays[f'{index}_{key}']=arr
    times.append(time.perf_counter()-t)
integral_s=time.perf_counter()-integral_start
assert len(rows)==2 and integral_s<freeze['limits']['integration_wall_s']
# Plane translated by one Z layer. Only compare supported node pairs, preserving ID mapping.
active=np.where(arrays['0_active_ROI'])[0]; shifted=active+1
assert np.allclose(X[shifted],X[active]+[0,0,DX],rtol=0,atol=1e-18)
ldiff=float(np.max(abs(arrays['0_lambda'][active]-arrays['1_lambda'][shifted])))
gdiff=rows[1]['G0_discrete_J_m2']-rows[0]['G0_discrete_J_m2']
same_set=np.array_equal(arrays['0_cross_plane_bonds'],arrays['1_cross_plane_bonds'])
assert abs(gdiff)/rows[0]['G0_discrete_J_m2']<1e-12 and ldiff<1e-12 and not same_set
sensitivity={'status':'PASS_SAME_PHASE_INTERIOR_TRANSLATION','main_z_m':0.,'control_z_m':DX,'different_actual_reference_bond_ID_sets':not same_set,'same_untranslated_set':same_set,'translation_equivalent_local_stencils':True,'same_grid_phase':True,'grid_phase_test':False,'mesh_convergence_test':False,'cross_plane_counts':[v['cross_plane_unique_bonds_whole_specimen'] for v in audits],'lambda_mapped_max_abs_difference':ldiff,'G_difference_J_m2':gdiff,'G_relative_difference':gdiff/rows[0]['G0_discrete_J_m2'],'K_difference_J_m2':rows[1]['K_contribution_J_m2']-rows[0]['K_contribution_J_m2'],'G_shear_difference_J_m2':rows[1]['G_contribution_J_m2']-rows[0]['G_contribution_J_m2'],'conclusion':'The two interior gaps sample translation-equivalent intact neighborhoods. No resolved finite-height effect in this comparison; not independent evidence of phase/grid convergence.'}
# Actual shell weights for one central active node; no alternate-weight energy computation.
node=int(active[len(active)//2]); incident=(b==node).any(1); shell=np.rint((r/DX)**2).astype(int); shells=[]
for ss in np.unique(shell[incident]):
    mask=incident&(shell==ss); w=r[mask]*V*h[mask]; un=r[mask]*V
    shells.append({'r_squared_over_dx2':int(ss),'neighbors':int(mask.sum()),'r_over_dx_mean':float((r[mask]/DX).mean()),'taper_mean':float(h[mask].mean()),'native_weighted_m_contribution_m4':float(w.sum()),'unweighted_geometric_sum_m4':float(un.sum()),'fraction_of_original_m':float(w.sum()/m[node]),'partial_volume_loss_m4':float((un-w).sum())})
dump(O/'TAPER_SHELL_AUDIT.json',{'node_atom_ID':node+1,'position_m':X[node],'native_m_m4':m[node],'continuous_no_taper_m_m4':np.pi*DELTA**4,'shells':shells,'scope':'Shell geometric weights only; no second taper model, extra energy integral or sensitivity scan'})
dump(O/'CROSS_PLANE_REFERENCE_TOPOLOGY_AUDIT.json',{'status':'PASS','cases':audits,'frozen_design_sha256':freeze_sha,'ROI_interpretation':'Both endpoint point-state energies. Incident bonds outside the X window remain in each full node family; counts distinguish whole, incident and both-in-window sets.'})
dump(O/'DISCRETE_SURFACE_ENERGY_RESULTS.json',{'status':'ACCEPTED_NOMINAL_ONLY','calibrated_Gc':False,'theory_G0_J_m2':theory['G0_theory_J_m2'],'cases':rows,'delta_over_dx':DELTA/DX,'interpretation':'Ideal formal point-state energy allocation on the actual pristine 3D topology; no natural crack advance or work of separation measured.','scope_comparison':'NOMINAL_ENERGY_SCALE_BELOW_EXPERIMENTAL_REFERENCES','frozen_design_sha256':freeze_sha,'formula':'sum_ROI(.5*s00**2*((9*K-15*G)*lambda**2+15*G*lambda)*V)/A; both sides once'})
dump(O/'INDEPENDENT_ENERGY_FORM_CROSSCHECK.json',{'status':'PASS','cases':identities,'form2':'theta=(3/m)*sum_all omega*r*e*V*h; W=.5*K*theta^2+(15*G/2m)*sum_all omega*(e-theta*r/3)^2*V*h','directed_endpoint_sums':True,'pristine_family_identity_only':True})
dump(O/'GEOMETRIC_SENSITIVITY.json',sensitivity); dump(O/'LAYER_PROFILES.json',layers)
np.savez_compressed(O/'LAMBDA_AND_ENERGY_DENSITY.npz',**arrays)
peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
assert peak<freeze['limits']['process_peak_RSS_MiB'] and sha(O/'DESIGN_FREEZE.json')==freeze_sha
dump(O/'INTEGRATION_RESOURCES.json',{'completed_utc':datetime.now(timezone.utc).isoformat(),'input_read_s':read_s,'integration_s':integral_s,'case_integration_s':times,'process_wall_s':time.perf_counter()-wall,'process_CPU_s':time.process_time()-cpu,'peak_RSS_MiB':peak,'formal_plane_integrals':2,'native_calls':0,'equilibrium_solves':0,'threads':{'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS')},'GPU_used':False,'irreproducibility_observed':False})
print(json.dumps({'cases':rows,'sensitivity':sensitivity,'resources':json.loads((O/'INTEGRATION_RESOURCES.json').read_text())},indent=2))
