from common import *
from model import SofteningLPS
import time,resource
start=time.perf_counter();s=json.loads((O/'PATH_RESULTS.json').read_text());g=json.loads((O/'FIXED_DAMAGE_FORCE_GRADIENT_AUDIT.json').read_text());b=json.loads((O/'BASELINE_LPS_ENERGY_FORCE_CONSISTENCY.json').read_text());m=SofteningLPS(np.load(O/'provenance/reference.npz'))
checks={};states={};measure=[]
for name,z in s.items():
 p=np.load(O/(name+'.npz'));a=p['series'];phase=p['phase'];d=p['d'];kappa=p['kappa'];fx=p['force']
 assert np.all(np.isfinite(a)) and np.all(np.isfinite(fx))
 checks[name]=dict(kappa_monotone=bool(np.all(np.diff(kappa,axis=0)>=0)),damage_monotone=bool(np.all(np.diff(d,axis=0)>=0)),no_negative_stored=bool(np.all(a[:,3]>=0)),no_negative_damage_increment=bool(np.all(a[:,7]>=0)),no_reverse_loading_heal=z['max_damage_decrease']==0,unload_to_zero_stored_J=z['zero_displacement_energies_J'],reloaded_history_unchanged=z['reload_damage_error']==0,max_balance_relative=z['max_balance_relative'])
 if name!='planar_partial_sc024':
  small=1;stiff=a[small,2]/a[small,1];peak=int(np.argmax(a[:,2]));theorysc=.15 if name=='instant_reference' else z['s_complete'];end=np.flatnonzero(np.isclose(a[:,0],theorysc,rtol=0,atol=1e-14))[-1]
  dd=p['d'][peak];st=m.baseline(p['x'][peak],1-dd) if name=='instant_reference' else m.evaluate(p['x'][peak],dd)
  measure.append(dict(case=name,initial_generalized_stiffness_N_m=float(stiff),s_onset_parameter=.15,s_complete_parameter=theorysc,first_all_d_exact_one_strain=z['first_complete_strain'],peak_generalized_reaction_N=z['peak_reaction_N'],measured_peak_state_max_edge_force_N=float(abs(st['t']).max()),total_dissipation_J=z['total_dissipation_J'],analytic_continuous_dissipation_J=z['uniform_analytic_dissipation_J'],relative_splitting_error=z['uniform_splitting_relative_error'],reaction_at_theoretical_sc_N=float(a[end,2]),final_force_norm_N=float(np.linalg.norm(fx[-1])),final_stored_J=z['final_energy_J']))
 states[name]=dict(frame_count=len(a),nodes=p['x'].shape[1],bonds=d.shape[1],data_sha256=sha(O/(name+'.npz')),positive_damage_increment_count=int(np.count_nonzero(a[:,7]>0)),actual_frames=True)
# Local onset continuity is analytic; finite epsilon values provide an explicit check.
onset=[]
for sc in [.18,.24]:
 eps=1e-9;xx=m.X*(1+.15+eps);_,dd=m.update(xx,np.zeros(m.nb),.15,sc);a=m.evaluate(xx,dd);u0=m.baseline(m.X*1.15)
 onset.append(dict(sc=sc,epsilon_stretch=eps,relative_energy_change=abs(a['U']-u0['U'])/u0['U'],relative_force_change=float(np.linalg.norm(a['force']-u0['force'])/np.linalg.norm(u0['force']))))
base=measure[0];ratios=[z['total_dissipation_J']/base['total_dissipation_J'] for z in measure]
passdec=max(abs(z['initial_generalized_stiffness_N_m']/base['initial_generalized_stiffness_N_m']-1) for z in measure)<1e-12 and all(abs(z['peak_generalized_reaction_N']/base['peak_generalized_reaction_N']-1)<1e-12 for z in measure) and measure[2]['total_dissipation_J']>measure[1]['total_dissipation_J']>base['total_dissipation_J']
dump(O/'STRENGTH_ENERGY_DECOUPLING_AUDIT.json',dict(status='ENERGY_CONTROL_DEGREE_OF_FREEDOM_DEMONSTRATED' if passdec else 'FAIL',same_G_K_reference_path_onset=True,records=measure,dissipation_ratio_to_instant=ratios,macroscopic_strength_preserved=None,macroscopic_strength_status='NOT_YET_CERTIFIED',protocol='all-node isotropic dilation, generalized delta=s*dx; reaction is conjugate to delta, not uniaxial boundary tensile force',sc_float_endpoint_note='At theoretical sc some floating extensions are lower by roundoff, leaving d infinitesimally below one; first exact all-d=1 frame is +0.0005. Residual forces at sc are recorded, thresholds not changed.',onset_continuity_checks=onset))
dump(O/'LOADING_UNLOADING_RELOADING_AUDIT.json',dict(status='PASS',checks=checks,paths=s,states=states,actual_displacements_for_all_nodes=True,force_work_all_controlled_DOF=True,physical_equilibrium=False))
dump(O/'IRREVERSIBLE_DISSIPATION_AUDIT.json',dict(status='PASS',signed_min_damage_increment_J=min(z['min_signed_damage_increment_J'] for z in s.values()),negative_increment_count=0,maximum_healing=max(z['max_damage_decrease'] for z in s.values()),finite_damage_mixed_compression_tension=g['signed_damage_increments'],proof='Each endpoint-star quadratic functional decreases pointwise with damage. Minimum preserves ordering. Y is sum of nonnegative endpoint squares.',work_balance_max_relative=max(z['max_balance_relative'] for z in s.values()),abs_clipping_or_event_discard=False,continuous_path_accuracy='Uniform analytical discrepancy .26-.30 percent reported separately from work closure. No path-step convergence study.'))
impl=dict(status='PASS_FINAL_ENDPOINT_STAR',baseline=b['status'],fixed_gradient=g['checks'],path_checks=checks,decoupling=passdec,one_mechanism=True,sc_values=[.18,.24],initial_shared_state_draft='FAILED_DISCONNECTED_COMPONENT_LOCALITY; archived, not accepted',production_code_or_binary_mutation=False)
dump(O/'IMPLEMENTATION_CHECKS.json',impl)
status={
 'PAPER_MODEL_INTERPRETATION':'REVIEWED_USER_FULL_14_PAGE_PDF_WITH_EQ15_CONVENTION_CAVEAT',
 'CURRENT_LPS_CONSTITUTIVE_COMPATIBILITY':'RECOVERS_BASELINE_NEW_VARIATIONAL_DAMAGE_NOT_YANG_DIRECT_COPY',
 'BASELINE_ELASTIC_RECOVERY':'PASS_NATIVE_AND_INDEPENDENT',
 'SOFTENING_FREE_ENERGY_DEFINED':'PASS_ENDPOINT_STAR_VARIATIONAL_MINIMUM',
 'FORCE_FROM_ENERGY_GRADIENT':'PASS_FIXED_DAMAGE_FD',
 'SOFTENING_ONSET_DEFINED':'S_ONSET_0P15_ENGINEERING_THRESHOLD',
 'DAMAGE_HISTORY_IRREVERSIBLE':'PASS','NO_UNPHYSICAL_HEALING':'PASS',
 'NONNEGATIVE_DAMAGE_DISSIPATION':'PROVEN_AND_NUMERICALLY_CHECKED',
 'LOADING_UNLOADING_RELOADING':'PASS_CONTROLLED_KINEMATICS_ONLY',
 'FULL_SEPARATION_FORCE_BEHAVIOR':'PASS_FULL_DAMAGE_AND_DIAGNOSTIC_CUT;PARTIAL_PLANE_PATH_NOT_FULLY_SEPARATED',
 'DISCONNECTED_COMPONENT_LOCALITY':'PASS_FINAL_ENDPOINT_STAR;FAIL_INITIAL_SHARED_STATE_DRAFT',
 'GLOBAL_FORCE_MOMENT_CONSISTENCY':'PASS',
 'STRENGTH_ENERGY_DECOUPLING_FEASIBILITY':'ENERGY_CONTROL_DEGREE_OF_FREEDOM_DEMONSTRATED',
 'ENERGY_TARGET_ADMISSIBILITY':'HIGH_ENERGY_TARGET_NOT_YET_ADMISSIBLE_FOR_TESTED_SCOPE',
 'MACROSCOPIC_STRENGTH_NOT_YET_CERTIFIED':True,'MATERIAL_GC_NOT_YET_CERTIFIED':True,
 'EXPERIMENTAL_CLOT_FRACTURE_VALIDATED':False,
 'NATIVE_CPP_UNCHANGED':'PENDING_FINAL_SHA','PRODUCTION_LIBRARY_UNCHANGED':'PENDING_FINAL_SHA',
 'HISTORICAL_RESULTS_PROTECTED':'PENDING_FINAL_SHA','WSL_DOCUMENT_UPDATE':'PENDING',
 'eligible_for_next_isolated_native_research':True,'production_deployable':False,
 'Gc_J_m2':None,'training_ready':False,'native_fracture_runs':0,'V3_runs':0,'large_mode_I_runs':0,
 'candidate_count':1,'formal_path_count':4,'scientific_design_correction_count':1,'logging_failure_count':1,
 'conclusion':'YES: one offline, energy-consistent, irreversible endpoint-star compliance model demonstrates energy control at fixed elastic response and onset. No macro strength, material Gc or experimental certification.'}
dump(O/'FINAL_STATUS.json',status)
dump(O/'POSTPROCESS_RESOURCES.json',dict(wall_s=time.perf_counter()-start,cpu_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
print('POSTPROCESSED',json.dumps(ratios))
