from analysis.dual_track_v1.common import *
from analysis.dual_track_v1.low_a import run
import traceback,time
assert json.loads((A/'PERIODIC_TOPOLOGY_AUDIT.json').read_text())['status']=='PASS'
d=json.loads((A/'DESIGN_FREEZE.json').read_text());c=d['criteria'];op=d['optimizer'];opt={k:op[k] for k in ['maxiter','maxfun','maxls','maxcor','ftol','gtol']}
freeze=dict(amplitude_m=2e-6,budget=dict(per_state_max_function_evaluations=60000),optimizer=opt,uncertainty='Empirical residual work plus last10 plateau and floating roundoff; not rigorous truncation error or material uncertainty',criteria=dict(relative_free_force_L2_over_controlled_force_L2=c['relative_free_force'],max_free_force_N=c['max_free_force_N'],energy_span_last10_relative=c['energy_plateau_relative'],max_position_change_last10_m=c['position_plateau_m'],boundary_reaction_imbalance=c['reaction_imbalance'],max_contact_energy_over_U=c['contact_energy_fraction'],max_contact_force_over_reaction=c['contact_force_fraction'],native_energy_force_relative=c['native_force_relative']))
assert not (A/'PVCE_DESIGN_FREEZE.json').exists();save(A/'PVCE_DESIGN_FREEZE.json',freeze);start=time.monotonic();results=[]
try:
 for i,tip in enumerate(d['tip_x_m']):
  case='A'+str(i);out=A/case;out.mkdir();save(out/'config.json',dict(tip_m=tip));run(case);r=json.loads((out/'EQUILIBRIUM.json').read_text());results.append(r)
  if not r['accepted']:break
 save(A/'PVCE_EQUILIBRIUM_AUDIT.json',dict(accepted_three=len(results)==3 and all(r['accepted'] for r in results),formal_states=len(results),states=results,wall_s=time.monotonic()-start))
 if len(results)==3 and all(r['accepted'] for r in results):
  u=np.array([r['U_J'] for r in results]);ep=np.array([r['empirical_uncertainty_J'] for r in results]);area=d['delta_A_m2'];g2=(u[0]-u[1])/area;g3=(3*u[0]-4*u[1]+u[2])/(2*area);e2=(ep[0]+ep[1])/area;e3=np.dot([3,4,1],ep)/(2*area)
  save(A/'PVCE_ENERGY_RELEASE_RESULTS.json',dict(status='PASS' if min(g2-e2,g3-e3)>0 else 'ENERGY_DIFFERENCE_BELOW_ERROR_BUDGET',identity=d['G_identity'],U_J=u,empirical_U_uncertainty_J=ep,G_two_J_m2=g2,G_three_J_m2=g3,empirical_error_G_two_J_m2=e2,empirical_error_G_three_J_m2=e3,delta_A_m2=area,scope=freeze['uncertainty']))
except BaseException as e:
 save(A/'LOW_STOP.json',dict(case=case,error=repr(e),traceback=traceback.format_exc(),wall_s=time.monotonic()-start));raise
