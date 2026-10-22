"""Post-discovery review: do not retune actions or relabel original final tests."""
import time,numpy as np,h5py
from scripts.straight_tube_exp1.f10ab.metrics import geometry
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import context
from .common import *
from .run import batch
def main():
    begin=time.perf_counter();rows=read(BASE/'ALL_RUNS.json');by={r['id']:r for r in rows};v=read(BASE/'VALIDATION.json');observations=[]
    for r in rows:
        if r['status']!='PASS' or not r['targets']['raw_extension_above_floor']:continue
        if r['id'].endswith(('_DT05','_FINE')):continue
        with h5py.File(BASE/r['trajectory_path']) as h:x0=h['positions_m'][0];x=h['positions_m'][-1]
        g=geometry(x);labels=g['component_labels'];ids=np.flatnonzero(labels==np.argmax(np.bincount(labels)));a=geometry(x0[ids]);b=geometry(x[ids]);rep=next((q for gg in v['groups'] for q in gg['replicates'] if q['id']==r['id']),None)
        observations.append(dict(id=r['id'],global_axis_delta_m=r['targets']['axis_delta_m'],global_components=g['proximity_components'],largest_component_ids=ids,largest_same_ID_axis_delta_m=b['principal_rms_lengths_3d_m']-a['principal_rms_lengths_3d_m'],largest_same_ID_Rg_change_m=b['Rg_3d_m']-a['Rg_3d_m'],orientation_valid=g['orientation_valid_3d'],numerically_resolved=bool(rep and rep['numerical_pass'] and np.any(np.array(rep['targets']['axis_delta_m'])>np.maximum(1e-9,5*np.array(rep['axis_noise_m'])))),axis_noise_m=rep['axis_noise_m'] if rep else None,scope='Post-hoc observation, no action retuning; final test observations have no new refinement and remain unconfirmed'))
    write(BASE/'SEED_SPECIFIC_EXTENSION.json',dict(rows=observations,discovery_r0_absolute_extension_found=False,cross_state_target_T3_confirmed=False,identity=identity()))
    # Exactly one new model-domain diagnostic. The original final test remains failed.
    old=by['final_T4_1'];t=old['task'].copy();t.update(id='joined_followup_final_T4_1',joined_diagnostic=True)
    with h5py.File(BASE/old['trajectory_path']) as h:x=h['positions_m'][-1]
    c=context();u=c['basis'].vector(t['amplitudes'],t['phases']);p=u@c['basis'].query_basis(x)[0]
    write(BASE/'LOW_DOMAIN_FOLLOWUP_PLAN.json',dict(task=t,original_failed_run=old['id'],original_failure=old['reasons'],reconstructed_last_accepted_pressure_min_Pa=float(abs(p).min()),reconstruction_scope='Independent CPU reconstruction, not original failed GPU kernel output',extension_gate_sha256=sha(BASE/'LOW_PRESSURE_ADMISSION.json'),purpose='Domain-extension diagnostic only; no optimizer feedback and no relabeling of frozen final-test results',max_new_physics_calls=1,identity=identity()))
    rr=batch([t],1);write(BASE/'LOW_DOMAIN_FOLLOWUP.json',dict(rows=rr,original_final_test_preserved=True,wall_s=time.perf_counter()-begin,identity=identity()))
    from . import audit,interpret,diagnostics
    audit.main();interpret.main()
    if not (BASE/'ROBUSTNESS_FIELD_FORCE.json').exists():diagnostics.main()
if __name__=='__main__':main()
