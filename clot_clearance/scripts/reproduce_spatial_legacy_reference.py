"""Re-run the unchanged 27-point branch tracker in a new output directory."""
from pathlib import Path
import json
import hashlib
import numpy as np
from pd_clot.pfpd_candidate import branch_stability_tracking as tracker
from pd_clot.pfpd_candidate.coupling_stability_audit import read_json, read_npz, save_json
from scripts import run_branch_stability_tracking as runner

ROOT = Path(__file__).resolve().parents[1]
TASK = ROOT/'verification/pfpd_at1_spatial_discretization_audit_v1'
OLD = ROOT/'verification/pfpd_at1_branch_stability_tracking_v1'
NEW = TASK/'level_0/legacy_reproduction'


def main():
    NEW.mkdir(exist_ok=False)
    (NEW/'provenance').mkdir()
    config = read_json(OLD/'provenance/PREDECLARED_CONFIG.json')
    config['schedules']['A'] = [.009, .0095, .00975, .009875, .0099375,
        .00996875, .009984375, .0099921875, .01, .010125, .01025]
    config['scope'] = 'Exact old geometry/volume; new requested strain schedule; no refined geometry'
    save_json(NEW/'provenance/PREDECLARED_CONFIG.json', config)
    (NEW/'provenance/PREDECLARED_CONFIG.sha256').write_text(hashlib.sha256(
        (NEW/'provenance/PREDECLARED_CONFIG.json').read_bytes()).hexdigest()+'\n')
    # Route writers only. Imported algorithms and original settings/tolerances
    # are unchanged; this process owns these globals and restores them below.
    old_tracker_out, old_runner_out = tracker.OUT, runner.OUT
    try:
        tracker.OUT = NEW
        runner.OUT = NEW
        result = runner.tracking_schedule('A')
    finally:
        tracker.OUT, runner.OUT = old_tracker_out, old_runner_out
    original = read_json(OLD/'critical_bracket/SCHEDULE_A.json')
    current = read_json(NEW/'critical_bracket/SCHEDULE_A.json')
    comparisons = {}
    for endpoint in ['lower','upper']:
        a, b = original[endpoint], current[endpoint]
        az = read_npz(OLD/'branch_states'/a['state']/'STATE.npz')
        bz = read_npz(NEW/'branch_states'/b['state']/'STATE.npz')
        av = read_npz(OLD/'hessian'/a['state']/'CONE_AND_MODE.npz')['critical_mode']
        bv = read_npz(NEW/'hessian'/b['state']/'CONE_AND_MODE.npz')['critical_mode']
        comparisons[endpoint] = dict(
            maximum_x_error_m=float(abs(az['x']-bz['x']).max()),
            maximum_phase_error=float(abs(az['s']-bz['s']).max()),
            reaction_relative_error=abs(a['reaction_N']-b['reaction_N'])/abs(a['reaction_N']),
            energy_relative_error=abs(a['metrics']['total_J']-b['metrics']['total_J'])/abs(a['metrics']['total_J']),
            force_residual_error=abs(a['metrics']['force_residual']-b['metrics']['force_residual']),
            KKT_residual_error=abs(a['metrics']['projected_phase_KKT']-b['metrics']['projected_phase_KKT']),
            curvature_error=abs(a['stability_indicator']-b['stability_indicator']),
            witness_absolute_cosine=float(abs(av@bv)/(np.linalg.norm(av)*np.linalg.norm(bv))),
            active_counts_equal=a['active']['counts']==b['active']['counts'],
            classification_equal=a['classification']==b['classification'])
    passed = result['bracket_status']=='BRACKETED' and current['final_interval']==original['final_interval']
    for row in comparisons.values():
        passed &= (row['maximum_x_error_m']<=1e-12 and row['maximum_phase_error']<=1e-8
                   and row['reaction_relative_error']<=1e-8 and row['energy_relative_error']<=1e-9
                   and row['force_residual_error']<=1e-9 and row['KKT_residual_error']<=1e-8
                   and row['curvature_error']<=1e-6 and row['witness_absolute_cosine']>=1-1e-6
                   and row['active_counts_equal'] and row['classification_equal'])
    report = dict(status='VERIFIED' if passed else 'COARSE_REFERENCE_REPRODUCTION_FAILED',
                  exact_original_geometry=True, original_interval=original['final_interval'],
                  reproduced_interval=current['final_interval'], width=current['width'],
                  method=current['method'], comparisons=comparisons,
                  scope='Legacy reference only; no physical-domain choice implied')
    save_json(TASK/'level_0/COARSE_REFERENCE_REPRODUCTION.json', report)
    print(json.dumps(report, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
