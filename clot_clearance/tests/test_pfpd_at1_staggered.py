from types import SimpleNamespace
import copy
import json
import numpy as np
from scripts.verify_at1_evolution import ROOT, OUT, config
from pd_clot.pfpd_candidate import staggered_at1 as solver

RESULTS = ROOT/'results/pfpd_at1_evolution_v1'


def test_saved_coupled_mechanics_phase_and_energy_all_converged():
    c = config()
    rows = json.loads((RESULTS/'loading_unloading/LOAD_PATH.json').read_text())
    assert len(rows) == len(c['schedules']['loading_unloading'])
    for row in rows:
        assert row['relative_force_residual'] <= c['mechanical_equilibrium']['relative_force_residual_tolerance']
        assert row['maximum_delta_s'] <= c['tolerances']['phase_step']
        assert row['relative_stagger_energy_change'] <= c['tolerances']['stagger_energy_relative']


def test_coupled_no_bond_heals_on_unload_and_reload():
    rows = json.loads((RESULTS/'loading_unloading/LOAD_PATH.json').read_text())
    assert len(rows) == len(config()['schedules']['loading_unloading']), 'Unload/reload was not reached: coupled acceptance gate FAILED'
    old = None
    for file in sorted((RESULTS/'loading_unloading').glob('load_*/accepted_state.npz')):
        with np.load(file) as state:
            s = state['s_pf_bond'].copy()
        if old is not None: assert np.all(s >= old)
        old = s
    assert rows[-1]['maximum_s'] > rows[10]['maximum_s']
    assert abs(rows[10]['maximum_s']-rows[12]['maximum_s']) <= 1e-10


def test_high_damage_stop_is_expected_and_does_not_accept_or_clip_proposal():
    summary = json.loads((RESULTS/'tensile_cube/SUMMARY.json').read_text())
    assert summary['status'] == 'HIGH_DAMAGE_GATE_REACHED'
    event = summary['event']; step = RESULTS/'tensile_cube'/f"load_{event['load_step']:03d}"
    assert event['proposed_maximum_s'] > config()['phase']['sc']
    assert event['last_accepted_maximum_s'] < config()['phase']['sc']
    assert not event['proposed_state_accepted'] and not event['continued_after_gate']
    assert not (step/'ACCEPTED.json').exists() and not (step/'FAILED.json').exists()
    with np.load(step/'high_damage_state.npz') as state:
        assert state['s_new'].max() > .95 and state['s_last_accepted'].max() < .95


def test_extra_stagger_iterations_change_no_accepted_state_materially():
    rows = json.loads((OUT/'staggered/EXTRA_ITERATIONS.json').read_text())
    c = config(); assert len(rows) == len(c['schedules']['loading_unloading'])
    assert max(r['max_delta_x_m'] for r in rows) <= c['tolerances']['extra_iteration_x_absolute_m']
    assert max(r['max_delta_s'] for r in rows) <= c['tolerances']['extra_iteration_state_absolute_s']


def test_mechanical_failure_prevents_phase_call_and_preserves_failed_state(monkeypatch,tmp_path):
    calls = []
    def failed(cloud,family,material,top,x,s,strain,c):
        return SimpleNamespace(x=x.copy(),converged=False,summary={'iterations':1},checks=[],history=[])
    monkeypatch.setattr(solver,'equilibrate_tensile',failed)
    monkeypatch.setattr(solver,'at1_constrained_update',lambda *a,**k:calls.append('phase'))
    target = tmp_path/'failure'
    summary = solver.run_path(config(),[0.],target)
    assert summary['status'] == 'FAILED' and calls == []
    assert (target/'FAILED.json').exists()
    with np.load(target/'load_000/failed_state.npz') as state:
        assert not np.any(state['s_current']) and not np.any(state['s_last_accepted'])


def test_stagger_limit_fails_without_accepting_or_continuing(tmp_path):
    c = config(); c['staggered']['maximum_iterations'] = 1
    summary = solver.run_path(c,[0.,.01],tmp_path/'limit')
    assert summary['status'] == 'FAILED' and summary['accepted_steps'] == 0
    assert not (tmp_path/'limit/load_001').exists()
    assert 'MAXIMUM_ITERATIONS' in summary['event']['reason']


def test_rank_loss_fails_loudly_and_saves_old_phase(monkeypatch,tmp_path):
    def rankloss(*a,**k): raise ValueError('Rank deficient family; outside candidate scope')
    monkeypatch.setattr(solver,'equilibrate_tensile',rankloss)
    summary = solver.run_path(config(),[0.,.01],tmp_path/'rank')
    assert summary['status'] == 'FAILED' and summary['accepted_steps'] == 0
    assert (tmp_path/'rank/load_000/failed_state.npz').exists()
    assert not (tmp_path/'rank/load_001').exists()


def test_external_work_is_raw_trapezoidal_diagnostic():
    rows = json.loads((RESULTS/'loading_unloading/LOAD_PATH.json').read_text())
    work = 0.
    for index,row in enumerate(rows):
        if index:
            prev = rows[index-1]
            work += .5*(row['reaction_N']+prev['reaction_N'])*(row['displacement_m']-prev['displacement_m'])
        assert work == row['external_work_J']
        assert row['energy_balance_residual_J'] == row['total_J']-work
