"""Recompute saved free residuals and fixed-x damage jumps without rerunning trajectories."""
import hashlib,json
from pathlib import Path
import numpy as np
from pd_clot.regularization.equilibrium_adr import make_fixture,FrozenForceEvaluator
from pd_clot.regularization.erosion_equilibrium_runner import write

ROOT=Path(__file__).resolve().parents[1]
V=ROOT/'verification/mechanical_equilibrium_v1/adr_completion'
R=ROOT/'results/streaming_erosion_equilibrium_v1'


def read(p):return json.loads(p.read_text())


def main():
    cases=[];paths=[V/'cases'/n for n in ['small_adr','small_adr_extra500','small_adr_resolved']]
    paths += [R/f'dt_erosion_{n}' for n in ['0p02','0p01','0p005']]
    for path in paths:
        c=read(path/'CONFIG.json');cl=make_fixture(c);states=sorted(path.glob('equilibrium*'));max_R=0.;rows=[]
        identity=read(path/'IDENTITY.json')
        assert identity['config_sha256']==hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()
        for p,sha in identity['source_sha256'].items():
            assert hashlib.sha256((path/'source_snapshot'/p).read_bytes()).hexdigest()==sha
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha
        for p in states:
            z=np.load(p/'state.npz');s=read(p/'SUMMARY.json');h=np.load(p/'history.npz')
            f=FrozenForceEvaluator(cl,z['integrity'],c)(z['x']);free=~cl.fixed
            residual=float(np.linalg.norm(f.residual[free])/max(np.linalg.norm(f.external_force_N[free]),c['mechanical_equilibrium']['force_floor_N']))
            np.testing.assert_array_equal(z['x'][cl.fixed],cl.X[cl.fixed])
            np.testing.assert_array_equal(z['bond_D'],1-z['integrity'])
            np.testing.assert_allclose(z['internal_force_N'],f.internal_force_N,rtol=0,atol=0)
            np.testing.assert_allclose(z['epsilon_plus'],f.epsilon_plus,rtol=0,atol=0)
            assert residual<=1e-4 and s['converged']
            np.testing.assert_array_equal(h['erosion_exposure_time_s'],float(z['erosion_exposure_time_s']))
            np.testing.assert_array_equal(h['represented_ultrasound_cycles'],0.)
            assert s['endpoint_mass_bound_ratio']<=1.
            assert all(ch['all_criteria_passed'] for ch in read(p/'checks.json')[-5:])
            max_R=max(max_R,residual);rows.append(dict(phase=s['phase'],exposure_s=s['erosion_exposure_time_s'],iterations=s['iterations'],relative_force_residual=residual))
        if (path/'damage_jumps.json').exists():
            s=read(path/'SUMMARY.json');assert s['status']=='COMPLETED_EQUILIBRATED_EROSION_DEVELOPMENT_RUN'
            assert s['accepted_damage_updates']==s['requested_damage_updates']
            for jp in sorted(path.glob('damage_jump_*')):
                z=np.load(jp/'state.npz');a=read(jp/'AUDIT.json')
                assert hashlib.sha256(z['x'].tobytes()).hexdigest()==a['coordinates_sha256']
                before=FrozenForceEvaluator(cl,1-z['D_before'],c)(z['x']).stored_energy_J
                after=FrozenForceEvaluator(cl,1-z['D_after'],c)(z['x']).stored_energy_J
                assert before==a['U_before_damage_J'] and after==a['U_after_damage_J']
                assert before-after==a['delta_U_damage_J']
                de=c['streaming_erosion']['erosion_exposure_step_s']
                np.testing.assert_array_equal(z['D_after'],np.minimum(1.,z['D_before']+de*z['rate_s_inv']))
                expected=c['streaming_erosion']['C_t_s_inv']*z['epsilon_plus_used']/c['streaming_erosion']['epsilon_ref']
                np.testing.assert_allclose(z['rate_s_inv'],expected,rtol=1e-15,atol=0)
            z=np.load(path/'states.npz');np.testing.assert_array_equal(z['bond_D'][0],0.)
            assert np.isclose(z['erosion_exposure_time_s'][-1],.1,rtol=0,atol=1e-15)
        cases.append(dict(name=path.name,equilibrium_states=len(states),maximum_relative_force_residual=max_R,checks=rows))
    write(V/'STATE_AUDIT.json',dict(status='PASS',cases=cases,scope='Recomputed existing kernel forces and signed damage energies; not an independent physical model'))

    metrics=['mean_damage','maximum_damage','maximum_bond_damage','maximum_epsilon_plus','maximum_displacement_m','stored_energy_J']
    histories={n:read(R/f'dt_erosion_{n}/history.json') for n in ['0p02','0p01','0p005']}
    reference=histories['0p005'];comparisons=[]
    for name,h in histories.items():
        by_metric={k:0. for k in metrics}
        for row in h:
            t=row['erosion_exposure_time_s'];ref=min(reference,key=lambda v:abs(v['erosion_exposure_time_s']-t))
            assert abs(ref['erosion_exposure_time_s']-t)<1e-14
            for k in metrics:
                delta=abs(row[k]-ref[k]);value=delta/abs(ref[k]) if ref[k] else (0. if delta==0 else float('inf'))
                by_metric[k]=max(by_metric[k],value)
        comparisons.append(dict(case=name,maximum_relative_differences_at_shared_exposures=by_metric,
            final=h[-1],equilibrium_iterations=sum(s['iterations'] for s in read(R/f'dt_erosion_{name}/SUMMARY.json')['equilibrations'])))
    worst=max(max(r['maximum_relative_differences_at_shared_exposures'].values()) for r in comparisons)
    tolerance=read(ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json')['mechanical_equilibrium']['slow_step_relative_acceptance']
    write(V/'SLOW_STEP_COMPARISON.json',dict(reference='0p005',method='Direct matching exposure times, no interpolation; max over all times and six predeclared metrics',
        acceptance=tolerance,maximum_relative_difference=worst,passed=worst<=tolerance,comparisons=comparisons))
    small=read(V/'SMALL_REFERENCE_GATE.json')['passed'];passed=small and worst<=tolerance
    write(V/'FINAL_STATUS.json',dict(mechanical_equilibrium='VERIFIED_FOR_CURRENT DEVELOPMENT_CASE' if passed else 'NOT VERIFIED',
        small_reference_gate=small,all_pre_relaxations_passed=True,all_damage_updates_equilibrium_gated=True,
        slow_step_gate=worst<=tolerance,physical_calibration='NOT PERFORMED',mesh_objectivity='STILL OPEN',
        wall_attachment='STILL OPEN',real_microbubble_CFD='NOT USED',erosion_energy='ENERGY_DISSIPATION_NOT_ESTABLISHED'))
    print(json.dumps(dict(state_audit='PASS',slow_step_maximum_relative_difference=worst,mechanical_pass=passed)))


if __name__=='__main__':main()
