"""Numerical energy calibration using cyclic opening of a symmetric weak plane.

All nodes are kinematically controlled in two rigid halves of a 0.5 mm cube.
Only plane-crossing bonds have the damageable weak-plane material. They fail by
the same cyclic energy increment used in the candidate, never by a cutting rule.
This coupon calibrates a DEMONSTRATION energy scale, not tissue toughness.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import brentq
from ..geometry import make_cloud
from ..fragment_mechanics import prepare_supported_shape, evaluate_supported
from .energy import potential, bond_quadrature_volume, EnergyLedger


def young_modulus(material):
    G, K = material['shear_modulus_Pa'], material['bulk_modulus_Pa']
    return 9*K*G/(3*K+G)


def cyclic_increment(driver_J, critical_J, cycles, coefficient, exponent, threshold=0.):
    ratio = np.maximum(np.asarray(driver_J)/critical_J-threshold, 0.)
    return cycles*coefficient*ratio**exponent


def coupon_cloud(base, spacing):
    cfg = copy.deepcopy(base['clot'])
    cells = int(round(.0005/spacing))
    if cells % 2 or not np.isclose(cells*spacing, .0005, rtol=0, atol=1e-12):
        raise ValueError('Coupon requires an even cell count over its fixed 0.5 mm dimensions')
    cfg.update(cells=[cells]*3, origin_m=[-.00025]*3, particle_spacing_m=spacing,
               fixed_base_thickness_m=spacing)
    return make_cloud(cfg)


def run_coupon(base, spacing, critical_density_Pa, increments=180, integrate_work=True):
    if critical_density_Pa <= 0:
        raise ValueError('Positive critical energy density required')
    cloud = coupon_cloud(base, spacing)
    pair = cloud.pairs
    side = np.where(cloud.X[:, 0] < 0, -.5, .5)
    weak = side[pair[:, 0]] != side[pair[:, 1]]
    B = len(pair)
    D = np.zeros(B)
    weights = bond_quadrature_volume(cloud)
    critical = critical_density_Pa*weights
    material, safety = base['material'], base['safety']
    dc = base['regularized_damage']
    ledger = EnergyLedger(0.)
    dissipation = np.zeros(B)
    history = []
    work = 0.
    # The opening ramp is a boundary protocol, identical in normalized strain.
    # It only supplies enough displacement for calibration, not a clot load.
    max_opening = 100*spacing*np.sqrt(2*critical_density_Pa/young_modulus(material))
    quadrature, qw = np.polynomial.legendre.leggauss(6)
    start = time.perf_counter()

    def evaluate(opening, g, shape):
        x = cloud.X.copy()
        x[:, 0] += opening*side
        mech = evaluate_supported(cloud, x, g, material, safety, shape)
        return x, mech, sum(potential(cloud, mech))

    def loading_work(opening, g, shape):
        val = 0.
        for q, w in zip(quadrature, qw):
            _, mech, _ = evaluate(opening*(q+1)/2, g, shape)
            val += w*float(-np.dot(mech[0][:, 0], side))*opening/2
        return val

    g = np.ones(B)
    shape = prepare_supported_shape(cloud, g, safety)
    # Continue the SAME opening-per-cycle ramp if a sensitivity case needs a
    # longer coupon protocol. Do not lower the failure threshold or cut bonds.
    for step in range(1, 3*increments+1):
        opening = max_opening*step/increments
        x, mech_before, before = evaluate(opening, g, shape)
        amplitude = .5*np.maximum(np.linalg.norm(x[pair[:, 1]]-x[pair[:, 0]], axis=1)/cloud.length-1, 0.)
        driver = .5*young_modulus(material)*amplitude**2*weights
        increment = cyclic_increment(driver, critical, 1, dc['C_E'], dc['m_E'], dc.get('threshold_ratio', 0.))
        oldD = D.copy()
        D[weak] = np.minimum(1., D[weak]+increment[weak])
        new_g = 1-D
        new_shape = prepare_supported_shape(cloud, new_g, safety)
        _, mech_after, after = evaluate(opening, new_g, new_shape)
        allocated = ledger.damage_drop(before, after, (D-oldD)*driver)
        dissipation += allocated
        # Independently integrate loading and unloading forces. No zero-energy
        # assumption is substituted for the force integral.
        if integrate_work:
            work += loading_work(opening, g, shape)-loading_work(opening, new_g, new_shape)
        elastic, stabilization = potential(cloud, mech_after)
        history.append(dict(coupon_cycle=step, opening_m=opening,
                            elastic_energy_at_peak_J=elastic, stabilization_energy_at_peak_J=stabilization,
                            external_work_J=work if integrate_work else None,
                            damage_dissipation_J=ledger.damage_dissipation_J,
                            failed_weak_bond_fraction=float(np.mean(D[weak] == 1)),
                            numerical_energy_residual_J=(ledger.damage_dissipation_J-work) if integrate_work else None))
        g, shape = new_g, new_shape
        if np.all(D[weak] == 1):
            break
    complete = bool(np.all(D[weak] == 1))
    area = .0005**2 if complete else 0.
    if not complete:
        raise RuntimeError('Coupon protocol ended before complete weak-plane failure; do not assign a crack area')
    value = ledger.damage_dissipation_J/area
    relative_residual = abs(work-ledger.damage_dissipation_J)/max(ledger.damage_dissipation_J, 1e-30) if integrate_work else None
    if integrate_work and relative_residual > .005:
        raise FloatingPointError(f'Coupon force-work integration residual {relative_residual}; refine quadrature')
    return dict(problem='pd_fracture_energy_coupon', spacing_m=spacing, horizon_m=cloud.horizon,
                horizon_ratio=base['clot']['horizon_ratio'], particle_count=len(cloud.X), initial_bonds=B,
                critical_energy_density_Pa=critical_density_Pa, Gc_measured_J_m2=value,
                created_crack_area_m2=area, crack_area_convention='Single projected interface area, both crack faces not double-counted',
                dissipated_fracture_work_J=ledger.damage_dissipation_J, external_work_J=work if integrate_work else None,
                final_stored_energy_after_unloading_J=sum(potential(cloud, evaluate_supported(cloud, cloud.X, g, material, safety, shape))),
                relative_work_residual=relative_residual, weak_plane_complete=complete,
                constitutive_material=material, fatigue_parameters=dc, increments=increments,
                maximum_coupon_cycles=3*increments, opening_increment_m=max_opening/increments,
                elapsed_s=time.perf_counter()-start, history=history,
                bond_critical_energy_J_min=float(critical.min()), bond_critical_energy_J_max=float(critical.max()),
                calibration_scope='Kinematic cyclic-opening weak-plane coupon; mesh transfer beyond this loading path must be verified')


def calibrate(base, spacing, Gc_demo, output, tolerance=.05, increments=180):
    out = Path(output)
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)
    import shutil
    project=Path(__file__).resolve().parents[2]
    source_paths=[Path(__file__),Path(__file__).with_name('energy.py'),project/'pd_clot/geometry.py',project/'pd_clot/fragment_mechanics.py',project/'pd_clot/mechanics.py']
    source_hashes={}
    for source in source_paths:
        name=source.relative_to(project);target=out/'source_snapshot'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        source_hashes[str(name)]=hashlib.sha256(source.read_bytes()).hexdigest()
    (out/'SOURCE_IDENTITY.json').write_text(json.dumps(source_hashes,indent=2)+'\n')
    if Gc_demo <= 0:
        raise ValueError('Gc_demo must be an explicit positive input in J/m^2')
    evaluations = []

    def objective(log_density):
        psi = float(np.exp(log_density))
        result = run_coupon(base, spacing, psi, increments, integrate_work=False)
        evaluations.append(dict(critical_density_Pa=psi, measured_Gc_J_m2=result['Gc_measured_J_m2']))
        print(json.dumps(dict(calibration_spacing_m=spacing, **evaluations[-1])), flush=True)
        return np.log(result['Gc_measured_J_m2']/Gc_demo)

    center = np.log(Gc_demo/spacing)
    low, high = center-12, center+2
    while objective(low) > 0:
        low -= 4
    while objective(high) < 0:
        high += 2
    solution = brentq(objective, low, high, xtol=min(tolerance/4, .005), maxiter=30)
    final = run_coupon(base, spacing, float(np.exp(solution)), increments, integrate_work=True)
    final.update(Gc_demo_J_m2=Gc_demo, relative_calibration_error=abs(final['Gc_measured_J_m2']/Gc_demo-1),
                 tolerance=tolerance, scalar_calibration_history=evaluations,
                 baseline_config_sha256=hashlib.sha256(json.dumps(base, sort_keys=True).encode()).hexdigest(),
                 experimental_calibration=False, D_break=1.)
    assert final['relative_calibration_error'] <= tolerance
    (out/'CALIBRATION.json').write_text(json.dumps(final, indent=2, allow_nan=False)+'\n')
    return final


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--spacing', type=float, required=True)
    p.add_argument('--Gc-demo', type=float, required=True, help='Demonstration target, J/m^2; never inferred')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--increments', type=int, default=180)
    args = p.parse_args()
    calibrate(json.loads(args.config.read_text()), args.spacing, args.Gc_demo, args.output, increments=args.increments)


if __name__ == '__main__':
    main()
