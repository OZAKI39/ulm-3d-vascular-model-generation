"""Independent scalar/energy gates. Run before coupled AT1 benchmarks."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
import importlib.util
import json
import os
import numpy as np
import psutil
from scipy.optimize import minimize_scalar
from pd_clot.pfpd_candidate.at1_evolution import at1_constrained_update, critical_driver, audit_update
from pd_clot.pfpd_candidate.adr_adapter_at1 import fixture, pair_weights, tensile_domain
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.normalization import c0, crack_energy

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'verification/pfpd_at1_evolution_v1'


def write_json(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def config():
    return json.loads((OUT/'provenance/PREDECLARED_CONFIG.json').read_text())


def worker_budget():
    memory = psutil.virtual_memory()
    budget = memory.available - .2*memory.total - 512*2**20
    if budget < 2**30:
        raise RuntimeError('Insufficient RAM for one worker while reserving 20% RAM')
    workers = max(1, min(3, psutil.cpu_count(logical=False) or 1, int(budget/2**30)))
    return dict(max_workers=workers, total_RAM_bytes=memory.total,
                available_RAM_bytes=memory.available, reserve_RAM_bytes=.2*memory.total,
                parent_budget_bytes=512*2**20, per_worker_budget_bytes=2**30,
                threads={k:os.environ.get(k) for k in
                         ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']})


@lru_cache(None)
def external_functions():
    base = OUT/'resources/external/phasefieldx/src/phasefieldx/Element'
    modules = []
    for index, file in enumerate(['Phase_Field/geometric_crack.py',
                                  'Phase_Field_Fracture/g_degradation_functions.py']):
        spec = importlib.util.spec_from_file_location(f'phasefieldx_pure_{index}', base/file)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        modules.append(module)
    return modules


def scalar_oracle(task):
    """Never calls the analytic solution inside either numerical objective."""
    ratio, old, yc, delta, gc = task
    y = ratio*yc
    def energy(s):
        return ((1-s)*(1-s)*y+2*yc*s)/yc
    result = minimize_scalar(energy, bounds=(old, 1.), method='bounded',
                             options=dict(xatol=1e-14, maxiter=1000))
    candidates = [(old, energy(old)), (float(result.x), energy(result.x)), (1., energy(1.))]
    selected = min(candidates, key=lambda a: a[1])
    geometric, degradation = external_functions()
    cw = geometric.geometric_crack_coefficient('AT1')
    gc_fem = gc*cw/c0('quadratic')['c0']
    def external_energy(s):
        return (degradation.quadratic_degradation_function(s)*y+
                gc_fem/(cw*delta)*geometric.geometric_crack_function(s, 'AT1'))/yc
    external = minimize_scalar(external_energy, bounds=(old, 1.), method='bounded',
                               options=dict(xatol=1e-14, maxiter=1000))
    ext_candidates = [(old, external_energy(old)), (float(external.x), external_energy(external.x)),
                      (1., external_energy(1.))]
    ext_selected = min(ext_candidates, key=lambda a: a[1])
    update = at1_constrained_update(old, y, yc)
    analytic = float(update.s_new)
    return dict(ratio=ratio, s_old=old, Y_Pa=y, Yc_Pa=yc,
                analytic_s=analytic, scipy_raw_s=float(result.x), scipy_success=bool(result.success),
                endpoint_candidates=candidates, selected_s=selected[0], selected_energy_over_Yc=selected[1],
                absolute_s_error=abs(analytic-selected[0]),
                external_raw_s=float(external.x), external_success=bool(external.success),
                external_endpoint_candidates=ext_candidates, external_selected_s=ext_selected[0],
                external_absolute_s_error=abs(analytic-ext_selected[0]),
                external_local_coefficient=gc_fem/(cw*delta),
                q=float(update.derivative_at_solution), active_set=str(update.active_set_status),
                audit=audit_update(update, y, yc))


def fixed_x_case(c):
    cloud, f, material, _ = fixture(c)
    x = f.X.copy(); x[:, 2] += .012*(f.X[:, 2]-f.X[:, 2].min())
    old = evaluate(f, x, material)
    minimum = tensile_domain(old, c['tolerances']['tensile_min_green_strain'])
    yc = critical_driver(c['phase']['Gc_demo_J_per_m2'], f.horizon)
    s0 = np.zeros(len(f.pairs)); update = at1_constrained_update(s0, old.W0, yc, c['phase']['sc'])
    new = evaluate(f, x, material, g=(1-update.s_new)**2, h=np.ones(len(f.pairs)))
    weights = pair_weights(f)
    energies = {name:float(np.dot(weights, getattr(update, name))) for name in
                ['elastic_before', 'elastic_after', 'crack_before', 'crack_after', 'total_before', 'total_after']}
    difference = {name:float(np.max(abs(getattr(new, name)-getattr(old, name))))
                  for name in ['K', 'Fbar', 'Ftilde', 'shape_i', 'shape_j']}
    bitwise = {name:bool(np.array_equal(getattr(new, name), getattr(old, name))) for name in difference}
    regression_crack = crack_energy(f, update.s_new, c['phase']['Gc_demo_J_per_m2'], c0()['c0'])
    summary = dict(energies_J=energies, maximum_changes=difference, bitwise_identical=bitwise,
                   minimum_green_strain=minimum, phase_maximum=float(update.s_new.max()),
                   audit=audit_update(update, old.W0, yc), c0=c0(), Yc_Pa=yc,
                   existing_crack_energy_J=regression_crack,
                   existing_crack_energy_difference_J=regression_crack-energies['crack_after'],
                   elastic_identity_error_J=new.energy-energies['elastic_after'])
    np.savez_compressed(OUT/'fixed_x/state.npz', x=x, X=f.X, pairs=f.pairs, weights=weights,
                        Y=old.W0, K_before=old.K, K_after=new.K, Fbar_before=old.Fbar,
                        Fbar_after=new.Fbar, Ftilde_before=old.Ftilde, Ftilde_after=new.Ftilde,
                        **{k:getattr(update,k) for k in update.__dataclass_fields__})
    write_json(OUT/'fixed_x/SUMMARY.json', summary)
    return summary


def main():
    c = config(); budget = worker_budget()
    write_json(OUT/'provenance/SCALAR_WORKER_BUDGET.json', budget)
    delta = c['geometry']['spacing_m']*c['geometry']['horizon_ratio']
    gc = c['phase']['Gc_demo_J_per_m2']; yc = critical_driver(gc, delta)
    tasks = [(r, old, yc, delta, gc) for old in c['s_old_values'] for r in c['scalar_ratios']]
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        rows = list(pool.map(scalar_oracle, tasks))
    write_json(OUT/'minimization/CASES.json', rows)
    history = []; s = 0.
    for step, ratio in enumerate(c['single_bond_path']):
        update = at1_constrained_update(s, ratio*yc, yc)
        row = {k:(str(getattr(update,k)) if k == 'active_set_status' else float(getattr(update,k)))
               for k in update.__dataclass_fields__}
        row.update(step=step, ratio=ratio, audit=audit_update(update, ratio*yc, yc))
        history.append(row); s = float(update.s_new)
    write_json(OUT/'analytic/SINGLE_BOND_PATH.json', history)
    fixed = fixed_x_case(c)
    tolerance = c['tolerances']['analytic_numeric_absolute_s']
    max_error = max(r['absolute_s_error'] for r in rows)
    max_external = max(r['external_absolute_s_error'] for r in rows)
    virgin_threshold = all(r['analytic_s'] == 0 for r in rows if r['s_old'] == 0 and r['ratio'] <= 1)
    energy = fixed['energies_J']
    passed = (max(max_error, max_external) <= tolerance and virgin_threshold and
              all(r['scipy_success'] and r['external_success'] and r['audit']['passed'] for r in rows) and
              all(r['audit']['passed'] for r in history) and all(fixed['bitwise_identical'].values()) and
              energy['total_after'] <= energy['total_before'] and
              energy['elastic_after'] <= energy['elastic_before'] and fixed['audit']['passed'])
    gate = dict(status='PASSED' if passed else 'FAILED', cases=len(rows), max_s_error=max_error,
                max_external_s_error=max_external, tolerance=tolerance, virgin_threshold=virgin_threshold,
                c0=c0(), Yc_Pa=yc, single_bond_irreversibility=all(r['s_new'] >= r['s_old'] for r in history),
                reload_beyond_peak_increases=history[-1]['s_new'] > history[-2]['s_new'])
    write_json(OUT/'analytic/UNIT_GATES.json', gate)
    print(json.dumps(gate, indent=2), flush=True)
    if not passed: raise SystemExit(1)


if __name__ == '__main__':
    main()
