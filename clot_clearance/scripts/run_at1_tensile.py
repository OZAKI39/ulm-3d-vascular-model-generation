"""Run the two predeclared paths plus an extra-iteration comparison."""
from concurrent.futures import ProcessPoolExecutor
import copy
import hashlib
import json
import numpy as np
from scripts.verify_at1_evolution import ROOT, OUT, config, worker_budget, write_json
from pd_clot.pfpd_candidate.staggered_at1 import run_path


def job(task):
    name, c, path = task
    return name, run_path(c, c['schedules'][path],
                         ROOT/'results/pfpd_at1_evolution_v1'/name)


def main():
    gate = json.loads((OUT/'analytic/UNIT_GATES.json').read_text())
    if gate['status'] != 'PASSED': raise RuntimeError('UNIT_GATES must pass first')
    unit_log = (OUT/'tests/UNIT.log').read_text()
    if '34 passed' not in unit_log or 'FAILED' in unit_log:
        raise RuntimeError('Scalar/fixed-x regression tests must pass first')
    c = config(); extra = copy.deepcopy(c)
    extra['staggered']['extra_converged_iterations'] = 2
    budget = worker_budget(); write_json(OUT/'provenance/COUPLED_WORKER_BUDGET.json',budget)
    tasks = [('tensile_cube',c,'tensile_cube'), ('loading_unloading',c,'loading_unloading'),
             ('loading_unloading_extra_iterations',extra,'loading_unloading')]
    results = {}
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        for name, summary in pool.map(job,tasks):
            results[name] = summary
            print(name, json.dumps(summary), flush=True)
    write_json(OUT/'staggered/ALL_RUNS.json',results)
    comparison = []
    base = ROOT/'results/pfpd_at1_evolution_v1'
    for step in range(len(c['schedules']['loading_unloading'])):
        paths = [base/name/f'load_{step:03d}/accepted_state.npz' for name in
                 ['loading_unloading','loading_unloading_extra_iterations']]
        if not all(p.exists() for p in paths): break
        with np.load(paths[0]) as a, np.load(paths[1]) as b:
            comparison.append(dict(step=step, max_delta_x_m=float(np.max(abs(a['x']-b['x']))),
                                   max_delta_s=float(np.max(abs(a['s_pf_bond']-b['s_pf_bond'])))))
    write_json(OUT/'staggered/EXTRA_ITERATIONS.json',comparison)
    passed = (results['tensile_cube']['status']=='HIGH_DAMAGE_GATE_REACHED' and
              results['loading_unloading']['status']=='CONVERGED' and
              results['loading_unloading_extra_iterations']['status']=='CONVERGED' and
              len(comparison)==len(c['schedules']['loading_unloading']) and
              all(r['max_delta_x_m']<=c['tolerances']['extra_iteration_x_absolute_m'] and
                  r['max_delta_s']<=c['tolerances']['extra_iteration_state_absolute_s'] for r in comparison))
    write_json(OUT/'staggered/COUPLED_GATE.json',dict(passed=passed,
        config_sha256=hashlib.sha256((OUT/'provenance/PREDECLARED_CONFIG.json').read_bytes()).hexdigest()))
    if not passed: raise SystemExit(1)


if __name__ == '__main__':
    main()
