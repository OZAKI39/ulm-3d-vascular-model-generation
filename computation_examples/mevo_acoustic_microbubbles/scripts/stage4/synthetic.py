"""Stage 4C: incompressible synthetic test field, same force-free RK path and IDs."""
import copy
import numpy as np
import yaml
from mevo_mb.bubbles.inlet_sampling import initialize_inlet
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from mevo_mb.coupling.bifurcation_experiment import BifurcationExperimentTracker
from mevo_mb.flow.bifurcation import SyntheticBifurcationFlowField
from mevo_mb.io.artifacts import new_run as particle_run
from mevo_mb.utils.config import load_config
from scripts.stage4.common import ROOT, new_run, write, finish


class AcousticOff:
    frequency_Hz = 2.5e6

    def interpolate_complex_pressure(self, x):
        return np.zeros(len(x), complex)

    def interpolate_complex_gradient(self, x):
        return np.zeros_like(x, dtype=complex)


def main():
    c = yaml.safe_load((ROOT / "configs/stage4/bifurcation.yaml").read_text())
    out = new_run("synthetic_flow", c, ["geometry"])
    g = BifurcationGeometry(c["geometry"])
    flow = SyntheticBifurcationFlowField(g, c["synthetic_flow"])
    rng = np.random.default_rng(901)
    points = rng.uniform([-.003,-.003,-.001], [.01,.003,.001], (1000,3))
    div = np.zeros(len(points))
    h = 1e-7
    for axis in range(3):
        delta = np.eye(3)[axis]*h
        div += (flow.velocity(points+delta,0)[:,axis]-flow.velocity(points-delta,0)[:,axis])/(2*h)
    velocity = flow.velocity(points,0)
    base = load_config(ROOT / "configs/stage3/cylinder.yaml")
    base["geometry"] = dict(c["geometry"], kind="bifurcation")
    base["experiment_scope"] = "SOFTWARE_TEST_ONLY"
    base["particles"].update(primary_enabled=False, secondary_enabled=False,
                              dt_s=c["synthetic_flow"]["dt_s"], duration_s=c["synthetic_flow"]["duration_s"],
                              max_displacement_per_radius=20., pair_backend="numpy")
    base["output"]["save_interval_s"] = .02
    rows = []
    for policy in ["uniform_area", "flux_weighted"]:
        sampling = dict(c["particles"], sampling_policy=policy)
        population = initialize_inlet(g, flow, sampling)
        repeat = initialize_inlet(g, flow, sampling)
        assert np.array_equal(population.positions_m, repeat.positions_m)
        cc = copy.deepcopy(base)
        cc["sampling"] = sampling
        run = particle_run(out, policy, cc)
        np.savez_compressed(run / "INITIAL_POPULATION.npz", positions_m=population.positions_m,
                            radii_m=population.radii_m, particle_id=population.ids)
        tracker = BifurcationExperimentTracker(cc, AcousticOff(), flow)
        summary = tracker.run(population, run)
        branch = summary["branch_metrics"]
        routed = [(int(pid), label) for pid,label in branch["assignment"].items() if label.startswith("BRANCH_")]
        correct = all((population.positions_m[pid,1] > 0) == (label == "BRANCH_A") for pid,label in routed)
        rows.append(dict(policy=policy, run=str(run.relative_to(ROOT)),
                         N_A=branch["N_A"], N_B=branch["N_B"], N_unassigned=branch["N_unassigned"],
                         exact_sampling_repeat=True, sign_routing_correct=correct,
                         crossings_recorded=len(branch["crossings"]), elapsed_s=summary["elapsed_s"]))
        write(out / "ROWS.json", rows)
    tests = dict(finite=np.isfinite(velocity).all().item(),
                 bounded=float(np.linalg.norm(velocity,axis=1).max()) < flow.maximum_velocity*1.2,
                 incompressible_test_field=float(np.abs(div).max()) < 1e-5,
                 bookkeeping=all(r["N_A"]+r["N_B"]+r["N_unassigned"] == c["particles"]["count"] for r in rows),
                 routes=all(r["sign_routing_correct"] and r["N_A"]>0 and r["N_B"]>0 for r in rows),
                 assigned_majority=all(r["N_unassigned"] < .1*c["particles"]["count"] for r in rows))
    finish(out, dict(status="PASS_WITH_SCOPE" if all(tests.values()) else "FAIL", tests=tests, rows=rows,
                     maximum_numerical_divergence_s_inverse=float(np.abs(div).max()),
                     scope="SOFTWARE_TEST_ONLY; analytic streamfunction field, not a CFD or physical steering result"))


if __name__ == "__main__":
    main()
