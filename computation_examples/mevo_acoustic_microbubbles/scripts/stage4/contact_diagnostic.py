"""Replay only the rejected step to diagnose iteration budget, not a resumed batch."""
import copy
import json
import numpy as np
from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.coupling.bifurcation_experiment import BifurcationExperimentTracker
from scripts.stage4.common import ROOT,new_run,write,finish
from scripts.stage4.experiments import setup


def main():
    source = ROOT / "results/stage4/primary_only/20261002T145553771005Z/20261002T145727_818459Z_B1_seed20261002"
    original = json.loads((source / "RUN_IDENTITY.json").read_text())["config"]
    state = np.load(source / "particles/last_accepted_state.npz")
    x,r,t = state["positions_m"],state["radii_m"],float(state["time_s"])
    flow = setup()[3]
    field = ComplexAcousticField.load(ROOT / original["pzt_control"]["field"])
    out = new_run("performance/contact_diagnostic",dict(source=str(source.relative_to(ROOT)),time_s=t,
          budgets=[100,300,1000,3000],contact_tolerance_m=original["particles"]["contact_tolerance_m"]),["phase_screening"])
    rows,values = [],[]
    for budget in [100,300,1000,3000]:
        config = copy.deepcopy(original)
        config["particles"]["contact_max_iterations"] = budget
        tracker = BifurcationExperimentTracker(config,field,flow)
        first = tracker.evaluate(x,r,t)
        dt = min(config["particles"]["dt_s"],float(np.min(config["particles"]["max_displacement_per_radius"]*r/np.maximum(np.linalg.norm(first["total_velocity_m_s"],axis=1),1e-30))))
        try:
            value,wall,contact,correction,iterations = tracker.step(x,r,t,dt,first)
            values.append(value)
            rows.append(dict(budget=budget,status="CONVERGED",dt_s=dt,maximum_correction_m=correction,
                             final_projection_iterations=iterations,wall_particles=int(wall.sum()),contact_pairs=len(contact)))
            np.save(out / ("step_budget_"+str(budget)+".npy"),value)
        except RuntimeError as error:
            rows.append(dict(budget=budget,status="REJECTED",dt_s=dt,error=repr(error)))
        write(out / "ROWS.json",rows)
    same = len(values)>=2 and all(np.array_equal(value,values[0]) for value in values[1:])
    finish(out,dict(status="PASS" if same else "FAIL",rows=rows,converged_step_values_bitwise_equal=same,
                   scope="Iteration-cap test only; tolerance, force law, time step and contact geometry unchanged"))


if __name__ == "__main__":
    main()
