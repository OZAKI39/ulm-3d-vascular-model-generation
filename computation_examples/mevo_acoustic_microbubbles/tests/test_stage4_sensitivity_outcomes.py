"""Numerical sensitivity must distinguish lost arrivals from conditional preference."""
import json
from pathlib import Path
import h5py
import numpy as np


def test_sensitivity_keeps_initial_population_denominators_and_censoring(tmp_path,monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]))
    from scripts.stage4 import convergence
    monkeypatch.setattr(convergence,"ROOT",tmp_path)
    for name,assignment,na,nb,nu,crossings in [
        ("reference",{"0":"BRANCH_A","1":"BRANCH_A","2":"BRANCH_B","3":"UNASSIGNED"},2,1,1,
         [dict(particle_id=i,crossing_time_s=t) for i,t in [(0,1.),(1,2.),(2,3.)]]),
        ("candidate",{"0":"BRANCH_A","1":"UNASSIGNED","2":"BRANCH_A","3":"UNASSIGNED"},2,0,2,
         [dict(particle_id=i,crossing_time_s=t) for i,t in [(0,1.2),(2,3.4)]])]:
        out=tmp_path/name
        (out/"particles").mkdir(parents=True)
        (out/"BRANCH_METRICS.json").write_text(json.dumps(dict(assignment=assignment,N_A=na,N_B=nb,N_unassigned=nu,f_A=na/(na+nb),crossings=crossings)))
        (out/"particles/SUMMARY.json").write_text(json.dumps(dict(accepted_steps=12000 if name=="reference" else 24000,
            nominal_dt_s=.0005 if name=="reference" else .00025,minimum_dt_s=1e-8,duration_s=6.)))
        (out/"metrics.json").write_text(json.dumps(dict(rows=[dict(time_s=t,mean_nearest_neighbor_m=1e-5,
            near_contact_fraction=0.,largest_cluster_fraction=.25,cluster_count=4) for t in [0.,.05,.1,6.]])))
        with h5py.File(out/"particles/trajectories.h5","w") as h:
            h["positions_m"]=np.zeros((2,4,3))
    row=convergence.compare(dict(run="reference",initial_population_sha256="same"),dict(run="candidate",initial_population_sha256="same"))
    assert row["assignment_agreement"]==.5
    assert row["changed_particle_ids"]==[1,2]
    assert row["f_A_difference"]>0
    assert row["A_over_initial_difference"]==0
    assert row["B_over_initial_difference"]==-.25
    assert row["branch_switch_matrix"]["BRANCH_A"]["UNASSIGNED"]==1
    assert row["branch_switch_matrix"]["BRANCH_B"]["BRANCH_A"]==1
    assert row["shared_outlet_crossings"]==2
    assert abs(row["mean_absolute_crossing_time_difference_s"]-.3)<1e-14
    assert all(v["relative_L2_difference"]==0 for v in row["early_cluster_metrics"].values())
    assert row["effective_timestep"]["candidate"]["mean_effective_dt_s"]==.00025
