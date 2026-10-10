"""Paired-seed confidence intervals and complete branch-transition tables."""
import csv
import json
import numpy as np
import yaml
from scripts.stage4.common import ROOT, new_run, require, write, finish


def interval(values,indices):
    a = np.asarray(values,float)
    if not np.isfinite(a).all():
        return dict(mean=None,std=None,bootstrap_95CI=None,
                    values=[float(v) if np.isfinite(v) else None for v in a],
                    defined_seed_count=int(np.isfinite(a).sum()),total_seed_count=len(a),
                    warning="A conditional branch fraction is undefined when no particles exit. No seed deletion or imputation; use counts/fractions of the initialized population.")
    return dict(mean=float(a.mean()),std=float(a.std(ddof=1)) if len(a)>1 else None,
                bootstrap_95CI=np.quantile(a[indices].mean(axis=1),[.025,.975]).tolist(),values=a.tolist())


def main():
    source = require("statistics/paired_runs")
    protocol = yaml.safe_load((ROOT / "configs/stage4/experiment.yaml").read_text())["statistics"]
    out = new_run("statistics",protocol,["statistics/paired_runs"])
    rows = source["rows"]
    seeds = sorted({r["seed"] for r in rows})
    assert len(seeds) >= 10
    lookup = {(r["seed"],r["case"]):r for r in rows}
    for seed in seeds:
        assert len({lookup[seed,label]["initial_population_sha256"] for label in ["B0","B1","B2","B3"]}) == 1
    indices = np.random.default_rng(protocol["bootstrap_seed"]).integers(0,len(seeds),(protocol["bootstrap_repetitions"],len(seeds)))
    cases,arrays = {},{}
    for label in ["B0","B1","B2","B3"]:
        values = [lookup[seed,label] for seed in seeds]
        arrays[label] = {}
        for key in ["f_A","f_B","N_A","N_B","N_unassigned","N_wall_loss"]:
            arrays[label][key] = np.array([r[key] for r in values],float)
        for key,count_key in [("A_fraction_of_initial","N_A"),("B_fraction_of_initial","N_B"),("unassigned_fraction","N_unassigned")]:
            arrays[label][key] = arrays[label][count_key]/np.array([sum(r[k] for k in ["N_A","N_B","N_unassigned","N_wall_loss"]) for r in values])
        cases[label] = {key:interval(a,indices) for key,a in arrays[label].items()}
    effects = {}
    for label,reference,key in [("B1","B0","f_A"),("B2","B0","f_A"),("B3","B0","f_B"),
                                 ("B2","B3","f_A"),("B2","B0","A_fraction_of_initial"),
                                 ("B3","B0","B_fraction_of_initial"),("B2","B3","A_fraction_of_initial")]:
        effects[label+"_minus_"+reference+"_"+key] = interval(arrays[label][key]-arrays[reference][key],indices)
    effects["mirror_fA_sum_minus_1"] = interval(arrays["B2"]["f_A"]+arrays["B3"]["f_A"]-1,indices)
    transition_rows = []
    outcomes = ["BRANCH_A","BRANCH_B","UNASSIGNED","WALL_LOSS"]
    pooled_transitions = {"B0->"+label:{origin:dict.fromkeys(outcomes,0) for origin in outcomes}
                          for label in ["B1","B2","B3"]}
    distribution_rows = []
    for seed in seeds:
        loaded = {label:json.loads((ROOT / lookup[seed,label]["run"] / "BRANCH_METRICS.json").read_text()) for label in cases}
        for label in ["B1","B2","B3"]:
            assert loaded["B0"]["assignment"].keys()==loaded[label]["assignment"].keys()
            matrix = {origin+"->"+destination:0 for origin in outcomes for destination in outcomes}
            for pid,baseline in loaded["B0"]["assignment"].items():
                control = loaded[label]["assignment"][pid]
                key = baseline+"->"+control
                matrix[key] += 1
                pooled_transitions["B0->"+label][baseline][control] += 1
            for transition,count in sorted(matrix.items()):
                transition_rows.append(dict(seed=seed,comparison="B0->"+label,transition=transition,count=count))
        for label,branch in loaded.items():
            junction = {e["particle_id"]:e for e in branch["bifurcation_crossings"]}
            crossings = {e["particle_id"]:e for e in branch["crossings"]}
            summary = json.loads((ROOT / lookup[seed,label]["run"] / "particles/SUMMARY.json").read_text())
            for pid,assignment in branch["assignment"].items():
                pid = int(pid)
                event = crossings.get(pid,{})
                bif = junction.get(pid,{})
                distribution_rows.append(dict(seed=seed,case=label,particle_id=pid,branch=assignment,
                    crossing_time_s=event.get("crossing_time_s"),residence_time_s=event.get("crossing_time_s",summary["duration_s"]),
                    right_censored=assignment == "UNASSIGNED",bifurcation_y_m=bif.get("position_m",[None]*3)[1],
                    bifurcation_z_m=bif.get("position_m",[None]*3)[2],
                    bifurcation_transverse_speed_m_s=float(np.linalg.norm(bif["transverse_velocity_m_s"])) if bif else None))
    for name,data in [("BRANCH_SWITCH_TABLE",transition_rows),("PARTICLE_OUTCOMES",distribution_rows)]:
        write(out / (name+".json"),data)
        with (out / (name+".csv")).open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    supports = all(effects[key]["bootstrap_95CI"] is not None and effects[key]["bootstrap_95CI"][0] > 0 for key in
                   ["B2_minus_B0_f_A","B3_minus_B0_f_B","B2_minus_B3_f_A"])
    finish(out,dict(status="PASS_WITH_SCOPE",cases=cases,paired_effects=effects,seed_count=len(seeds),seeds=seeds,
                    bootstrap=protocol,paired_identity_verified=True,statistical_directional_control_supported=supports,
                    pooled_transition_counts=pooled_transitions,
                    transition_scope="All A/B/U/W combinations including zero counts. Pooled matched-particle counts are descriptive, not independent particle replicates for inference.",
                    residence_time_origin="Elapsed from the common initialization plane at t=0, not from the upstream CFD inlet; right-censored rows stop at the observation cutoff",
                    control_effect_definition="Requested f_target=N_target/N_assigned, conditional on observed exit. Fractions of the initial ensemble and unassigned fractions are separate outcomes; no inference of increased absolute target delivery from conditional enrichment alone.",
                    scope="Resample independent seeds as paired blocks; particles within a run are not treated as independent seed replicates. Conditional-on-assignment and fractions of all initialized bubbles both reported; numerical/model robustness remains separate."))


if __name__ == "__main__":
    main()
