#!/usr/bin/env python3
"""Generate P2 synthetic or real validation data; every dt is VALIDATION_ONLY."""
import argparse
import csv
import json
from pathlib import Path
import sys
import time
import numpy as np

PACKAGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKAGE/"src"))
from particle_3d.rbc_distribution import sample_rbc_geometries,load_contract,quantile_indices,stratified_indices,digest
from particle_3d.rbc import RBCGeometry
from particle_3d.particle2_cases import static_case,rigid_rotation_case,simple_shear_case,real_rbc_trajectory,SHEAR_ANGLE_STEPS,NORM_ATOL,IDENTITY_RELATIVE_BUDGET,PERIOD_RELATIVE_BOUND_FACTOR
from particle_3d.audit import read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.particle1_cases import select_inlet_centroid
from particle_3d.validation_boundary import ValidationBoundaryClassifier

DATA=PACKAGE/"reports/particle2/data"


def write_json(name,value):
    (DATA/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False,
                                    default=lambda x:x.item() if isinstance(x,np.generic) else x)+"\n")


def write_csv(name,rows):
    with (DATA/name).open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--stage",choices=["synthetic","real"],required=True);args=parser.parse_args()
    c=load_contract();pop=sample_rbc_geometries(c["validation"]["N"],c["validation"]["seed"])
    meta=json.loads((DATA/"C57BL6_RBC_GEOMETRY_VALIDATION_100000.metadata.json").read_text())
    if pop.metadata["sample_structured_array_sha256"]!=meta["sample_structured_array_sha256"]:
        raise ValueError("saved geometry population differs from reproducible sampler")
    selected=quantile_indices(pop.samples);geometries=[RBCGeometry.from_population(pop,i) for i in selected]
    if args.stage=="synthetic":
        # Persist budgets before results, never fitted to observed errors.
        write_json("synthetic_validation_budgets.json",dict(norm_atol=NORM_ATOL,identity_relative_budget=IDENTITY_RELATIVE_BUDGET,
            shear_gamma_dt=list(SHEAR_ANGLE_STEPS),period_relative_bound_rule=f"{PERIOD_RELATIVE_BOUND_FACTOR}*gamma*dt",
            production_particle_timestep_frozen=False))
        for prefix,fn in [("04_static",static_case),("05_rotation",rigid_rotation_case)]:
            rows,metrics=fn(geometries);write_csv(prefix+".csv",rows);write_json(prefix+"_metrics.json",metrics)
            print(prefix,metrics,flush=True)
        all_metrics=[]
        for i,angle in enumerate(SHEAR_ANGLE_STEPS):
            rows,metrics=simple_shear_case(geometries,angle,progress=lambda n,total:print(f"shear dt{i}: {n}/{total}",flush=True))
            write_csv(f"06_shear_dt{i}.csv",rows);write_json(f"06_shear_dt{i}_metrics.json",metrics);all_metrics.extend(metrics)
            print("shear",i,"max period relative error",max(m["relative_error"] for m in metrics),flush=True)
        wide=[RBCGeometry.from_population(pop,i) for i in stratified_indices(pop.samples,c["validation"]["orientation_selection_seed"])]
        rows,metrics=simple_shear_case(wide,SHEAR_ANGLE_STEPS[-1],progress=lambda n,total:print(f"64 geometries: {n}/{total}",flush=True))
        write_csv("07_distribution_wide_shear.csv",rows);write_csv("07_distribution_wide_metrics.csv",metrics);write_json("07_distribution_wide_metrics.json",metrics)
        if not all(m["passed"] for m in all_metrics+metrics):raise SystemExit("synthetic validation failed; preserve evidence")
    else:
        audit,mesh,flow,boundaries=read_frozen(PACKAGE.parent/"formal_3D_flow_solver/FEM_SimVascular")
        field=FrozenFEMField.from_grids(mesh,flow);classifier=ValidationBoundaryClassifier(boundaries)
        initial,candidates=select_inlet_centroid(field,boundaries["INLET"])
        # Reuse the three P1 geometry/velocity-defined validation values unchanged.
        initial["dt_provenance"]="same three Particle-1 validated dt values; no new production selection"
        write_json("08_real_initialization.json",initial);write_csv("08_inlet_candidates.csv",candidates)
        write_json("08_frozen_provenance.json",audit)
        metrics=[]
        for gi,geometry in enumerate(geometries):
            for di,dt in enumerate(initial["validation_timesteps_s"]):
                started=time.perf_counter()
                rows,summary=real_rbc_trajectory(field,classifier,initial,geometry,dt,
                    progress=lambda step:print(f"real geometry={gi} dt={di}: step {step}",flush=True))
                summary.update(geometry_index=gi,dt_index=di,elapsed_wall_time_s=time.perf_counter()-started)
                write_csv(f"08_real_g{gi}_dt{di}.csv",rows);write_json(f"08_real_g{gi}_dt{di}.json",summary);metrics.append(summary)
                print("real",gi,di,summary["exit_boundary"],"PASS" if summary["passed"] else "FAIL",flush=True)
        write_json("08_real_sweep.json",metrics)
        if not all(m["passed"] for m in metrics):raise SystemExit("real validation failed; preserve evidence")


if __name__=="__main__":main()
