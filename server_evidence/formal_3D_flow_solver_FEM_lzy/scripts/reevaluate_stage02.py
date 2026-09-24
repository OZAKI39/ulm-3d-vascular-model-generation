#!/usr/bin/env python3
"""Re-evaluate saved fields after validating curved-cell point location."""
import json
import shutil
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from mpi4py import MPI
from fem3d.result_io import load_primary
from fem3d.benchmark_geometry import load_benchmark_config
from fem3d.diagnostics import analytic_diagnostics
from fem3d.audit import write_json

config,_=load_benchmark_config(ROOT)
for base in sorted((ROOT/"outputs/stage02/cases").iterdir()):
    if not (base/"metadata/run.json").is_file():
        continue
    domain,tags,fields,metadata,displacements=load_primary(ROOT,base.name,MPI.COMM_WORLD)
    result=analytic_diagnostics(domain,fields["velocity"],fields["pressure"],fields["lambda_pa"],config["geometry"],
        metadata["mu_pa_s"],metadata["Q_m3_s"],metadata["rho_kg_m3"],ROOT/"outputs/stage02/jit_cache")
    result["exact_solution_for_this_boundary_model"]=metadata["mode"]=="reference"
    result["point_location_method"]="Native nonlinear pull-back verifies reference-tetra membership; eval tol=1e-12"
    if MPI.COMM_WORLD.rank==0:
        path=base/"qc/analytic_comparison.json"
        old=path.with_name("analytic_initial_convex_hull_sampling.json")
        if not old.exists():
            shutil.copy2(path,old)
        write_json(path,result)
        print(base.name,{k:result[k] for k in ("lambda_relative_error","velocity_L2_relative_error","velocity_global_L2_relative_error","pressure_profile_error")},flush=True)
