#!/usr/bin/env python3
"""One small canonical case; no vascular input path is accepted."""
import argparse
import json
import os
import platform
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
from dolfinx import io,mesh
from mpi4py import MPI
from fem3d.audit import sha256,timestamp,write_json
from fem3d.benchmark_geometry import load_benchmark_config
from fem3d.solver import solve_stokes,_solve
from fem3d.diagnostics import flux_diagnostics,analytic_diagnostics


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--profile",required=True,choices=("pipe_coarse","pipe_medium","pipe_fine"))
    p.add_argument("--case",required=True)
    p.add_argument("--mode",choices=("natural","reference"),default="natural")
    p.add_argument("--q-factor",type=float,default=1.)
    p.add_argument("--mu-factor",type=float,default=1.)
    p.add_argument("--rho-factor",type=float,default=1.)
    args=p.parse_args()
    if not args.case.replace("_","").isalnum():
        raise ValueError("Invalid case name")
    config,config_hash=load_benchmark_config(ROOT)
    comm=MPI.COMM_WORLD
    base=ROOT/"outputs/stage02/cases"/args.case
    if (base/"metadata/run.json").exists():
        raise RuntimeError("Existing case evidence must be retained; use a new case name")
    if comm.rank==0:
        for folder in ("metadata","qc","solution"):
            (base/folder).mkdir(exist_ok=True,parents=True)
    comm.barrier()
    start=time.perf_counter()
    meshbase=ROOT/"outputs/stage02/meshes"/args.profile
    with io.XDMFFile(comm,str(meshbase/"pipe.xdmf"),"r") as file:
        domain=file.read_mesh(name="pipe_mesh",ghost_mode=mesh.GhostMode.shared_facet)
        domain.topology.create_connectivity(2,3)
        facet_tags=file.read_meshtags(domain,name="facet_tags")
        cell_tags=file.read_meshtags(domain,name="cell_tags")
    assert np.all(cell_tags.values==100)
    physics=config["physics"]
    mu=physics["dynamic_viscosity_pa_s"]*args.mu_factor
    Q=physics["target_inflow_m3_s"]*args.q_factor
    rho=physics["density_kg_m3"]*args.rho_factor
    cache=ROOT/"outputs/stage02/jit_cache"
    radius=config["geometry"]["radius_m"]
    if args.mode=="natural" and Q>0:
        solved=solve_stokes(domain,facet_tags,mu,Q,radius,cache)
    else:
        solved=_solve(domain,facet_tags,mu,Q,radius,cache,reference_geometry=config["geometry"] if args.mode=="reference" else None)
    flux=flux_diagnostics(domain,facet_tags,solved["velocity"],Q,cache)
    if comm.rank==0:
        write_json(base/"qc/solver.json",solved["solver"])
        write_json(base/"qc/flux.json",flux)
    assert flux["status"]=="PASS",flux
    # Core finite/sign/flux checks passed before analytic comparison begins.
    analytic=analytic_diagnostics(domain,solved["velocity"],solved["pressure"],solved["lambda_pa"],
        config["geometry"],mu,Q,rho,cache)
    analytic["exact_solution_for_this_boundary_model"]=args.mode=="reference"
    checkpoint={}
    for name in ("velocity","pressure"):
        field=solved[name]
        space=field.function_space
        count=space.dofmap.index_map.size_local
        blocksize=space.dofmap.index_map_bs
        local=(space.tabulate_dof_coordinates()[:count].copy(),field.x.array[:count*blocksize].reshape(count,blocksize).copy())
        data=comm.gather(local,root=0)
        if comm.rank==0:
            checkpoint[name+"_coordinates_m"]=np.concatenate([a[0] for a in data])
            checkpoint[name+"_values"]=np.concatenate([a[1] for a in data])
    if comm.rank==0:
        checkpoint["lambda_pa"]=np.array(solved["lambda_pa"])
        np.savez_compressed(base/"solution/primary_checkpoint.npz",**checkpoint)
        write_json(base/"qc/analytic_comparison.json",analytic)
        write_json(base/"metadata/config_snapshot.json",config)
        write_json(base/"metadata/run.json",{"status":"PASS","timestamp":timestamp(),"run_id":os.environ.get("FEM3D_RUN_ID"),
            "hostname":platform.node(),"case":args.case,"profile":args.profile,"mode":args.mode,"purpose":config["purpose"],
            "Q_m3_s":Q,"mu_pa_s":mu,"rho_kg_m3":rho,"mpi_ranks":comm.size,"gpu_used":False,
            "config_sha256":config_hash,"mesh_sha256":sha256(meshbase/"pipe.h5"),
            "source_manifest_sha256":sha256(ROOT/"remote/source_manifest.json"),
            "checkpoint_sha256":sha256(base/"solution/primary_checkpoint.npz"),
            "stdout":os.environ.get("FEM3D_STDOUT_LOG"),"stderr":os.environ.get("FEM3D_STDERR_LOG"),
            "wall_time_s":time.perf_counter()-start,"exit_status":0,"derived_fields_exported":False})
        print(json.dumps({"case":args.case,"solver":solved["solver"],"flux":flux,
            "analytic":{k:v for k,v in analytic.items() if k not in ("velocity_profiles","pressure_sections","velocity_sections")}},indent=2))


if __name__=="__main__":
    main()
