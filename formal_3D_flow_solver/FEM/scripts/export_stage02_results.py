#!/usr/bin/env python3
"""After the core gate, restore coefficients and export reviewable fields."""
import argparse
import json
import os
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
from mpi4py import MPI
from dolfinx import fem,io,plot
import basix.ufl
from fem3d.result_io import load_primary
from fem3d.derived_fields import derive
from fem3d.audit import sha256,timestamp,write_json

p=argparse.ArgumentParser()
p.add_argument("--case")
a=p.parse_args()
gate=json.loads((ROOT/"inputs/stage02/core_gate.json").read_text())
assert gate["status"]=="PASS" and MPI.COMM_WORLD.size==1
cases=[a.case] if a.case else gate["case_names"]
for case in cases:
    base=ROOT/"outputs/stage02/cases"/case
    if (base/"solution/export_manifest.json").exists():
        raise RuntimeError(f"Already exported {case}; preserve prior results")
    domain,tags,primary,metadata,displacements=load_primary(ROOT,case,MPI.COMM_WORLD)
    fields,definitions=derive(primary["velocity"],ROOT/"outputs/stage02/jit_cache")
    P2=fem.functionspace(domain,basix.ufl.element("Lagrange",domain.basix_cell(),2))
    pressure_vis=fem.Function(P2,name="pressure_gauge_pa")
    pressure_vis.interpolate(primary["pressure"])
    pressure_vis.x.scatter_forward()
    xdmf=base/"solution/fields.xdmf"
    with io.XDMFFile(domain.comm,str(xdmf),"w") as file:
        file.write_mesh(domain)
        file.write_meshtags(tags,domain.geometry)
        file.write_function(primary["velocity"])
        file.write_function(pressure_vis)
        for field in fields.values():
            file.write_function(field)
    derived={name:field.x.array.reshape(-1,field.function_space.dofmap.index_map_bs).copy() for name,field in fields.items()}
    coordinates=next(iter(fields.values())).function_space.tabulate_dof_coordinates()
    np.savez_compressed(base/"solution/derived_checkpoint.npz",coordinates_m=coordinates,**derived)
    topology,types,points=plot.vtk_mesh(primary["velocity"].function_space)
    # Tabulate scalar pressure at exactly the velocity/VTK nodes for rendering.
    scalar_points=P2.tabulate_dof_coordinates()
    from scipy.spatial import cKDTree
    distances,indices=cKDTree(scalar_points).query(points)
    assert distances.max()==0
    np.savez_compressed(base/"solution/visualization_mesh.npz",topology=topology,cell_types=types,points_m=points,
        velocity_m_s=primary["velocity"].x.array.reshape(-1,3),pressure_gauge_pa=pressure_vis.x.array[indices])
    manifest={"timestamp":timestamp(),"export_pid":os.getpid(),"run_id":os.environ.get("FEM3D_RUN_ID"),
        "case":case,"status":"PASS","lambda_pa":primary["lambda_pa"],"lambda_units":"Pa",
        "primary_restart":"primary_checkpoint.npz + the profile's pipe.xdmf/pipe.h5; exact P2/P1 nodal coefficients",
        "xdmf_fields":"velocity P2, pressure P1 embedded into P2 for geometry-compatible output, centroid DG0 derived fields",
        "definitions":definitions,"coordinate_match_displacement_m":displacements,
        "files":{path.name:sha256(path) for path in (base/"solution").iterdir() if path.is_file()}}
    write_json(base/"solution/export_manifest.json",manifest)
    print(case,"exported primary and derived fields",flush=True)
