#!/usr/bin/env python3
"""New-process coefficient reload, derived-value reload and FEM flux recheck."""
import argparse
import json
import os
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
import numpy as np
from mpi4py import MPI
from scipy.spatial import cKDTree
from dolfinx import io,mesh
from fem3d.result_io import load_primary
from fem3d.derived_fields import derive
from fem3d.diagnostics import flux_diagnostics
from fem3d.audit import sha256,timestamp,write_json

p=argparse.ArgumentParser(); p.add_argument("--case"); a=p.parse_args()
gate=json.loads((ROOT/"inputs/stage02/core_gate.json").read_text())
cases=[a.case] if a.case else gate["case_names"]
assert MPI.COMM_WORLD.size==1
for case in cases:
    base=ROOT/"outputs/stage02/cases"/case
    manifest=json.loads((base/"solution/export_manifest.json").read_text())
    assert manifest["export_pid"]!=os.getpid()
    for name,expected in manifest["files"].items():
        assert sha256(base/"solution"/name)==expected
    domain,tags,primary,metadata,displacements=load_primary(ROOT,case,MPI.COMM_WORLD)
    # The XDMF/HDF5 mesh is independently opened as well as the exact checkpoint.
    with io.XDMFFile(domain.comm,str(base/"solution/fields.xdmf"),"r") as file:
        exported_mesh=file.read_mesh(name="pipe_mesh",ghost_mode=mesh.GhostMode.shared_facet)
        exported_mesh.topology.create_connectivity(2,3)
        exported_tags=file.read_meshtags(exported_mesh,name="facet_tags")
    assert exported_mesh.topology.index_map(3).size_global==domain.topology.index_map(3).size_global
    assert set(exported_tags.values)=={1,2,4}
    flux=flux_diagnostics(domain,tags,primary["velocity"],metadata["Q_m3_s"],ROOT/"outputs/stage02/jit_cache")
    original=json.loads((base/"qc/flux.json").read_text())
    flux_diff=max(abs(flux[key]-original[key]) for key in ("actual_Q_in_m3_s","actual_Q_out_m3_s"))/abs(metadata["Q_m3_s"])
    assert flux_diff<1e-11 and flux["status"]=="PASS"
    derived,_=derive(primary["velocity"],ROOT/"outputs/stage02/jit_cache")
    stored=np.load(base/"solution/derived_checkpoint.npz")
    coordinates=next(iter(derived.values())).function_space.tabulate_dof_coordinates()
    distance,index=cKDTree(stored["coordinates_m"]).query(coordinates)
    assert distance.max()<1e-17
    errors={}
    for name,field in derived.items():
        current=field.x.array.reshape(-1,field.function_space.dofmap.index_map_bs)
        saved=stored[name][index]
        errors[name]=float(np.linalg.norm(current-saved)/max(np.linalg.norm(saved),1e-300))
        assert errors[name]<1e-12
    grad=stored["velocity_gradient_s_inv"].reshape(-1,3,3)
    strain=stored["strain_rate_tensor_s_inv"].reshape(-1,3,3)
    curl=np.c_[grad[:,2,1]-grad[:,1,2],grad[:,0,2]-grad[:,2,0],grad[:,1,0]-grad[:,0,1]]
    np.testing.assert_allclose(strain,.5*(grad+grad.transpose(0,2,1)),atol=1e-12,rtol=1e-12)
    np.testing.assert_allclose(stored["vorticity_s_inv"],curl,atol=1e-12,rtol=1e-12)
    result={"status":"PASS","timestamp":timestamp(),"case":case,"new_process":True,"reload_pid":os.getpid(),
        "export_pid":manifest["export_pid"],"relative_flux_difference":flux_diff,
        "primary_coefficients_equal":True,"derived_values_equal":True,"derived_relative_errors":errors,
        "coordinate_match_displacements_m":displacements,"tensor_component_order_verified":True,
        "xdmf_hdf5_mesh_reload":"PASS","lambda_pa":primary["lambda_pa"]}
    write_json(base/"qc/result_reload.json",result)
    print(case,"new-process reload PASS",flush=True)
