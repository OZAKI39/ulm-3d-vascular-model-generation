#!/usr/bin/env python3
"""Descriptive mesh comparison only; no convergence claim."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fem3d.audit import timestamp, write_json

rows={}
for profile in ("coarse","medium"):
    base=ROOT/"outputs/stage01"/profile
    qc=json.loads((base/"qc/geometry_qc.json").read_text())
    metadata=json.loads((base/"metadata/meshing.json").read_text())
    rows[profile]={"mesh_profile":metadata["mesh_profile"],"requested_bulk_size_m":metadata["requested_mesh_controls"]["bulk_target_size_m"],
        "vertex_count":qc["vertex_count"],"tetrahedron_count":qc["tetrahedron_count"],
        "boundary_facet_count":qc["topology"]["exterior_facet_count"],"boundary_counts":qc["topology"]["boundary_counts"],
        "quality":qc["quality"],"wall_time_s":metadata["wall_time_s"],"peak_rss_kib":metadata["peak_rss_kib"],
        "hard_gate_status":qc["hard_gate_status"],"status":qc["status"],"gpu_used":False,
        "mesh_sha256":qc["mesh_sha256"],"run_id":metadata["run_id"],"configuration_sha256":metadata["configuration_sha256"]}
comparison={"timestamp":timestamp(),"profiles":rows,
    "tetrahedron_ratio_medium_to_coarse":rows["medium"]["tetrahedron_count"]/rows["coarse"]["tetrahedron_count"],
    "mesh_convergence_assessed":False,"interpretation":"Two development meshes share exactly the same frozen boundary; the finer interior target changes interior density only. No physical solution or convergence study."}
write_json(ROOT/"reports/stage01/mesh_comparison.json",comparison)
print(json.dumps({k:{key:value for key,value in v.items() if key!='quality'} for k,v in rows.items()},indent=2))
