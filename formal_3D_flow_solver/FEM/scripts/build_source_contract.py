#!/usr/bin/env python3
"""Resolve configured lineage and measure existing tagged triangles, without meshing."""
import csv
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import pyvista as pv
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from fem3d.audit import sha256, timestamp, write_json
from fem3d.source import length_to_m, resolve_source_path, validate_contract


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def binary_stl_triangles(path):
    with path.open("rb") as f:
        f.read(80)
        count = int.from_bytes(f.read(4), "little")
        dtype = np.dtype([("normal", "<f4", (3,)), ("vertices", "<f4", (3,3)), ("attribute", "<u2")])
        records = np.fromfile(f, dtype=dtype, count=count)
    if path.stat().st_size != 84 + count*50 or len(records) != count:
        raise ValueError("Not the expected binary STL")
    return records["vertices"]


def build():
    if sys.flags.optimize:
        raise RuntimeError("Run this audit without Python -O so all evidence checks stay enabled")
    config = read_json(ROOT/"configs/stage00.json")
    ref = Path(config["reference_3d"])
    evidence = {}
    mappings = []

    def resolve(value):
        path = resolve_source_path(value, ref, config["allowed_historical_prefix"])
        if str(value) != str(path):
            mappings.append({"recorded":str(value), "wsl":str(path)})
        return path

    def register(path, role):
        path = Path(path).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        evidence[str(path)] = {"path":str(path), "sha256":sha256(path), "role":role}
        return path

    def load(path, role):
        return read_json(register(path, role))

    source_configs = []
    for key in ("flow_config", "surface_config", "preprocess_config"):
        path = register(resolve(config[key]), key)
        source_configs.append(str(path))
    flow, surface_cfg, pre_cfg = [yaml.safe_load(Path(p).read_text()) for p in source_configs]
    run = resolve(flow["paths"]["source_surface_run"])
    schema = config["geometry_schema"]
    loader = register(ref/schema["loader_source"], "Official loader fixes surface and BC filenames")
    for key in ("tagged_surface", "meter_surface", "boundary_conditions", "boundary_manifest"):
        if Path(schema[key]).name not in loader.read_text():
            raise ValueError(f"Configured schema no longer matches official input loader: {key}")
    tagged = register(run/schema["tagged_surface"], "Canonical tag-bearing surface, um")
    meter = register(run/schema["meter_surface"], "CFD SI STL companion, tags absent")
    manifest_path = register(run/schema["boundary_manifest"], "Authoritative port_id/role/entity mapping")
    manifest = rows(manifest_path)
    bc_path = run/schema["boundary_conditions"]
    bc = load(bc_path, "Original BC values; roles reusable, physical values not adopted")
    final_qc_path = run/"qc/final_surface_qc.json"
    final_qc = load(final_qc_path, "Final surface topology and labelled cap metadata")
    summary = load(run/"qc/run_summary.json", "Final surface acceptance and geometry history")
    units_qc = load(run/"qc/meter_scale_qc.json", "Recorded um to m coordinate relation")
    assert units_qc["scale_factor"] == 1e-6 and units_qc["status"] == "PASS"
    assert final_qc["status"] == "PASS"
    frozen = load(run/"input/frozen_open_geometry_reference.json", "Frozen open surface lineage")
    open_run = resolve(frozen["source_run_path"])
    for key, hk in (("source_open_vtp","source_vtp_sha256"),("source_open_stl","source_stl_sha256")):
        p = register(resolve(frozen[key]), "Unmodified frozen open surface")
        assert sha256(p) == frozen[hk], (p, "Frozen geometry hash differs")
    original = load(open_run/"input/original_surface_reference.json", "Original reconstruction hashes and lineage")
    for path, expected in original["immutable_hashes_before"].items():
        p = register(resolve(path), "Read-only frozen lineage evidence")
        assert sha256(p) == expected, (p, "Lineage hash differs")
    frozen_surface_config = register(run/"input/cfd_surface_prepare.yaml", "Frozen surface configuration")
    source_configs.append(str(frozen_surface_config))
    frozen_cfg = yaml.safe_load(frozen_surface_config.read_text())
    pre_run = resolve(frozen_cfg["paths"]["cfd_preprocess_run"])
    assert pre_run == resolve(surface_cfg["paths"]["cfd_preprocess_run"])
    pre_manifest = load(pre_run/"input/input_manifest.json", "Actual preprocess run resolution; not reselecting newest")
    geometry_ref = load(pre_run/"input/geometry_reference.json", "Reconstruction reference")
    model_run = resolve(pre_manifest["model_run"])
    assert model_run == resolve(geometry_ref["run_root"])
    assert model_run == resolve(original["geometry_reference"]["run_root"])
    sampling_run = resolve(pre_manifest["sampling_run"])
    assert sampling_run == resolve(pre_cfg["paths"]["sampling_run"])
    model_meta = load(model_run/"input/metadata.json", "ROI identity, units and canonical SWC SHA256")
    canonical = register(model_run/"input/roi_core.swc", "Canonical local ROI SWC")
    assert sha256(canonical) == model_meta["canonical_swc_sha256"]
    assert model_meta["coordinate_unit"] == "um" and model_meta["radius_unit"] == "um"
    for name, value in geometry_ref["files"].items():
        p = register(resolve(value), "Ultraliser source geometry")
        assert sha256(p) == geometry_ref["sha256"][name]
    roi_id = pre_manifest["roi_id"]
    selected = register(sampling_run/"manifests/selected_rois.csv", "Selected ROI identity, global nodes and edges")
    matching = [r for r in rows(selected) if r["roi_id"] == roi_id]
    assert len(matching) == 1 and int(matching[0]["anchor_id"]) == model_meta["anchor_id"]
    roi_path = register(sampling_run/"roi_library"/f"{roi_id}.npz", "ROI archive schema from sampling_io.write_roi_library")
    with np.load(roi_path, allow_pickle=False) as archive:
        assert str(archive["roi_id"]) == roi_id
        assert str(archive["source_model_id"]) == pre_manifest["source_model_id"]
        roi_box = [archive["bbox_min_um"].tolist(),archive["bbox_max_um"].tolist()]
    rodent_run = resolve(pre_manifest["rodent_run"])
    rodent = load(rodent_run/"normalized_index.json", "SWC source and normalized analysis data")
    samples = [s for s in rodent["samples"] if s["record"]["sample_id"] == pre_manifest["source_model_id"]]
    assert len(samples) == 1
    sample = samples[0]
    source_swc = register(resolve(sample["record"]["swc_path"]), "Original source SWC")
    analysis_swc = register(resolve(sample["analysis_swc_text_path"]), "Validated single component SWC")
    for key in ("normalized_swc_path","reference_normalized_swc_path"):
        register(resolve(sample[key]), key)
    pre_ports_path = register(pre_run/"roi/port_classification.csv", "Original ROI planes/normals/areas and port roles")
    pre_ports = {p["port_id"]:p for p in rows(pre_ports_path)}
    register(pre_run/"roi/port_planes.vtp", "Original ROI plane visualization, not distal FEM cap")
    plane_contract_path = register(resolve(flow["paths"]["physical_plane_contract"]), "Interior LBM measurement planes, not FEM boundaries")
    plane_contract = read_json(plane_contract_path)
    # Historical completed flow runs are reported by their metadata, not guessed from names.
    production = []
    for path in (ref/flow["paths"]["output_root"]).glob("*/qc/run_summary.json"):
        data = read_json(path)
        if "geometry_source" in data and resolve(data["geometry_source"]) == run and "PRODUCTION" in data.get("status", ""):
            register(path, "Completed production evidence matching configured source geometry")
            production.append({"metadata_path":str(path), "run":str(path.parents[1]),
                               "status":data["status"], "execution_mode":data.get("execution_mode"),
                               "fresh_full_production_steady_solve":data.get("fresh_full_production_steady_solve")})
    surf = pv.read(tagged)
    raw_faces = np.asarray(surf.faces).reshape(-1,4)
    assert (raw_faces[:,0] == 3).all(), "Source must already be triangulated"
    points = np.asarray(surf.points, dtype=np.float64)
    faces = raw_faces[:,1:]
    tri = points[faces]
    tags = np.asarray(surf.cell_data["CellEntityIds"]).astype(int)
    cross = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
    areas = np.linalg.norm(cross,axis=1)/2
    assert np.isfinite(points).all() and (areas > 0).all()
    center = points.mean(axis=0)
    signed_volume = float(np.einsum("ij,ij->i",tri[:,0]-center,cross).sum()/6)
    assert abs(signed_volume) > 0
    direction = float(np.sign(signed_volume))
    # The source records literal STL triangle sequence equality after float32 scaling.
    meter_tri = binary_stl_triangles(meter)
    scale_equal = np.array_equal(length_to_m(tri,"um").astype(np.float32), meter_tri)
    assert scale_equal, "Tagged VTP and meter STL do not describe identical triangle sequences"
    clean = surf.clean(tolerance=0.0, absolute=True)
    edges = clean.extract_feature_edges(boundary_edges=True, non_manifold_edges=True,
                                       feature_edges=False, manifold_edges=False)
    assert edges.n_cells == 0, "Source surface has open or nonmanifold edges"
    region = clean.connectivity()
    assert len(np.unique(region.cell_data["RegionId"])) == 1
    wall_tag = int(final_qc["boundary_mapping"]["wall_entity_id"])
    assert set(tags) == {wall_tag} | {int(p["vmtk_cap_entity_id"]) for p in manifest}
    boundary_qc = {p["port_id"]:p for p in final_qc["boundary_mapping"]["boundaries"]}
    bc_roles = {bc["inlet"]["port_id"]:bc["inlet"]["role"], **{p["port_id"]:p["role"] for p in bc["outlets"]}}
    ports = []
    for row in manifest:
        pid = row["port_id"]
        assert row["role"] == bc_roles[pid] == pre_ports[pid]["role"] == boundary_qc[pid]["role"]
        assert row["role"] in ("ASSUMED_INLET","ASSUMED_OUTLET")
        tag = int(row["vmtk_cap_entity_id"])
        mask = tags == tag
        assert int(mask.sum()) == int(row["triangle_count"]) == int(boundary_qc[pid]["triangle_count"])
        area = float(areas[mask].sum())
        assert np.isclose(area,float(row["area_um2"]),rtol=1e-6,atol=0)
        centroid = np.average(tri[mask].mean(axis=1),axis=0,weights=areas[mask])
        normal = cross[mask].sum(axis=0)*direction
        normal /= np.linalg.norm(normal)
        deviation = float(np.abs((tri[mask]-centroid)@normal).max())
        assert deviation <= 1e-4, "Distal cap is not planar to float32 source precision"
        pre = pre_ports[pid]
        pre_n = np.array([float(pre["outward_normal_"+axis]) for axis in "xyz"])
        # This check corroborates sign only. It never assigns roles by coordinate direction.
        alignment = float(normal@pre_n/np.linalg.norm(pre_n))
        assert alignment > .999, "Cap normal disagrees with the same port's original metadata normal"
        matching_labels = [k for k,v in plane_contract["ports"].items() if v["source_port_id"] == pid]
        assert len(matching_labels) == 1
        name = matching_labels[0]
        interior_dot = float(np.dot(normal, plane_contract["ports"][name]["planes"]["central"]["unit_normal"]))
        ports.append({"name":name, "port_id":pid, "role":"inlet" if row["role"] == "ASSUMED_INLET" else "outlet",
                      "source_role":row["role"], "boundary_origin":row["boundary_origin"],
                      "surface_entity_id":tag,"triangle_count":int(mask.sum()),
                      "area_m2":area*1e-12,"source_area_um2":float(row["area_um2"]),
                      "plane_origin_m":length_to_m(centroid,"um").tolist(), "outward_normal":normal.tolist(),
                      "plane_equation":"dot(x_m - plane_origin_m, outward_normal) = 0",
                      "max_plane_deviation_m":deviation*1e-6,
                      "normal_source":"Area-weighted oriented triangles for the metadata-assigned CellEntityIds; outward sign from closed signed volume",
                      "role_source":str(manifest_path), "geometry_source":str(tagged),
                      "upstream_role_source":str(pre_ports_path), "upstream_plane_origin_m":[float(pre[a+"_um"])*1e-6 for a in "xyz"],
                      "upstream_outward_normal":pre_n.tolist(), "upstream_normal_dot":alignment,
                      "legacy_interior_plane_normal_dot_cap_normal":interior_dot,
                      "legacy_interior_plane_role":"Interior flux measurement only; its normal is not the FEM cap normal",
                      "source_patch_path":str(register(resolve(row["stl_path"]),"Legacy cap patch for cross-reference")),
                      "normal_recomputed_from_tagged_cap":True})
    contract = {"schema_version":1,"status":"PASS","created_utc":timestamp(),
                "geometry_path":str(tagged),"geometry_units":"um","geometry_sha256":sha256(tagged),
                "geometry_role":"Canonical existing tagged closed surface; inlet/outlet caps are artificial boundary patches",
                "meter_geometry_path":str(meter),"meter_geometry_units":"m","meter_geometry_sha256":sha256(meter),
                "solver_units":"SI","length_conversion_to_m":1e-6,"area_conversion_to_m2":1e-12,
                "input_adapter":"src/fem3d/source.py:length_to_m; never pass raw micrometer coordinates to solver",
                "source_run":str(run),"selection_evidence":str(ref/config["flow_config"])+"#paths.source_surface_run",
                "current_flow_mode":flow["execution"]["mode"],"completed_production_evidence":production,
                "inlets":[p for p in ports if p["role"] == "inlet"],"outlets":[p for p in ports if p["role"] == "outlet"],
                "wall_definition":f"All and only CellEntityIds={wall_tag}; includes core and artificial extension walls; future PDMS no-slip",
                "wall_entity_id":wall_tag,"tag_array":"CellEntityIds","wall_triangle_count":int((tags==wall_tag).sum()),
                "source_metadata_files":[k for k in evidence if Path(k).suffix in (".json",".csv")],
                "source_configs":source_configs,"source_files":list(evidence.values()),
                "historical_path_mappings":mappings,"simulation_direction":bc["simulation_direction"],
                "lineage":{"source_swc":str(source_swc),"analysis_swc":str(analysis_swc),
                           "rodent_run":str(rodent_run),"sampling_run":str(sampling_run),"roi_archive":str(roi_path),
                           "roi_id":roi_id,"roi_anchor":model_meta["anchor_id"],"roi_bbox_um":roi_box,
                           "canonical_roi_swc":str(canonical),"reconstruction_run":str(model_run),
                           "radius_scale":model_meta["radius_scale_for_ultraliser"],"preprocess_run":str(pre_run),
                           "frozen_open_run":str(open_run)},
                "geometry_checks":{"triangles":len(faces),"points_in_file":len(points),"bounds_um":[points.min(axis=0).tolist(),points.max(axis=0).tolist()],
                                   "signed_volume_um3":signed_volume,"meter_stl_matches_scaled_tagged_triangles":scale_equal,
                                   "boundary_or_nonmanifold_edges_after_exact_vertex_merge":edges.n_cells,"connected_components":1},
                "legacy_physics_not_adopted":{"preprocess_inlet_Q_m3_s":bc["inlet"]["flow_rate_m3_s"],
                                              "current_lbm_inlet_Q_m3_s":flow["boundary_conditions"]["target_volume_flow_m3_s"],
                                              "outlet_pressures_pa":flow["boundary_conditions"]["outlet_gauge_pressures_pa"],
                                              "experimental_Q_m3_s":None,"reason":"User has specified BC type, not an experimental Q value; old parabolic profile and pressure corrections are not FEM boundary conditions"},
                "warnings":["Roles are explicit structural assumptions, not measured physiological directions.",
                            "Artificial extensions and original ROI planes differ from distal caps.",
                            "Upstream final surface status includes PENDING_MANUAL_REVIEW; no human approval is inferred.",
                            "Completed production evidence is a validated replay, not a new FRESH_STEADY solution."]}
    validate_contract(contract)
    write_json(ROOT/"reports/stage00/source_contract.json",contract)
    # Copy only two necessary surfaces and small defining metadata into this project.
    staging = ROOT/"inputs/stage00"
    staging.mkdir(parents=True,exist_ok=True)
    copies = []
    for src, name in ((tagged,"source_surface_um.vtp"),(meter,"source_surface_m.stl"),
                      (manifest_path,"boundary_manifest.csv"),(bc_path,"legacy_boundary_conditions.json")):
        dest = staging/name
        shutil.copy2(src,dest)
        copies.append({"original_path":str(src),"project_path":str(dest.relative_to(ROOT)),"sha256":sha256(dest)})
    shutil.copy2(ROOT/"reports/stage00/source_contract.json",staging/"source_contract.json")
    write_json(staging/"input_manifest.json",{"source_of_truth":str(ROOT),"copies":copies})
    print(json.dumps({"status":"PASS","geometry":str(tagged),"inlets":len(contract['inlets']),"outlets":len(contract['outlets']),"source_files":len(evidence)},indent=2))


if __name__ == "__main__":
    try:
        build()
    except Exception as exc:
        write_json(ROOT/"reports/stage00/source_contract.json",{"status":"FAIL","error":repr(exc),"created_utc":timestamp()})
        raise
