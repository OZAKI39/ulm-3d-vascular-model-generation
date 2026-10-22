"""Inspect, extract, and optionally invoke the existing VascularMD model pipeline."""
import logging
from pathlib import Path
import traceback

import numpy as np
import pyvista as pv

from . import qc
from .mevo_export import export_roi, source_polydata, write_csv, write_json
from .mevo_graph import graph_counts, load_exact_graph
from .mevo_roi import extract_graph
from .pipeline import Options, process_file, provenance
from .roi_landmarks import LandmarkError, landmark_template, load_landmarks, sha256, write_landmarks

LOG = logging.getLogger(__name__)


def inspect_file(source: Path, output: Path, *, units="mm", subject_id=None):
    source, output = source.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "FAILED", "source_file": str(source), "warnings": [], "outputs": {}}
    try:
        digest = sha256(source)
        exact = load_exact_graph(source)
        rows = []
        for node, data in exact.topo.nodes(data=True):
            original = int(data["full_id"])
            x, y, z, radius = map(float, data["coords"])
            rows.append({"original_swc_id": original, "topo_node_id": node, "x": x, "y": y, "z": z,
                "radius": radius, "in_degree": exact.topo.in_degree(node), "out_degree": exact.topo.out_degree(node),
                "topology_type": data["type"], "is_root": exact.topo.in_degree(node) == 0,
                "is_terminal": exact.topo.out_degree(node) == 0, "is_bifurcation": exact.topo.out_degree(node) > 1,
                "swc_type": exact.graph.nodes[original]["swc_type"], "component_id": exact.locations[original]["component_id"]})
        csv = output / f"{source.stem}_mevo_topology_nodes.csv"
        mapping = output / f"{source.stem}_mevo_full_topo_mapping.csv"
        vtp = output / f"{source.stem}_mevo_topology.vtp"
        template = output / f"{source.stem}_mevo_landmarks.template.yaml"
        write_csv(csv, rows)
        write_csv(mapping, list(exact.locations.values()))
        mesh = source_polydata(exact)
        mesh.field_data["units"] = [units]
        mesh.save(vtp)
        write_landmarks(template, landmark_template(source, template, units=units, subject_id=subject_id))
        report.update(source_sha256=digest, units=units, source=graph_counts(exact.graph),
                      topology_node_count=len(exact.topo), topology_branch_count=exact.topo.number_of_edges(),
                      vascularmd_check_full_graph=exact.native_valid, upstream=provenance(),
                      swc_type_values=sorted({int(d["swc_type"]) for _, d in exact.graph.nodes(data=True)}),
                      automatic_anatomical_classification_performed=False,
                      outputs={"topology_nodes": csv.name, "full_topo_mapping": mapping.name,
                               "topology_vtp": vtp.name, "landmark_template": template.name})
        if not exact.native_valid:
            report["warnings"].append("Native check_full_graph rejects some source degrees; inspection/exact extraction remain separate from modelability")
        report["warnings"].append("M2/M3, A2/A3 and P2/P3 labels are not inferred from SWC TYPE, branch order, coordinates or diameter; explicit anatomical landmarks are required")
        report["original_file_preserved"] = sha256(source) == digest
        if not report["original_file_preserved"]:
            raise ValueError("Source changed during inspection")
        report["status"] = "PASS_WITH_WARNINGS"
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"
        report["traceback"] = traceback.format_exc()
    write_json(output / f"{source.stem}_mevo_inspect_qc.json", report)
    return report


def extract_file(source: Path, landmarks: Path, output: Path, *, allow_source_mismatch=False,
                 allow_natural_terminal=False, proximal_context="none", model=False, surface=False,
                 accept_native_merges=False, territory=None, side=None):
    source, landmarks, output = source.resolve(), landmarks.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "FAILED", "needs_manual_review": True, "source_file": str(source),
              "landmarks_file": str(landmarks), "warnings": [], "rois": {}, "outputs": {}}
    digest = None
    try:
        digest = sha256(source)
        report.update(source_sha256=digest, landmarks_sha256=sha256(landmarks), upstream=provenance())
        doc, warnings = load_landmarks(landmarks, source, allow_source_mismatch=allow_source_mismatch)
        report["warnings"].extend(warnings)
        for message in warnings:
            LOG.warning(message)
        report["landmark_provenance"] = {k: v for k, v in doc.items() if k != "rois"}
        units = doc["source"].get("units", "mm")
        report["units"] = units
        with (output / "landmarks.used.yaml").open("xb") as stream:
            stream.write(landmarks.read_bytes())
        exact = load_exact_graph(source)
        report["source"] = graph_counts(exact.graph)
        report["topology_node_count"] = len(exact.topo)
        report["topology_branch_count"] = exact.topo.number_of_edges()
        enabled = {name: roi for name, roi in doc["rois"].items() if roi["enabled"]}
        if not enabled:
            raise LandmarkError("No enabled anatomical ROI. 尚未标注：请通过 annotate 或 topology VTP 确认边界后填写 landmark YAML；不会使用空间采样替代 MeVO")
        roi_data = []
        for name, roi in enabled.items():
            item = {"name": name, "side": roi["side"], "territory": roi["territory"], "segments": roi["segments"],
                    "source_file": source.name, "source_sha256": digest,
                    "proximal_nodes": roi["proximal_nodes"], "distal_nodes": roi["distal_nodes"],
                    "excluded_subtrees": roi["exclude_subtree_roots"], "manually_verified": roi["manually_verified"],
                    "proximal_boundary_status": roi["proximal_boundary_status"],
                    "internal_landmarks": roi.get("internal_landmarks", {}), "notes": roi.get("notes", ""),
                    "status": "FAILED", "needs_manual_review": True, "warnings": [], "outputs": {}}
            report["rois"][name] = item
            try:
                if (territory and roi["territory"] != territory) or (side and roi["side"] != side):
                    raise LandmarkError("CLI territory/side disagrees with enabled ROI metadata; no automatic relabeling")
                result = extract_graph(exact, roi, allow_natural_terminal=allow_natural_terminal,
                                       proximal_context=proximal_context)
                item["warnings"] = warnings + result.warnings
                status = "MANUALLY_VERIFIED" if roi["manually_verified"] and not warnings else "UNVERIFIED"
                item["anatomical_status"] = status
                item["anatomical_completeness"] = result.anatomical_completeness
                item["coverage"] = result.coverage
                item.update(graph_counts(result.graph))
                item["statistics"] = qc.graph_statistics(result.graph, item["branch_count"])
                radii = [d["coords"][3] for _, d in result.graph.nodes(data=True)]
                item["statistics"].update(diameter_min=float(min(radii)*2), diameter_median=float(np.median(radii)*2),
                                          diameter_max=float(max(radii)*2))
                item["strict_roi_node_count"] = len(result.graph)
                item["context_node_count"] = len(result.context_nodes)
                item["context_original_node_ids"] = sorted(result.context_nodes)
                item["context_source"] = "original upstream topology edge (possibly partial when proximal lies inside it)" if result.context_nodes else None
                item["context_branch_length"] = float(sum(np.linalg.norm(exact.graph.nodes[a]["coords"][:3] - exact.graph.nodes[b]["coords"][:3])
                    for a, b in (result.modelable_graph.edges if result.modelable_graph is not None else [])
                    if a in result.context_nodes or b in result.context_nodes))
                strict = output / f"{source.stem}_{name}_roi.swc"
                item["outputs"]["strict"] = export_roi(strict, result.graph, exact, name, roi, digest, status)
                modeling_input = strict
                if result.modelable_graph is not None:
                    modeling_input = output / f"{source.stem}_{name}_roi_modelable.swc"
                    item["outputs"]["modelable"] = export_roi(modeling_input, result.modelable_graph, exact, name, roi,
                                                               digest, status, context_nodes=result.context_nodes)
                roi_data.append((name, roi, result))
                if model or surface:
                    modeled = process_file(modeling_input, output / f"{name}_vascularmd", Options(
                        input_units=units, surface=surface, qc_plot=True, accept_native_merges=accept_native_merges))
                    item["model_status"] = modeled["status"]
                    item["model_outputs"] = {k: str(Path(v).relative_to(output)) for k, v in modeled["outputs"].items()}
                    item["model_warnings"] = modeled["warnings"]
                    if modeled["status"] != "success":
                        item["model_failures"] = modeled["failed_components"]
                        raise ValueError("VascularMD modeling failed; exact strict ROI is retained, no fallback surface generated")
                    if surface:
                        vtk = Path(modeled["outputs"]["surface_vtk"])
                        stl = vtk.with_suffix(".stl")
                        pv.read(vtk).triangulate().save(stl)
                        item["model_outputs"]["surface_stl"] = str(stl.relative_to(output))
                item["needs_manual_review"] = status != "MANUALLY_VERIFIED"
                item["status"] = "PASS_WITH_WARNINGS" if item["warnings"] else "PASS"
            except Exception as exc:
                item["status"] = "NEEDS_MANUAL_REVIEW" if isinstance(exc, LandmarkError) else "FAILED"
                item["needs_manual_review"] = True
                item["error"] = f"{type(exc).__name__}: {exc}"
                item["failure_context"] = getattr(exc, "context", {})
                item["traceback"] = traceback.format_exc()
            LOG.info("[%s] Source nodes=%s; topology branches=%s; proximal=%s; distal=%s; excluded=%s; "
                     "extracted nodes=%s; branches=%s; terminals=%s; bifurcations=%s; unbounded=%s; manually_verified=%s; status=%s",
                     name, len(exact.graph), exact.topo.number_of_edges(), roi["proximal_nodes"], roi["distal_nodes"],
                     roi["exclude_subtree_roots"], item.get("node_count"), item.get("branch_count"),
                     item.get("terminal_count"), item.get("bifurcation_count"),
                     item.get("coverage", item.get("failure_context", {}).get("coverage", {})).get("unbounded_path_count"),
                     roi["manually_verified"], item["status"])
        preview = output / f"{source.stem}_mevo_preview.vtp"
        mesh = source_polydata(exact, roi_data)
        mesh.field_data["source_sha256"] = [digest]
        mesh.field_data["units"] = [units]
        mesh.field_data["anatomical_status"] = [f"{name}={item.get('anatomical_status', 'FAILED')}" for name, item in report["rois"].items()]
        mesh.save(preview)
        report["outputs"]["qc_vtp"] = preview.name
        statuses = [v["status"] for v in report["rois"].values()]
        report["status"] = ("FAILED" if "FAILED" in statuses else "NEEDS_MANUAL_REVIEW" if "NEEDS_MANUAL_REVIEW" in statuses
                            else "PASS_WITH_WARNINGS" if "PASS_WITH_WARNINGS" in statuses else "PASS")
        report["needs_manual_review"] = any(v["needs_manual_review"] for v in report["rois"].values())
    except Exception as exc:
        report["status"] = "NEEDS_MANUAL_REVIEW" if isinstance(exc, LandmarkError) else "FAILED"
        report["error"] = f"{type(exc).__name__}: {exc}"
        report["traceback"] = traceback.format_exc()
    report["original_file_preserved"] = digest is not None and sha256(source) == digest
    if digest is not None and not report["original_file_preserved"]:
        report["status"], report["needs_manual_review"] = "FAILED", True
    report["manifest"] = f"{source.stem}_mevo_manifest.json"
    write_json(output / report["manifest"], report)
    return report
