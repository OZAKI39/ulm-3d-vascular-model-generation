"""One native load/model/sample/export per input; explicit failure reports."""
from __future__ import annotations

import contextlib
import hashlib
import importlib.metadata
import json
import logging
import platform
import subprocess
import time
import traceback
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path

import networkx as nx
import numpy as np
from threadpoolctl import threadpool_limits

from . import qc
from .branch_aic import JUNCTION_POLICY as BRANCH_JUNCTION_POLICY, fit_branches
from .swc_export import read_source, validate_tree, write_swc
from .vascularmd_adapter import COMMIT, JUNCTION_POLICY, TopologyChanged, cache_native_distance_points, guard_exhausted_apex_search, load_tree, prepare_export_topology, sample_model

LOG = logging.getLogger(__name__)


def provenance():
    root = Path(__file__).resolve().parents[1]
    vendor = root / "third_party/vascularmd"
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=vendor, text=True, capture_output=True, check=True).stdout.strip()
    if commit != COMMIT:
        raise ValueError(f"VascularMD checkout differs from reviewed commit: {commit}")
    return {"repository": "https://github.com/megdec/vascularmd", "commit": commit,
            "source_files_sha256": {name: hashlib.sha256((vendor / name).read_bytes()).hexdigest()
                                    for name in ("ArterialTree.py", "Spline.py", "Model.py", "Nfurcation.py", "utils.py")}}


@dataclass
class Options:
    sample_mode: str = "original-count"
    spacing_mm: float | None = None
    auto_resample: bool | None = None
    surface: bool = True
    circumferential_n: int = 24
    longitudinal_density: float = 0.2
    qc_plot: bool = False
    plot_branches: int = 6
    blas_threads: int = 1
    accept_native_merges: bool = False
    input_units: str = "mm"
    model_mode: str = "network"

    def validate(self):
        if self.model_mode not in {"network", "branches"}:
            raise ValueError("model_mode must be network or branches")
        if self.model_mode == "branches" and (self.surface or self.sample_mode != "original-count" or self.accept_native_merges):
            raise ValueError("branches mode requires no surface, original-count sampling and no native merges")
        if self.input_units not in {"mm", "um"}:
            raise ValueError("input_units must be mm or um; XYZ and radius must share this physical unit")
        if self.input_units != "mm" and self.spacing_mm is not None:
            raise ValueError("--spacing-mm requires --input-units mm; use original-count for micrometre inputs")
        if self.sample_mode not in ("original-count", "spacing"):
            raise ValueError("Unsupported sample mode")
        if self.sample_mode == "spacing" and (self.spacing_mm is None or not np.isfinite(self.spacing_mm) or self.spacing_mm <= 0):
            raise ValueError("spacing mode requires a finite positive --spacing-mm")
        if self.sample_mode != "spacing" and self.spacing_mm is not None:
            raise ValueError("--spacing-mm requires --sample-mode spacing")
        if self.circumferential_n < 8 or self.circumferential_n % 4:
            raise ValueError("--circumferential-n must be a multiple of four, at least eight")
        if not np.isfinite(self.longitudinal_density) or self.longitudinal_density <= 0:
            raise ValueError("--longitudinal-density must be positive and finite")
        if self.plot_branches < 1 or self.blas_threads < 1:
            raise ValueError("Plot count and BLAS thread count must be positive")


def save_native_surface(tree, path: Path, options: Options):
    tree.compute_cross_sections(options.circumferential_n, options.longitudinal_density, parallel=False)
    mesh = tree.mesh_surface()
    if mesh.n_cells == 0 or not np.isfinite(mesh.points).all():
        raise ValueError("Native surface is empty or contains nonfinite points")
    mesh.save(str(path))
    return {"points": mesh.n_points, "cells": mesh.n_cells, "open_edges": mesh.n_open_edges,
            "post_filter": "none added; native Nfurcation mesh relaxation/projection is retained",
            "cfd_ready_certified": False}


def process_file(input_path: Path, output_dir: Path, options: Options) -> dict:
    options.validate()
    input_path = input_path.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = input_path.stem + "_vmd_"
    paths = {"swc": output_dir / (prefix + "smooth.swc"), "surface": output_dir / (prefix + "surface.vtk"),
             "qc": output_dir / (prefix + "qc.json"), "radius": output_dir / (prefix + "radius.csv"),
             "plot": output_dir / (prefix + "radius.png"), "log": output_dir / (prefix + "run.log")}
    existing = [str(p) for p in paths.values() if p.exists()]
    if existing:
        raise FileExistsError(f"Output already exists; choose a new output directory: {existing}")
    start = time.monotonic()
    source_hash = hashlib.sha256(input_path.read_bytes()).hexdigest()
    report = {"status": "failed", "needs_manual_review": True, "source": str(input_path), "source_sha256": source_hash,
              "vascularmd_commit": COMMIT, "options": asdict(options), "warnings": [], "failed_components": [],
              "units": f"XYZ and radius in {options.input_units}; original numerical units preserved; no scaling",
              "bifurcation_swc_policy": BRANCH_JUNCTION_POLICY if options.model_mode == "branches" else JUNCTION_POLICY,
              "native_model_call": {"radius_model": True, "criterion": "AIC", "akaike": False,
                                    "note": "criterion chooses native smoothing lambdas; akaike controls a different control-point-count selection option"},
              "versions": {"python": platform.python_version(), **{p: importlib.metadata.version(p) for p in ("numpy", "scipy", "pyvista", "vtk", "networkx", "geomdl", "threadpoolctl")}},
              "outputs": {}}
    tree, source, branches = None, None, None
    stage = "input_validation"
    paths["log"].touch(exist_ok=False)
    file_handler = logging.FileHandler(paths["log"], mode="a", encoding="utf-8")
    logging.getLogger().addHandler(file_handler)
    captured = []
    try:
        with paths["log"].open("a") as native_log, contextlib.redirect_stdout(native_log), warnings.catch_warnings(record=True) as captured, threadpool_limits(limits=options.blas_threads):
            warnings.simplefilter("always")
            source = read_source(input_path)
            report["upstream"] = provenance()
            stage = "native_load_and_preprocessing"
            tree, original_topo, branches = load_tree(source, options.auto_resample)
            report["preprocessing"] = {"automatic_resampling": options.auto_resample is not False,
                                       "original_node_count": len(source.graph), "vascularmd_preprocessing_node_count": len(tree.get_full_graph()),
                                       "official_density_range": [0.4, 0.6], "density_unit": f"points/{options.input_units}",
                                       "minimum_points_when_resampled": 4}
            report["input"] = qc.graph_statistics(source.graph, len(branches))
            LOG.info("%s: roots=%s terminals=%s bifurcations=%s branches=%s; nodes %s -> %s after native preprocessing",
                     input_path.name, report["input"]["root_count"], report["input"]["terminal_count"], report["input"]["bifurcation_count"],
                     len(branches), len(source.graph), len(tree.get_full_graph()))
            stage = "native_model"
            raw_branches = branches
            if options.model_mode == "branches":
                report["native_model_call"].update({"entrypoint": "Spline.approximation per original branch",
                    "end_constraint": [True, False, False, True], "max_distance": 6,
                    "nfurcation_reconstruction": False})
                with cache_native_distance_points():
                    graph, report["branch_fits"] = fit_branches(tree, branches, source)
                topology = {"policy": "fixed-original-branch-endpoints", "input_topology_preserved": True,
                            "merged_bifurcation_swc_ids": [], "merge_groups": [],
                            "all_original_terminals_preserved": True, "original_branch_count": len(branches),
                            "native_branch_count": len(branches), "fixed_endpoint_xyzr_preserved": True}
            else:
                with cache_native_distance_points(), guard_exhausted_apex_search():
                    tree.model_network(radius_model=True, criterion="AIC")
                stage = "model_sampling_and_topology_validation"
                export_topo, branches, topology = prepare_export_topology(tree, original_topo, branches, source, options.accept_native_merges)
                graph = sample_model(tree, export_topo, branches, source, options.sample_mode, options.spacing_mm)
            report["model_seconds"] = time.monotonic() - start
            graph.graph["topology_policy"] = topology
            report["comparison"] = qc.compare(branches, source.graph, graph, raw_branches=raw_branches)
            report["topology"] = topology
            report["topology_preserved"] = topology["input_topology_preserved"]
            report["native_model_topology_preserved"] = True
            if topology["merged_bifurcation_swc_ids"]:
                report["warnings"].append(f"Explicitly accepted native bifurcation merges: {topology['merge_groups']}; all original terminals preserved.")
            report["warnings"].append(report["bifurcation_swc_policy"])
            report["warnings"].append("Input parent direction is treated as flow direction; SWC alone cannot confirm measured physiological flow.")
            mixed = [list(b.key) for b in branches if len(set(b.types[1:-1])) > 1]
            report["type_mapping"] = {"policy": "original node TYPE retained at the same ordered branch index" if options.model_mode == "branches" else "original junction/end types; internal type from nearest original normalized arc position inside the same branch", "mixed_type_branches": mixed}
            if mixed and options.model_mode != "branches":
                report["warnings"].append(f"Mixed original TYPE labels within branches {mixed}: nearest arc-position inheritance used; no artery classification inferred.")
            stage = "swc_export"
            mapping = write_swc(paths["swc"], graph, input_path, COMMIT, options.sample_mode)
            if options.model_mode == "branches":
                report["source_to_output_node_ids"] = mapping
            round_trip = read_source(paths["swc"])
            validate_tree(round_trip.graph)
            if len(round_trip.graph) != len(graph) or set(round_trip.graph.edges) != {(mapping[a], mapping[b]) for a, b in graph.edges}:
                raise ValueError("SWC round-trip topology differs from sampled graph")
            for node, output_id in mapping.items():
                np.testing.assert_array_equal(graph.nodes[node]["coords"], round_trip.graph.nodes[output_id]["coords"])
                if graph.nodes[node]["swc_type"] != round_trip.graph.nodes[output_id]["swc_type"]:
                    raise ValueError(f"SWC type changed during serialization at node {node}")
            report["round_trip_valid"] = True
            report["outputs"]["swc"] = str(paths["swc"].resolve())
            qc.write_radius_csv(paths["radius"], branches)
            report["outputs"]["radius_csv"] = str(paths["radius"].resolve())
            if options.qc_plot:
                qc.plot_diagnostics(paths["plot"], branches, options.plot_branches, input_units=options.input_units,
                                    bifurcation_label="Near fixed original bifurcation" if options.model_mode == "branches" else "Near native bifurcation")
                report["outputs"]["qc_plot"] = str(paths["plot"].resolve())
            if options.surface:
                stage = "native_surface"
                report["surface"] = save_native_surface(tree, paths["surface"], options)
                report["outputs"]["surface_vtk"] = str(paths["surface"].resolve())
            report["status"], report["needs_manual_review"] = "success", False
    except Exception as exc:
        report["status"] = "partial" if "swc" in report["outputs"] else "failed"
        failure = {"stage": stage, "error_type": type(exc).__name__, "message": str(exc),
                   "component": getattr(exc, "context", getattr(tree, "active_component", None)), "traceback": traceback.format_exc()}
        report["failed_components"].append(failure)
        if tree is not None and tree._model_graph is not None:
            report["completed_model_edges"] = [list(e) for e in tree.get_model_graph().edges if tree.get_model_graph().edges[e].get("spline") is not None]
        LOG.error("%s: %s; context=%s", input_path.name, exc, failure["component"])
        # Native bifurcation merging must NEVER be disguised as topology-preserving SWC.
        # Still retain its native surface for inspection when that stage is usable.
        if isinstance(exc, TopologyChanged) and options.surface:
            report["topology_preserved"] = False
            report["warnings"].append("Native topology changed. No smooth SWC exported; any VTK is diagnostic only and needs manual review.")
            surface_warnings = []
            try:
                with paths["log"].open("a") as native_log, contextlib.redirect_stdout(native_log), threadpool_limits(limits=options.blas_threads), warnings.catch_warnings(record=True) as surface_warnings:
                    warnings.simplefilter("always")
                    report["surface"] = save_native_surface(tree, paths["surface"], options)
                report["surface"]["diagnostic_only"] = True
                report["outputs"]["surface_vtk"] = str(paths["surface"].resolve())
                report["status"] = "partial"
            except Exception as surface_error:
                report["failed_components"].append({"stage": "diagnostic_surface", "error_type": type(surface_error).__name__,
                                                     "message": str(surface_error), "component": tree.active_component,
                                                     "traceback": traceback.format_exc()})
            finally:
                captured.extend(surface_warnings)
    finally:
        report["warnings"].extend(sorted({f"{w.category.__name__}: {w.message} ({Path(w.filename).name}:{w.lineno})" for w in captured}))
        report["original_file_preserved"] = hashlib.sha256(input_path.read_bytes()).hexdigest() == source_hash
        if not report["original_file_preserved"]:
            report["status"], report["needs_manual_review"] = "failed", True
        report["elapsed_seconds"] = time.monotonic() - start
        with paths["qc"].open("x", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, ensure_ascii=False, allow_nan=False)
            handle.write("\n")
        logging.getLogger().removeHandler(file_handler)
        file_handler.close()
    LOG.info("%s: %s; QC=%s", input_path.name, report["status"], paths["qc"])
    return report
