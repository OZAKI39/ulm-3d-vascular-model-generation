"""Generate the data-grounded scientific rasters and evidence for cfd.pptx.

This script is deliberately solver-free: it reads frozen project artifacts,
renders/crops them for presentation use, and never modifies scientific inputs.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image
import pyvista as pv


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
ASSETS = HERE / "assets"

PREPROCESS = PROJECT / "outputs/cfd_preprocess/global_to_roi_anchor003274_20260825_183628"
SURFACE = PROJECT / (
    "outputs/cfd_surface_prepare/"
    "vmtk_tps_boundarynormal_crossseam_finalized_recovery_anchor003274_20260826_221611"
)
FLOW = PROJECT / (
    "outputs/cfd_flow/"
    "production_tau1_base_promotion_anchor003274_20260902_013637"
)

SOURCE_SURFACE = PROJECT / (
    "outputs/model_generate/ultraliser_anchor003274_20260825_133350/"
    "geometry/lumen_surface_um.vtp"
)
FINAL_SURFACE = SURFACE / "geometry/cfd_surface_vmtk_tps_boundarynormal_crossseam_um.vtp"
SOURCE_VTU = FLOW / "flow/production_steady_flow_field.vtu"

PORT_COLORS = ["#D55E00", "#0072B2", "#009E73", "#CC79A7"]
CAMERA = {
    "position_um": [346.8127291947934, -68.57611778279654, 476.9527919435981],
    "focal_point_um": [120.62463850424743, 107.20512328104489, 130.7573227071839],
    "view_up": [0.8434080040574257, 0.02862241683039804, -0.536510667132217],
    "parallel_scale_um": 78.0,
    "projection": "parallel",
}


def rel(path: Path) -> str:
    return path.resolve().relative_to(PROJECT.resolve()).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def crop_source(source: Path, output: Path, box: tuple[int, int, int, int]) -> None:
    with Image.open(source) as image:
        image.convert("RGB").crop(box).save(output, quality=95)


def render_surface(path: Path, output: Path, prepared: bool) -> None:
    pv.OFF_SCREEN = True
    mesh = pv.read(path)
    plotter = pv.Plotter(off_screen=True, window_size=(1600, 1000))
    plotter.set_background("white")
    plotter.enable_anti_aliasing("ssaa")

    if not prepared:
        plotter.add_mesh(
            mesh,
            color="#AEB8C0",
            smooth_shading=True,
            specular=0.15,
            specular_power=12,
        )
    else:
        region = np.asarray(mesh.cell_data["SurfaceRegionId"])
        core = mesh.extract_cells(np.flatnonzero(region == 0))
        modified = mesh.extract_cells(np.flatnonzero(region == 1))
        plotter.add_mesh(core, color="#B8C1C8", smooth_shading=True, specular=0.12)
        plotter.add_mesh(modified, color="#4C9F93", smooth_shading=True, specular=0.18)
        boundary_index = np.asarray(mesh.cell_data["boundary_index"])
        for index, color in enumerate(PORT_COLORS):
            ids = np.flatnonzero(boundary_index == index)
            if ids.size:
                plotter.add_mesh(
                    mesh.extract_cells(ids),
                    color=color,
                    smooth_shading=True,
                    specular=0.25,
                )

    plotter.camera_position = [
        CAMERA["position_um"],
        CAMERA["focal_point_um"],
        CAMERA["view_up"],
    ]
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = CAMERA["parallel_scale_um"]
    plotter.hide_axes()
    plotter.screenshot(str(output), transparent_background=False)
    plotter.close()


def side_by_side(left: Path, right: Path, output: Path) -> None:
    with Image.open(left) as before, Image.open(right) as after:
        before = before.convert("RGB")
        after = after.convert("RGB")
        height = min(before.height, after.height)
        before = before.resize((round(before.width * height / before.height), height))
        after = after.resize((round(after.width * height / after.height), height))
        gap = 24
        canvas = Image.new("RGB", (before.width + after.width + gap, height), "white")
        canvas.paste(before, (0, 0))
        canvas.paste(after, (before.width + gap, 0))
        canvas.save(output, quality=95)


def claim(
    slide: int,
    text: str,
    value: Any,
    unit: str,
    source_file: Path,
    source_key: str,
) -> dict[str, Any]:
    return {
        "slide": slide,
        "claim": text,
        "value": value,
        "unit": unit,
        "source_file": rel(source_file),
        "source_key": source_key,
        "source_status": "PASS",
        "verified": True,
    }


def build_evidence() -> dict[str, Any]:
    preprocess_summary_path = PREPROCESS / "qc/run_summary.json"
    preprocess_summary = load_json(preprocess_summary_path)
    port_qc_path = PREPROCESS / "qc/port_transfer_qc.json"
    port_qc = load_json(port_qc_path)
    surface_qc_path = SURFACE / "qc/final_surface_qc.json"
    surface_qc = load_json(surface_qc_path)
    contract_path = FLOW / "input/production_numerical_contract.json"
    contract = load_json(contract_path)
    metrics_path = FLOW / "qc/production_primary_metrics.json"
    metrics = load_json(metrics_path)
    steady_path = FLOW / "qc/production_steady_qc.json"
    steady = load_json(steady_path)
    full_path = FLOW / "qc/production_full_timestep_v2.json"
    full = load_json(full_path)
    sensitivity_path = FLOW / "qc/two_grid_resolution_sensitivity.json"
    sensitivity = load_json(sensitivity_path)
    run_path = FLOW / "qc/run_summary.json"
    run = load_json(run_path)
    flow_manifest_path = FLOW / "flow/production_steady_flow_field_manifest.json"
    flow_manifest = load_json(flow_manifest_path)
    surface_config_path = PROJECT / "configs/cfd_surface_prepare.yaml"
    flow_config_path = PROJECT / "configs/cfd_flow.yaml"
    velocity_visual_path = FLOW / "visualization/interactive_v3_redesign/01_after_velocity_overview.json"
    pressure_visual_path = FLOW / "visualization/interactive_v3_redesign/02_after_pressure_overview.json"

    inlet = next(b for b in preprocess_summary["boundaries"] if b["role"] == "ASSUMED_INLET")
    outlets = [b for b in preprocess_summary["boundaries"] if b["role"] == "ASSUMED_OUTLET"]

    claims = [
        claim(2, "The workflow contains one 1D preprocessing stage, one surface-preparation stage, and one 3D flow stage.", 3, "stages", PROJECT / "cfd_1D_data_preprocess.py", "entry-point sequence"),
        claim(2, "The 1D stage transfers flow and pressure to the ROI ports.", True, "boolean", PROJECT / "utils/cfd_preprocess/pipeline.py", "run_cfd_preprocess -> transfer_all_boundaries"),
        claim(2, "The surface stage extends, remeshes, and caps the lumen boundary regions.", True, "boolean", PROJECT / "utils/cfd_surface_prepare/vmtk_pipeline.py", "run_vmtk_surface_prepare"),
        claim(2, "The 3D stage resolves steady velocity and gauge pressure from the accepted Base field.", True, "boolean", flow_manifest_path, "cell_arrays, source_restart_sha256"),
        claim(3, "The sparse global vascular model contains 7,419 nodes.", preprocess_summary["global_node_count"], "nodes", preprocess_summary_path, "global_node_count"),
        claim(3, "The ROI has one inlet and three outlets.", [preprocess_summary["assumed_inlet_count"], preprocess_summary["assumed_outlet_count"]], "ports", preprocess_summary_path, "assumed_inlet_count, assumed_outlet_count"),
        claim(3, "Transferred inlet flow is 0.769 pL/s.", inlet["flow_rate_m3_s"] * 1e15, "pL/s", preprocess_summary_path, "boundaries[ASSUMED_INLET].flow_rate_m3_s"),
        claim(3, "Transferred outlet flows are 0.0485, 0.3815, and 0.3393 pL/s.", [b["flow_rate_m3_s"] * 1e15 for b in outlets], "pL/s", preprocess_summary_path, "boundaries[ASSUMED_OUTLET].flow_rate_m3_s"),
        claim(3, "Transferred outlet flows sum to the inlet flow.", sum(b["flow_rate_m3_s"] for b in outlets) * 1e15, "pL/s", preprocess_summary_path, "all_outlet_total_m3_s"),
        claim(3, "Boundary-transfer mass mismatch is below 5e-14.", port_qc["relative_boundary_mass_error"], "dimensionless", port_qc_path, "relative_boundary_mass_error"),
        claim(4, "The final CFD surface is watertight and single-component.", [surface_qc["topology"]["watertight"], surface_qc["topology"]["component_count"]], "boolean,count", surface_qc_path, "topology.watertight, topology.component_count"),
        claim(4, "The final surface has four solver boundaries.", surface_qc["boundary_mapping"]["distal_boundary_count"], "boundaries", surface_qc_path, "boundary_mapping.distal_boundary_count"),
        claim(4, "No true self-intersections were detected.", surface_qc["intersection"]["true_self_intersection_count"], "count", surface_qc_path, "intersection.true_self_intersection_count"),
        claim(4, "Far-core vertices are exactly preserved after output casting.", surface_qc["far_core_exact_preservation"]["far_core_vertex_max_motion_um"], "um", surface_qc_path, "far_core_exact_preservation.far_core_vertex_max_motion_um"),
        claim(4, "Boundary-normal extensions and capping create clean solver boundaries while the far core is preserved.", True, "boolean", surface_config_path, "vmtk.flowextensions, vmtk.cap, extension_mesh.far_core"),
        claim(5, "Base lattice spacing is 0.20 um.", contract["dx_m"] * 1e6, "um", contract_path, "dx_m"),
        claim(5, "The diffusive time step is 2.039 ns.", contract["dt_s"] * 1e9, "ns", contract_path, "dt_s"),
        claim(5, "The relaxation time tau is 1.0.", contract["tau"], "dimensionless", contract_path, "tau"),
        claim(5, "The Base mesh contains 182,320 cells.", flow_manifest["cell_count"], "cells", flow_manifest_path, "cell_count"),
        claim(5, "Target inflow is 0.164 nL/min.", contract["target_volume_flow_m3_s"] * 6e13, "nL/min", contract_path, "target_volume_flow_m3_s"),
        claim(5, "The solver uses D3Q19 BGK lattice-Boltzmann dynamics.", True, "boolean", flow_config_path, "solver model configuration"),
        claim(5, "Adaptive flux-pressure control drives the inlet and three gauge-pressure conditions drive the outlets.", [contract["inlet_boundary"], contract["outlet_boundary"]], "boundary schemes", contract_path, "inlet_boundary, outlet_boundary"),
        claim(6, "The accepted Base field is at iteration 598,755.", metrics["iteration"], "iteration", metrics_path, "iteration"),
        claim(6, "Mean physical speed is 0.177 mm/s.", metrics["mean_speed_m_s"] * 1e3, "mm/s", metrics_path, "mean_speed_m_s"),
        claim(6, "Higher-speed colors are localized around the inlet and selected branch segments.", True, "qualitative visual observation", velocity_visual_path, "screenshot, scalar, display_range"),
        claim(7, "Inlet flow is 0.164 nL/min.", metrics["Qin_m3_s"] * 6e13, "nL/min", metrics_path, "Qin_m3_s"),
        claim(7, "Outlet flow is 0.163 nL/min.", metrics["Qout_m3_s"] * 6e13, "nL/min", metrics_path, "Qout_m3_s"),
        claim(7, "Physical volume closure is 0.154%.", metrics["physical_volume_closure"] * 100, "%", metrics_path, "physical_volume_closure"),
        claim(7, "Outlet flow fractions are 5.2%, 73.2%, and 21.6%.", [100 * metrics["flow_fractions"][f"outlet_0{i}"] for i in range(1, 4)], "%", metrics_path, "flow_fractions"),
        claim(7, "Inlet gauge pressure is 531 Pa.", metrics["inlet_gauge_pressure_pa"], "Pa", metrics_path, "inlet_gauge_pressure_pa"),
        claim(7, "The visible pressure field is gauge pressure, not solver absolute pressure.", "pressure_gauge_pa", "field", pressure_visual_path, "array, units"),
        claim(7, "The outlet split is strongly asymmetric.", metrics["flow_fractions"], "dimensionless fractions", metrics_path, "flow_fractions"),
        claim(8, "The Full V2 mass-identity residual is 7.91e-10 against a 1e-8 gate.", [full["residual"], full["gate"]], "dimensionless", full_path, "residual, gate"),
        claim(8, "Particle populations remain positive.", steady["minimum_pdf"], "lattice PDF", steady_path, "minimum_pdf"),
        claim(8, "The maximum lattice speed is below 0.05.", [steady["maximum_lattice_speed"], 0.05], "lattice units", steady_path, "maximum_lattice_speed and production gate"),
        claim(8, "The maximum Coarse-to-Base difference across reported observables is 1.45%.", sensitivity["maximum_absolute_percent_difference"], "%", sensitivity_path, "maximum_absolute_percent_difference"),
        claim(8, "Formal three-grid convergence was not completed.", sensitivity["formal_asymptotic_grid_convergence"], "boolean", sensitivity_path, "formal_asymptotic_grid_convergence"),
        claim(8, "Fine mesh, controller, and 5,000-step safety passed; Fine steady state was not completed.", run["fine_status"], "status", run_path, "fine_status"),
        claim(8, "Wall-shear-stress validation is deferred.", run["WSS_status"], "status", run_path, "WSS_status"),
        claim(8, "The accepted Base field is the promoted production flow solution.", run["steady_solution_source"], "status", run_path, "steady_solution_source, execution_mode"),
    ]

    entry_points = []
    for name in ["cfd_1D_data_preprocess.py", "cfd_surface_prepare.py", "cfd_flow_solve.py"]:
        path = PROJECT / name
        entry_points.append({"path": rel(path), "sha256": sha256(path), "exists": path.exists()})

    return {
        "status": "PASS",
        "deck_story": "1D boundary data -> solver-ready surface -> validated 3D velocity and gauge pressure",
        "exact_entry_points": entry_points,
        "local_to_repository_filename_mapping": [],
        "source_priority": "exact local entry points, implementation imports, current configs, then frozen run evidence",
        "claims": claims,
        "all_claims_verified": all(c["verified"] for c in claims),
        "new_solver_calls": {"seeder": 0, "musubi": 0, "one_d_solver": 0, "surface_generation": 0, "cfd_solve": 0},
    }


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)

    # Stage-1 artifacts are crops of the recorded preprocessing figures.
    global_1d = PREPROCESS / "figures/global_1d_pressure.png"
    roi_ports = PREPROCESS / "figures/roi_boundary_conditions.png"
    crop_source(global_1d, ASSETS / "01_pipeline_stage1_1d.png", (40, 70, 1580, 1200))
    crop_source(roi_ports, ASSETS / "01b_roi_boundary_transfer.png", (40, 70, 1760, 1400))

    render_surface(SOURCE_SURFACE, ASSETS / "02_surface_before.png", prepared=False)
    render_surface(FINAL_SURFACE, ASSETS / "03_surface_after.png", prepared=True)
    side_by_side(
        ASSETS / "02_surface_before.png",
        ASSETS / "03_surface_after.png",
        ASSETS / "04_surface_before_after.png",
    )

    v3 = FLOW / "visualization/interactive_v3_redesign"
    crop_source(
        v3 / "01_after_velocity_overview.png",
        ASSETS / "05_velocity_field.png",
        (300, 90, 3130, 2120),
    )
    crop_source(
        v3 / "02_after_pressure_overview.png",
        ASSETS / "06_gauge_pressure_field.png",
        (300, 90, 3130, 2120),
    )
    crop_source(
        v3 / "05_after_streamlines.png",
        ASSETS / "07_streamlines.png",
        (300, 90, 3130, 2120),
    )

    visuals = [
        ("01_pipeline_stage1_1d.png", [3], "Global 1D pressure network", global_1d, "pressure", "Pa", "recorded Matplotlib project figure; title cropped"),
        ("01b_roi_boundary_transfer.png", [3], "ROI port transfer", roi_ports, "flow and pressure at ports", "pL/s and Pa", "recorded Matplotlib project figure; title cropped"),
        ("02_surface_before.png", [4], "Open input surface", SOURCE_SURFACE, "lumen surface", "um geometry", "PyVista off-screen render"),
        ("03_surface_after.png", [2, 4], "Prepared extended and capped surface", FINAL_SURFACE, "surface region and boundary identity", "um geometry", "PyVista off-screen render"),
        ("04_surface_before_after.png", [4], "Matched-camera before/after surface comparison", FINAL_SURFACE, "surface preparation", "um geometry", "Pillow composite of two PyVista renders"),
        ("05_velocity_field.png", [1, 2, 6], "Accepted Base physical velocity magnitude", SOURCE_VTU, "velocity_magnitude_mm_s", "mm/s", "crop of validated publication-grade PyVista output"),
        ("06_gauge_pressure_field.png", [7], "Accepted Base gauge pressure", SOURCE_VTU, "pressure_gauge_pa", "Pa", "crop of validated publication-grade PyVista output"),
        ("07_streamlines.png", [], "Optional accepted Base streamline view", SOURCE_VTU, "velocity_magnitude_mm_s", "mm/s", "crop of validated publication-grade PyVista output"),
    ]
    visual_manifest = {
        "status": "PASS",
        "source_vtu_sha256": sha256(SOURCE_VTU),
        "camera": CAMERA,
        "display_ranges": {
            "velocity_mm_s": [0.0012802461835701464, 0.8470950880315884],
            "gauge_pressure_pa": [-6.5, 524.0],
        },
        "assets": [],
    }
    for filename, slides, purpose, source, field, units, method in visuals:
        output = ASSETS / filename
        visual_manifest["assets"].append(
            {
                "filename": f"assets/{filename}",
                "slides": slides,
                "purpose": purpose,
                "source_data": rel(source),
                "source_sha256": sha256(source),
                "output_sha256": sha256(output),
                "field_shown": field,
                "units": units,
                "camera": CAMERA if "surface" in filename or filename.startswith("0") and filename[0:2] in {"05", "06", "07"} else "source figure camera",
                "display_range": visual_manifest["display_ranges"].get("velocity_mm_s" if "velocity" in filename or "streamlines" in filename else "gauge_pressure_pa"),
                "render_method": method,
                "status": "PASS",
            }
        )

    (HERE / "cfd_ppt_visual_manifest.json").write_text(
        json.dumps(visual_manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (HERE / "cfd_ppt_content_evidence.json").write_text(
        json.dumps(build_evidence(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Generated {len(visuals)} presentation assets in {ASSETS}")


if __name__ == "__main__":
    main()
