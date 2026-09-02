"""Interactive Trame selector for topology-based molecular-target candidates."""

from __future__ import annotations

import argparse
import asyncio
import json
from pathlib import Path

import numpy as np
import yaml

from ...molecular.molecular_target_auto_selection import (
    AutomaticInfluenceAnchorResult,
    select_automatic_influence_anchor,
)
from ...molecular.molecular_target_candidate_io import (
    load_candidate_catalog,
    save_selected_target_mask,
    save_spatially_heterogeneous_target_mask,
    save_spatially_heterogeneous_target_report,
)
from ...molecular.molecular_target_candidates import simplify_candidate_selection
from ...molecular.molecular_target_spatial_heterogeneity import (
    SpatiallyHeterogeneousTargetResult,
    build_spatially_heterogeneous_target,
)
from ..vtk.pyvista_flow import validate_cfd_flow_dependencies
from ..vtk.vtk_flow_grid import LUMEN_ARRAY, SPEED_ARRAY, VELOCITY_ARRAY


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select one or more topology-based candidate vessel beds and export "
            "the coordinate-aware Boolean NPZ used by molecular_target.mask_npz."
        )
    )
    parser.add_argument(
        "--result-dir",
        type=Path,
        required=True,
        help="Candidate-preparation result directory.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output NPZ path. Defaults to <result-dir>/selected_molecular_target_mask.npz.",
    )
    parser.add_argument("--host", default="127.0.0.1", help="Trame server host.")
    parser.add_argument("--port", type=int, default=8080, help="Trame server port.")
    parser.add_argument(
        "--open-browser",
        action="store_true",
        help="Open the selector in the default browser.",
    )
    parser.add_argument(
        "--exit-after-save",
        action="store_true",
        help="Stop the selector after a target has been saved successfully.",
    )
    return parser.parse_args()


def create_target_selector_server(
    result_dir: Path,
    output_path: Path | None = None,
    exit_after_save: bool = False,
):
    """Create the candidate tree, live highlight, and server-side NPZ export."""

    validate_cfd_flow_dependencies()
    import pyvista as pv
    from pyvista.trame.views import PyVistaRemoteView
    from trame.app import get_server
    from trame.ui.vuetify3 import VAppLayout
    from trame.widgets import html, vuetify3

    result_dir = Path(result_dir).resolve()
    candidate_path = result_dir / "molecular_target_candidates.npz"
    field_path = result_dir / "final_flow_field.vti"
    streamline_path = result_dir / "final_streamlines.vtp"
    if not candidate_path.is_file():
        raise ValueError(f"Missing candidate catalog: {candidate_path}")
    if not field_path.is_file():
        raise ValueError(f"Missing final CFD field: {field_path}")

    catalog = load_candidate_catalog(candidate_path)
    image = pv.read(field_path)
    if not isinstance(image, pv.ImageData):
        raise ValueError("The target selector requires final_flow_field.vti to be ImageData.")
    expected_shape = (int(image.dimensions[0] - 1), int(image.dimensions[1] - 1))
    if expected_shape != catalog.shape:
        raise ValueError(
            f"Candidate grid {catalog.shape} does not match final CFD grid {expected_shape}."
        )

    output_path = (
        result_dir / "selected_molecular_target_mask.npz"
        if output_path is None
        else Path(output_path).resolve()
    )
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 900))
    fluid = image.threshold((0.5, 1.5), scalars=LUMEN_ARRAY, preference="cell")
    fluid = fluid.cell_data_to_point_data(pass_cell_data=True)
    plotter.set_background("white")
    plotter.add_mesh(
        fluid,
        scalars=SPEED_ARRAY,
        preference="point",
        cmap="turbo",
        opacity=0.82,
        lighting=False,
        nan_opacity=0.0,
        show_scalar_bar=True,
        scalar_bar_args={"title": "Speed [um/s]"},
        name="accepted_flow",
    )
    centers = fluid.cell_centers()
    if VELOCITY_ARRAY in centers.point_data and centers.n_points > 0:
        vectors = np.asarray(centers.point_data[VELOCITY_ARRAY], dtype=float)
        magnitudes = np.linalg.norm(vectors, axis=1)
        valid = np.isfinite(magnitudes) & (magnitudes > 0.0)
        directions = np.zeros_like(vectors)
        directions[valid] = vectors[valid] / magnitudes[valid, None]
        centers.point_data["flow_direction"] = directions
        valid_indices = np.flatnonzero(valid)
        if valid_indices.size:
            arrow_count = min(600, valid_indices.size)
            selected_indices = valid_indices[
                np.unique(
                    np.rint(
                        np.linspace(0, valid_indices.size - 1, arrow_count)
                    ).astype(int)
                )
            ]
            arrow_points = centers.extract_points(
                selected_indices,
                adjacent_cells=False,
                include_cells=False,
            )
            arrow_points.points[:, 2] = 0.65 * float(image.spacing[0])
            arrows = arrow_points.glyph(
                orient="flow_direction",
                scale=False,
                factor=10.0 * float(image.spacing[0]),
            )
            plotter.add_mesh(
                arrows,
                color="#111827",
                opacity=0.72,
                lighting=False,
                name="flow_direction_arrows",
            )
    if streamline_path.is_file():
        streamlines = pv.read(streamline_path)
        if isinstance(streamlines, pv.PolyData) and streamlines.n_lines > 0:
            plotter.add_mesh(
                streamlines,
                color="white",
                line_width=1.2,
                opacity=0.75,
                name="flow_streamlines",
            )
    plotter.view_xy()
    plotter.enable_parallel_projection()
    plotter.reset_camera()
    fluid_bounds = fluid.bounds
    fluid_width = float(fluid_bounds[1] - fluid_bounds[0])
    fluid_height = float(fluid_bounds[3] - fluid_bounds[2])
    plotter.camera.parallel_scale = 0.55 * max(
        fluid_height,
        fluid_width / 1.20,
    )
    plotter.render()

    server = get_server(name=f"molecular_target_selector_{result_dir.name}", client_type="vue3")
    state = server.state
    controller = server.controller
    state.candidate_tree_items = _candidate_tree_items(catalog)
    automatic_defaults = _automatic_defaults(result_dir)
    state.selected_candidate_ids = []
    state.active_candidate_ids = []
    state.selection_workflow = "manual"
    state.influence_wall_area_fraction = automatic_defaults[0]
    state.positive_wall_fraction = automatic_defaults[1]
    state.target_correlation_length_um = automatic_defaults[2]
    state.target_random_seed = automatic_defaults[3]
    state.target_random_field_modes = automatic_defaults[4]
    state.automatic_metrics_available = catalog.automatic_metrics_available
    state.automatic_status = (
        "Enter influence size, positive-wall fraction, correlation length, and seed, then "
        "preview the spatially heterogeneous synthetic target."
        if catalog.automatic_metrics_available
        else "Automatic selection is unavailable for this legacy catalog. Rebuild it first."
    )
    state.selection_summary = "No candidate selected."
    state.candidate_detail = "Select a tree item to inspect its steady-flow descriptors."
    state.save_status = f"The final Boolean target mask will be saved to:\n{output_path}"
    automatic_anchor: AutomaticInfluenceAnchorResult | None = None
    automatic_result: SpatiallyHeterogeneousTargetResult | None = None

    def update_highlight(selected_candidate_ids=None, **_kwargs) -> None:
        selected = [str(value) for value in (selected_candidate_ids or [])]
        simplified = simplify_candidate_selection(catalog, selected)
        mask = catalog.mask_for_candidate_ids(simplified)
        image.cell_data["selected_target_candidate"] = mask.ravel(order="F").astype(np.uint8)
        plotter.remove_actor("selected_candidate_bed", reset_camera=False)
        plotter.remove_actor("selected_candidate_wall", reset_camera=False)
        plotter.remove_actor("automatic_influence_region", reset_camera=False)
        plotter.remove_actor("automatic_influence_wall", reset_camera=False)
        plotter.remove_actor("automatic_target_wall", reset_camera=False)
        if np.any(mask):
            selected_cells = image.threshold(
                (0.5, 1.5),
                scalars="selected_target_candidate",
                preference="cell",
            )
            selected_cells = selected_cells.copy(deep=True)
            selected_cells.points[:, 2] = 0.35 * float(image.spacing[0])
            plotter.add_mesh(
                selected_cells,
                color="#ff1744",
                opacity=0.32,
                lighting=False,
                show_edges=False,
                name="selected_candidate_bed",
            )
            selected_wall = mask & catalog.solid_wall_mask
            wall_indices = np.argwhere(selected_wall)
            if wall_indices.size:
                points = np.column_stack(
                    [
                        catalog.x_coordinates_um[wall_indices[:, 0]],
                        catalog.z_coordinates_um[wall_indices[:, 1]],
                        np.full(
                            wall_indices.shape[0],
                            0.75 * float(image.spacing[0]),
                        ),
                    ]
                )
                plotter.add_points(
                    points,
                    color="#ffea00",
                    point_size=7.0,
                    render_points_as_spheres=True,
                    name="selected_candidate_wall",
                )
        state.selection_summary = _selection_summary(catalog, simplified, mask)
        if hasattr(controller, "target_view_update"):
            controller.target_view_update()

    def update_automatic_highlight(
        result: SpatiallyHeterogeneousTargetResult,
    ) -> None:
        plotter.remove_actor("selected_candidate_bed", reset_camera=False)
        plotter.remove_actor("selected_candidate_wall", reset_camera=False)
        plotter.remove_actor("automatic_influence_region", reset_camera=False)
        plotter.remove_actor("automatic_influence_wall", reset_camera=False)
        plotter.remove_actor("automatic_target_wall", reset_camera=False)
        image.cell_data["automatic_influence_region"] = (
            result.influence_region_mask.ravel(order="F").astype(np.uint8)
        )
        influence_cells = image.threshold(
            (0.5, 1.5),
            scalars="automatic_influence_region",
            preference="cell",
        ).copy(deep=True)
        influence_cells.points[:, 2] = 0.30 * float(image.spacing[0])
        plotter.add_mesh(
            influence_cells,
            color="#7c3aed",
            opacity=0.10,
            lighting=False,
            show_edges=False,
            name="automatic_influence_region",
        )
        for mask, color, size, actor_name, height in (
            (
                result.influence_wall_mask,
                "#facc15",
                5.0,
                "automatic_influence_wall",
                0.68,
            ),
            (
                result.target_positive_wall_mask,
                "#ff1744",
                9.0,
                "automatic_target_wall",
                0.82,
            ),
        ):
            indices = np.argwhere(mask)
            if not indices.size:
                continue
            points = np.column_stack(
                (
                    catalog.x_coordinates_um[indices[:, 0]],
                    catalog.z_coordinates_um[indices[:, 1]],
                    np.full(indices.shape[0], height * float(image.spacing[0])),
                )
            )
            plotter.add_points(
                points,
                color=color,
                point_size=size,
                render_points_as_spheres=True,
                name=actor_name,
            )
        state.selection_summary = _automatic_selection_summary(result)
        if hasattr(controller, "target_view_update"):
            controller.target_view_update()

    @state.change("selected_candidate_ids")
    def on_selection_change(selected_candidate_ids, **_kwargs) -> None:
        nonlocal automatic_anchor
        nonlocal automatic_result
        if state.selection_workflow == "automatic" and not selected_candidate_ids:
            return
        if selected_candidate_ids:
            state.selection_workflow = "manual"
            automatic_anchor = None
            automatic_result = None
            state.automatic_status = (
                "Manual candidate mode is active. Preview automatic selection again to switch "
                "back to spatial heterogeneity."
            )
        update_highlight(selected_candidate_ids)

    @state.change("active_candidate_ids")
    def on_active_candidate(active_candidate_ids, **_kwargs) -> None:
        active = [str(value) for value in (active_candidate_ids or [])]
        state.candidate_detail = (
            _candidate_detail(catalog, active[-1])
            if active
            else "Select a tree item to inspect its steady-flow descriptors."
        )

    def save_target() -> None:
        if state.selection_workflow == "automatic" and automatic_result is not None:
            if automatic_anchor is None:
                state.save_status = "Automatic preview has no influence anchor."
                return
            save_spatially_heterogeneous_target_mask(
                output_path,
                catalog,
                automatic_result,
            )
            report_path = output_path.with_name(
                "automatic_molecular_target_selection.json"
            )
            save_spatially_heterogeneous_target_report(
                report_path,
                catalog,
                automatic_anchor,
                automatic_result,
            )
            state.save_status = (
                f"Saved the spatially heterogeneous target to:\n{output_path}\n"
                f"Audit report:\n{report_path}\n\n"
                "Use this file with:\n"
                "molecular_target.enabled: true\n"
                "molecular_target.region_mode: mask_npz\n"
                f"molecular_target.mask_npz_path: {output_path}"
            )
            if exit_after_save:
                asyncio.create_task(server.stop())
            return
        selected = [str(value) for value in state.selected_candidate_ids]
        try:
            simplified = save_selected_target_mask(
                output_path,
                catalog,
                selected,
                selection_mode="manual",
                achieved_wall_area_fraction=_selected_wall_area_fraction(
                    catalog,
                    simplify_candidate_selection(catalog, selected),
                ),
            )
        except ValueError as exc:
            state.save_status = str(exc)
            return
        state.save_status = (
            f"Saved {len(simplified)} non-redundant candidate selection(s) to:\n"
            f"{output_path}\n\n"
            "Use this file with:\n"
            "molecular_target.enabled: true\n"
            "molecular_target.region_mode: mask_npz\n"
            f"molecular_target.mask_npz_path: {output_path}"
        )
        if exit_after_save:
            asyncio.create_task(server.stop())

    def preview_automatic_target() -> None:
        nonlocal automatic_anchor
        nonlocal automatic_result
        try:
            influence_fraction = float(state.influence_wall_area_fraction)
            positive_fraction = float(state.positive_wall_fraction)
            correlation_length = float(state.target_correlation_length_um)
            seed_value = float(state.target_random_seed)
            modes_value = float(state.target_random_field_modes)
            if not seed_value.is_integer() or not modes_value.is_integer():
                raise ValueError("Random seed and random-field modes must be whole numbers.")
            anchor = select_automatic_influence_anchor(
                catalog,
                influence_fraction,
            )
            result = build_spatially_heterogeneous_target(
                catalog,
                anchor,
                influence_wall_area_fraction=influence_fraction,
                positive_wall_fraction_within_influence=positive_fraction,
                correlation_length_um=correlation_length,
                random_seed=int(seed_value),
                random_field_modes=int(modes_value),
            )
        except (TypeError, ValueError) as exc:
            state.automatic_status = str(exc)
            return
        automatic_anchor = anchor
        automatic_result = result
        state.selection_workflow = "automatic"
        state.selected_candidate_ids = []
        state.active_candidate_ids = [anchor.anchor_candidate_id]
        update_automatic_highlight(result)
        state.automatic_status = (
            "Spatially heterogeneous synthetic target previewed.\n"
            f"Coarse anchor: {anchor.anchor_candidate_id}\n"
            f"Influence radius: {result.influence_radius_um:.6g} um\n"
            "Positive wall fraction: "
            f"{result.achieved_positive_wall_fraction_within_influence:.6g}\n"
            f"Grid-connected patches: {result.patch_count}\n"
            f"Random seed: {result.random_seed}"
        )

    controller.save_target_mask = save_target
    controller.preview_automatic_target = preview_automatic_target
    controller.clear_target_selection = lambda: state.update(
        {"selected_candidate_ids": [], "selection_workflow": "manual"}
    )
    if (
        (result_dir / "automatic_molecular_target_selection.json").is_file()
        and catalog.automatic_metrics_available
        and all(value is not None for value in automatic_defaults[:3])
    ):
        preview_automatic_target()

    with VAppLayout(server):
        with vuetify3.VMain(style="padding:0;"):
            with html.Div(style="display:flex;width:100%;height:100vh;overflow:hidden;"):
                with html.Div(
                    style=(
                        "width:410px;min-width:410px;height:100%;overflow:auto;"
                        "padding:12px;background:#f7f8fa;border-right:1px solid #d9dde3;"
                    )
                ):
                    html.H2("Candidate vessel beds", style="margin:0 0 6px 0;")
                    html.P(
                        "These are selectable vascular regions, not automatically predicted tumour regions.",
                        style="font-size:13px;color:#5f6368;margin:0 0 10px 0;",
                    )
                    with vuetify3.VCard(variant="outlined", style="margin-bottom:10px;"):
                        vuetify3.VCardTitle("Automatic heterogeneous target")
                        with vuetify3.VCardText():
                            vuetify3.VTextField(
                                label="Influence-region network wall fraction",
                                model_value=("influence_wall_area_fraction",),
                                update_modelValue="influence_wall_area_fraction = $event",
                                type="number",
                                step="0.001",
                                density="compact",
                            )
                            vuetify3.VTextField(
                                label="Positive wall fraction inside influence region",
                                model_value=("positive_wall_fraction",),
                                update_modelValue="positive_wall_fraction = $event",
                                type="number",
                                step="0.01",
                                density="compact",
                            )
                            vuetify3.VTextField(
                                label="Patch correlation length [um]",
                                model_value=("target_correlation_length_um",),
                                update_modelValue="target_correlation_length_um = $event",
                                type="number",
                                step="1",
                                density="compact",
                            )
                            vuetify3.VTextField(
                                label="Random realization seed",
                                model_value=("target_random_seed",),
                                update_modelValue="target_random_seed = $event",
                                type="number",
                                step="1",
                                density="compact",
                            )
                            vuetify3.VTextField(
                                label="Random-field modes (numerical)",
                                model_value=("target_random_field_modes",),
                                update_modelValue="target_random_field_modes = $event",
                                type="number",
                                step="64",
                                density="compact",
                            )
                            vuetify3.VBtn(
                                "Preview automatic selection",
                                color="secondary",
                                click=controller.preview_automatic_target,
                                disabled=(
                                    "!automatic_metrics_available || "
                                    "influence_wall_area_fraction === null || "
                                    "positive_wall_fraction === null || "
                                    "target_correlation_length_um === null",
                                ),
                            )
                            html.Pre(
                                "{{ automatic_status }}",
                                style=(
                                    "white-space:pre-wrap;margin-top:8px;"
                                    "font:11px/1.4 Consolas,monospace;"
                                ),
                            )
                    vuetify3.VTreeview(
                        items=("candidate_tree_items",),
                        item_title="title",
                        item_value="value",
                        item_children="children",
                        selectable=True,
                        select_strategy="independent",
                        activatable=True,
                        active_strategy="single-independent",
                        model_value=("selected_candidate_ids",),
                        activated=("active_candidate_ids",),
                        update_modelValue="selected_candidate_ids = $event",
                        update_activated="active_candidate_ids = $event",
                        density="compact",
                        open_all=True,
                    )
                    with vuetify3.VCard(variant="outlined", style="margin-top:10px;"):
                        vuetify3.VCardTitle("Candidate details")
                        with vuetify3.VCardText():
                            html.Pre(
                                "{{ candidate_detail }}",
                                style="white-space:pre-wrap;font:12px/1.4 Consolas,monospace;",
                            )
                    with vuetify3.VCard(variant="outlined", style="margin-top:10px;"):
                        vuetify3.VCardTitle("Current target")
                        with vuetify3.VCardText():
                            html.Pre(
                                "{{ selection_summary }}",
                                style="white-space:pre-wrap;font:12px/1.4 Consolas,monospace;",
                            )
                    with html.Div(style="display:flex;gap:8px;margin-top:10px;"):
                        vuetify3.VBtn(
                            "Save target NPZ",
                            color="primary",
                            click=controller.save_target_mask,
                        )
                        vuetify3.VBtn(
                            "Clear",
                            variant="outlined",
                            click=controller.clear_target_selection,
                        )
                    html.Pre(
                        "{{ save_status }}",
                        style=(
                            "white-space:pre-wrap;margin-top:10px;padding:8px;"
                            "background:#ffffff;border:1px solid #d9dde3;font:11px/1.4 Consolas,monospace;"
                        ),
                    )
                with html.Div(style="flex:1;height:100%;position:relative;"):
                    view = PyVistaRemoteView(
                        plotter,
                        interactive_ratio=1.0,
                        still_ratio=1.0,
                        style="width:100%;height:100%;",
                    )
                    controller.target_view_update = view.update
                    with vuetify3.VCard(
                        elevation=5,
                        style=(
                            "position:absolute;left:14px;top:14px;z-index:10;"
                            "background:rgba(255,255,255,.90);"
                        ),
                    ):
                        html.Div(
                            "Accepted steady CFD field with molecular-target preview",
                            style="padding:8px 12px;font:600 14px Arial,sans-serif;",
                        )

    server.controller.on_server_exited.add(lambda **_: plotter.close())
    server._target_selector_plotter = plotter
    server._target_selector_catalog = catalog
    server._target_selector_image = image
    server._target_selector_view = view
    return server


def _candidate_tree_items(catalog) -> list[dict[str, object]]:
    items = {
        candidate.candidate_id: {
            "title": _candidate_tree_title(candidate),
            "value": candidate.candidate_id,
            "children": [],
        }
        for candidate in catalog.candidates
    }
    roots: list[dict[str, object]] = []
    for candidate in catalog.candidates:
        item = items[candidate.candidate_id]
        parent_id = candidate.parent_candidate_id
        if parent_id is not None and parent_id in items:
            items[parent_id]["children"].append(item)
        else:
            roots.append(item)
    return roots


def _candidate_tree_title(candidate) -> str:
    kind = "Subtree" if candidate.kind == "downstream_subtree" else "Unit"
    return (
        f"{kind} {candidate.candidate_id.split(':', 1)[1]} | "
        f"flow {candidate.root_flow_fraction:.3g} | "
        f"T {candidate.residence_time_s:.3g} s"
    )


def _candidate_detail(catalog, candidate_id: str) -> str:
    candidate = catalog.candidate_by_id(candidate_id)
    shear = (
        f"{candidate.mean_wall_shear_pa:.6g} Pa"
        if np.isfinite(candidate.mean_wall_shear_pa)
        else "not available"
    )
    return (
        f"ID: {candidate.candidate_id}\n"
        f"Type: {candidate.kind}\n"
        f"Vessel-bed units: {candidate.member_unit_ids}\n"
        f"Inlet flow: {candidate.inlet_flow_um3_s:.6g} um^3/s\n"
        f"Root-flow fraction: {candidate.root_flow_fraction:.6g}\n"
        f"Network-flow fraction: {candidate.network_flow_fraction:.6g}\n"
        f"Expected bubble visits: {candidate.expected_bubble_visits:.6g}\n"
        f"Endothelial wall area: {candidate.endothelial_wall_area_um2:.6g} um^2\n"
        f"Network wall-area fraction: {candidate.endothelial_wall_area_fraction:.6g}\n"
        f"Topology depth: {candidate.topology_depth}\n"
        f"Wall-area radius of gyration: {candidate.radius_of_gyration_um:.6g} um\n"
        f"Volume: {candidate.volume_um3:.6g} um^3\n"
        f"Residence-time scale: {candidate.residence_time_s:.6g} s\n"
        f"Wall-area-weighted mean wall shear: {shear}"
    )


def _selection_summary(catalog, selected: tuple[str, ...], mask: np.ndarray) -> str:
    if not selected:
        return "No candidate selected."
    wall_sites = int(np.count_nonzero(mask & catalog.solid_wall_mask))
    return (
        f"Non-redundant candidates: {len(selected)}\n"
        f"Candidate IDs: {', '.join(selected)}\n"
        f"Selected grid cells: {int(np.count_nonzero(mask))}\n"
        f"Selected eligible wall sites: {wall_sites}"
    )


def _automatic_selection_summary(
    result: SpatiallyHeterogeneousTargetResult,
) -> str:
    return (
        "Mode: automatic spatial heterogeneity\n"
        f"Coarse anchor: {result.anchor_candidate_id}\n"
        f"Influence radius: {result.influence_radius_um:.6g} um\n"
        "Influence network wall fraction: "
        f"{result.achieved_influence_wall_area_fraction:.6g}\n"
        "Positive fraction inside influence: "
        f"{result.achieved_positive_wall_fraction_within_influence:.6g}\n"
        f"Positive wall sites: {int(np.count_nonzero(result.target_positive_wall_mask))}\n"
        f"Grid-connected patches: {result.patch_count}"
    )


def _selected_wall_area_fraction(catalog, selected: tuple[str, ...]) -> float:
    if not selected or not np.isfinite(catalog.network_endothelial_wall_area_um2):
        return float("nan")
    unit_ids: set[int] = set()
    for candidate_id in selected:
        unit_ids.update(catalog.candidate_by_id(candidate_id).member_unit_ids)
    selected_area = sum(
        catalog.topology.units[unit_id].endothelial_wall_area_um2
        for unit_id in unit_ids
    )
    return float(selected_area / catalog.network_endothelial_wall_area_um2)


def _automatic_defaults(
    result_dir: Path,
) -> tuple[float | None, float | None, float | None, int, int]:
    report_path = result_dir / "automatic_molecular_target_selection.json"
    if report_path.is_file():
        document = json.loads(report_path.read_text(encoding="utf-8"))
        if str(document.get("schema_version")) == "v2":
            influence = document["influence_region"]
            patches = document["target_positive_patches"]
            random_field = document["random_field"]
            return (
                float(influence["requested_network_wall_area_fraction"]),
                float(patches["requested_positive_fraction_within_influence"]),
                float(random_field["correlation_length_um"]),
                int(random_field["seed"]),
                int(random_field["mode_count"]),
            )
    run_config_path = result_dir / "run_config.yaml"
    if run_config_path.is_file():
        document = yaml.safe_load(run_config_path.read_text(encoding="utf-8-sig")) or {}
        section = document.get("molecular_target_selection", {})
        influence_fraction = section.get(
            "influence_region_endothelial_wall_area_fraction"
        )
        positive_fraction = section.get(
            "target_positive_wall_fraction_within_influence"
        )
        correlation_length = section.get("target_correlation_length_um")
        return (
            None if influence_fraction is None else float(influence_fraction),
            None if positive_fraction is None else float(positive_fraction),
            None if correlation_length is None else float(correlation_length),
            int(section.get("random_seed", 42)),
            int(section.get("random_field_modes", 512)),
        )
    return None, None, None, 42, 512


def main() -> None:
    args = parse_args()
    server = create_target_selector_server(
        args.result_dir,
        args.output,
        exit_after_save=args.exit_after_save,
    )
    server.start(
        host=args.host,
        port=int(args.port),
        open_browser=bool(args.open_browser),
    )


if __name__ == "__main__":
    try:
        main()
    except (ImportError, ValueError) as exc:
        raise SystemExit(f"Molecular-target selector stopped: {exc}") from None
