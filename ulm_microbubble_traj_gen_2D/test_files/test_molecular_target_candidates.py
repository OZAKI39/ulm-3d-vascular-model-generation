from __future__ import annotations

import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_candidate_io import (
    load_candidate_catalog,
    save_candidate_catalog,
    save_selected_target_mask,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_candidates import (
    _resolve_junction_ownership,
    build_molecular_target_candidates,
    simplify_candidate_selection,
)
from ulm_microbubble_traj_gen_2D.utils.molecular import molecular_target_candidates
from ulm_microbubble_traj_gen_2D.utils.particles.particle_hydrodynamic_fields import (
    build_particle_hydrodynamic_fields,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.vessel_bed_topology import (
    build_vessel_bed_topology,
)
from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import straight_channel_case
from ulm_microbubble_traj_gen_2D.utils.visualization.apps.trame_molecular_target_selector import (
    create_target_selector_server,
)
from ulm_microbubble_traj_gen_2D.utils.visualization.vtk.vtk_flow_grid import LUMEN_ARRAY, SPEED_ARRAY


class VesselBedTopologyTests(unittest.TestCase):
    def test_natural_units_stop_at_branches_and_form_downstream_subtrees(self) -> None:
        vessels = _branching_vessels()

        topology = build_vessel_bed_topology(vessels)

        self.assertEqual(topology.units[0].segment_ids, (0, 1))
        self.assertEqual(topology.units[0].child_unit_ids, (1, 2))
        self.assertEqual(topology.units[1].child_unit_ids, (3, 4))
        self.assertEqual(topology.descendants(1), (3, 4))
        self.assertEqual(topology.segment_id_to_unit_id[1], 0)
        self.assertAlmostEqual(topology.units[1].volume_um3, np.pi * 4.0**2)
        self.assertAlmostEqual(
            topology.units[1].endothelial_wall_area_um2,
            2.0 * np.pi * 4.0,
        )
        self.assertEqual(topology.units[1].topology_depth, 1)

    def test_negative_single_segment_flow_reverses_direction_without_changing_input(self) -> None:
        vessel = _vessel(7, -1, [], 0.0, 1.0, -5.0, 2.0)

        topology = build_vessel_bed_topology([vessel])

        self.assertEqual(topology.root_unit_ids, (0,))
        self.assertEqual(topology.units[0].segment_ids, (7,))
        self.assertEqual(topology.units[0].flow_rate_um3_s, 5.0)
        self.assertEqual(vessel.flow_rate, -5.0)


class MolecularTargetCandidateTests(unittest.TestCase):
    def test_junction_cells_follow_the_accepted_flow_into_downstream_branch(self) -> None:
        domain, raster, flow = straight_channel_case()
        junction = np.zeros(domain.shape, dtype=bool)
        junction[6:8] = raster.lumen_mask[6:8]
        vessel_ids = raster.vessel_id.copy()
        vessel_ids[8:][raster.lumen_mask[8:]] = 2
        raster = replace(
            raster,
            vessel_id=vessel_ids,
            junction_core_mask=junction,
        )
        vessels = [
            _vessel(0, -1, [1, 2], 0.0, 7.0, 50.0, 3.0),
            _vessel(1, 0, [], 7.0, 19.0, 25.0, 2.0),
            _vessel(2, 0, [], 7.0, 19.0, 25.0, 2.0),
        ]
        fields = _straight_hydrodynamic_fields(domain, raster, flow)

        catalog = build_molecular_target_candidates(
            domain,
            raster,
            flow,
            fields,
            vessels,
        )

        downstream_unit = catalog.topology.segment_id_to_unit_id[2]
        np.testing.assert_array_equal(
            catalog.unit_id_grid[junction],
            np.full(np.count_nonzero(junction), downstream_unit),
        )

    @unittest.skipIf(
        molecular_target_candidates.njit is None,
        "Numba is not installed in this environment.",
    )
    def test_numba_junction_trace_matches_python_reference_with_stagnant_cells(self) -> None:
        domain, raster, flow = straight_channel_case()
        junction = np.zeros(domain.shape, dtype=bool)
        junction[6:9] = raster.lumen_mask[6:9]
        raster = replace(raster, junction_core_mask=junction)
        ownership = np.zeros(domain.shape, dtype=np.int32)
        ownership[9:] = 1
        velocity = np.asarray(flow.velocity_xz_um_s, dtype=np.float64).copy()
        velocity[:, :2] = 0.0
        flow = replace(flow, velocity_xz_um_s=velocity)

        reference, reference_unresolved = _resolve_junction_ownership(
            domain,
            raster,
            flow,
            raster.lumen_mask,
            ownership,
            use_numba=False,
        )
        accelerated, accelerated_unresolved = _resolve_junction_ownership(
            domain,
            raster,
            flow,
            raster.lumen_mask,
            ownership,
            use_numba=True,
        )

        self.assertGreater(reference_unresolved, 0)
        self.assertEqual(accelerated_unresolved, reference_unresolved)
        np.testing.assert_array_equal(accelerated, reference)

    def test_candidate_mask_is_coordinate_aligned_and_excludes_open_caps(self) -> None:
        domain, raster, flow = straight_channel_case()
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 3.0]),
            x_d=np.asarray([19.0, 0.0, 3.0]),
            radius=3.0,
            flow_rate=50.0,
        )
        fields = build_particle_hydrodynamic_fields(
            domain,
            raster,
            flow,
            continuous_geometry=build_continuous_vessel_geometry(
                [vessel], domain
            ),
        )

        catalog = build_molecular_target_candidates(
            domain,
            raster,
            flow,
            fields,
            [vessel],
        )
        mask = catalog.mask_for_candidate_ids(["segment:0"])

        self.assertEqual(catalog.shape, domain.shape)
        self.assertEqual(mask.dtype, np.bool_)
        self.assertFalse(np.any(mask[fields.open_boundary_mask]))
        self.assertFalse(np.any(mask & ~(catalog.lumen_mask | catalog.solid_wall_mask)))
        self.assertGreater(np.count_nonzero(mask & catalog.solid_wall_mask), 0)
        candidate = catalog.candidate_by_id("segment:0")
        self.assertAlmostEqual(candidate.volume_um3, np.pi * 3.0**2 * 19.0)
        self.assertAlmostEqual(candidate.residence_time_s, candidate.volume_um3 / 50.0)

    def test_parent_selection_removes_redundant_child_and_round_trips_npz(self) -> None:
        catalog = _branching_catalog()
        selected = simplify_candidate_selection(
            catalog,
            ["subtree:2", "segment:4", "subtree:2"],
        )
        self.assertEqual(selected, ("subtree:2",))

        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            catalog_path = directory_path / "candidates.npz"
            json_path = directory_path / "candidates.json"
            target_path = directory_path / "selected_target.npz"
            save_candidate_catalog(catalog_path, json_path, catalog)
            loaded = load_candidate_catalog(catalog_path)
            saved_ids = save_selected_target_mask(
                target_path,
                loaded,
                ["subtree:2", "segment:4"],
            )

            self.assertEqual(saved_ids, ("subtree:2",))
            self.assertTrue(json_path.is_file())
            with np.load(target_path, allow_pickle=False) as data:
                self.assertEqual(data["target_mask"].dtype, np.bool_)
                np.testing.assert_array_equal(data["x_um"], loaded.x_coordinates_um)
                np.testing.assert_array_equal(data["z_um"], loaded.z_coordinates_um)
                self.assertFalse(np.any(data["target_mask"][loaded.open_boundary_mask]))

    def test_candidate_accessibility_uses_network_flow_and_formal_observation_time(self) -> None:
        catalog = _branching_catalog(injection_rate_per_s=20.0, observation_time_s=2.5)
        candidate = catalog.candidate_by_id("subtree:2")

        self.assertAlmostEqual(catalog.network_inlet_flow_um3_s, 100.0)
        self.assertAlmostEqual(candidate.network_flow_fraction, 0.6)
        self.assertAlmostEqual(candidate.expected_bubble_visits, 30.0)
        self.assertAlmostEqual(
            candidate.endothelial_wall_area_fraction,
            candidate.endothelial_wall_area_um2
            / catalog.network_endothelial_wall_area_um2,
        )
        self.assertTrue(catalog.automatic_metrics_available)
        self.assertAlmostEqual(
            catalog.mapped_endothelial_wall_area_um2
            + catalog.unmapped_endothelial_wall_area_um2,
            catalog.network_endothelial_wall_area_um2,
        )
        self.assertAlmostEqual(
            float(np.sum(catalog.wall_area_weight_um2)),
            catalog.mapped_endothelial_wall_area_um2,
        )

    def test_v1_catalog_remains_loadable_for_manual_selection(self) -> None:
        catalog = _branching_catalog()
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            v2_path = directory_path / "v2.npz"
            v1_path = directory_path / "v1.npz"
            save_candidate_catalog(v2_path, directory_path / "v2.json", catalog)
            new_keys = {
                "unit_endothelial_wall_area_um2",
                "unit_wall_area_centroid_x_um",
                "unit_wall_area_centroid_z_um",
                "unit_wall_area_second_moment_um4",
                "unit_topology_depth",
                "candidate_topology_depth",
                "candidate_network_flow_fraction",
                "candidate_endothelial_wall_area_um2",
                "candidate_endothelial_wall_area_fraction",
                "candidate_wall_area_centroid_x_um",
                "candidate_wall_area_centroid_z_um",
                "candidate_radius_of_gyration_um",
                "candidate_expected_bubble_visits",
                "network_endothelial_wall_area_um2",
                "network_inlet_flow_um3_s",
                "injection_rate_per_s",
                "observation_time_s",
                "wall_area_weight_um2",
                "wall_segment_id_grid",
                "accessible_wall_mask",
                "expected_bubble_visits_by_unit",
                "mapped_endothelial_wall_area_um2",
                "unmapped_endothelial_wall_area_um2",
            }
            with np.load(v2_path, allow_pickle=False) as data:
                legacy = {
                    key: np.asarray(data[key])
                    for key in data.files
                    if key not in new_keys
                }
            legacy["schema_version"] = np.asarray("v1")
            np.savez_compressed(v1_path, **legacy)

            loaded = load_candidate_catalog(v1_path)

            self.assertEqual(loaded.schema_version, "v1")
            self.assertFalse(loaded.automatic_metrics_available)
            self.assertGreater(
                np.count_nonzero(loaded.mask_for_candidate_ids(["subtree:2"])),
                0,
            )

    def test_v2_catalog_remains_loadable_for_manual_selection_only(self) -> None:
        catalog = _branching_catalog()
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            current_path = directory_path / "v3.npz"
            v2_path = directory_path / "v2.npz"
            save_candidate_catalog(current_path, directory_path / "v3.json", catalog)
            v3_only_keys = {
                "wall_area_weight_um2",
                "wall_segment_id_grid",
                "accessible_wall_mask",
                "expected_bubble_visits_by_unit",
                "mapped_endothelial_wall_area_um2",
                "unmapped_endothelial_wall_area_um2",
            }
            with np.load(current_path, allow_pickle=False) as data:
                legacy = {
                    key: np.asarray(data[key])
                    for key in data.files
                    if key not in v3_only_keys
                }
            legacy["schema_version"] = np.asarray("v2")
            np.savez_compressed(v2_path, **legacy)

            loaded = load_candidate_catalog(v2_path)

            self.assertEqual(loaded.schema_version, "v2")
            self.assertFalse(loaded.automatic_metrics_available)
            self.assertGreater(
                np.count_nonzero(loaded.mask_for_candidate_ids(["subtree:2"])),
                0,
            )

    def test_trame_selector_builds_from_candidate_and_vti_artifacts(self) -> None:
        import pyvista as pv

        catalog = _branching_catalog(injection_rate_per_s=100.0, observation_time_s=1.0)
        with tempfile.TemporaryDirectory() as directory:
            result_dir = Path(directory)
            save_candidate_catalog(
                result_dir / "molecular_target_candidates.npz",
                result_dir / "molecular_target_candidates.json",
                catalog,
            )
            spacing = float(catalog.x_coordinates_um[1] - catalog.x_coordinates_um[0])
            image = pv.ImageData(
                dimensions=(catalog.shape[0] + 1, catalog.shape[1] + 1, 1),
                origin=(
                    float(catalog.x_coordinates_um[0]) - 0.5 * spacing,
                    float(catalog.z_coordinates_um[0]) - 0.5 * spacing,
                    0.0,
                ),
                spacing=(spacing, spacing, spacing),
            )
            image.cell_data[LUMEN_ARRAY] = catalog.lumen_mask.ravel(order="F").astype(np.uint8)
            image.cell_data[SPEED_ARRAY] = np.ones(catalog.shape, dtype=float).ravel(order="F")
            image.save(result_dir / "final_flow_field.vti")

            server = create_target_selector_server(result_dir)

            self.assertGreater(len(server.state.candidate_tree_items), 0)
            self.assertTrue(hasattr(server.controller, "save_target_mask"))
            self.assertTrue(hasattr(server.controller, "preview_automatic_target"))
            server.state.influence_wall_area_fraction = 0.38
            server.state.positive_wall_fraction = 0.50
            server.state.target_correlation_length_um = 2.0
            server.state.target_random_seed = 42
            server.state.target_random_field_modes = 64
            server.controller.preview_automatic_target()
            self.assertEqual(server.state.selection_workflow, "automatic")
            self.assertEqual(server.state.active_candidate_ids, ["subtree:2"])
            self.assertEqual(server.state.selected_candidate_ids, [])
            server.controller.save_target_mask()
            with np.load(
                result_dir / "selected_molecular_target_mask.npz",
                allow_pickle=False,
            ) as data:
                self.assertEqual(
                    str(data["selection_mode"].item()),
                    "automatic_spatial_heterogeneity",
                )
                self.assertGreater(np.count_nonzero(data["target_mask"]), 0)
            self.assertEqual(
                type(server._target_selector_view).__name__,
                "PyVistaRemoteView",
            )
            screenshot = server._target_selector_plotter.screenshot(return_img=True)
            self.assertGreater(int(np.ptp(screenshot)), 0)
            server._target_selector_plotter.close()


def _branching_catalog(
    *,
    injection_rate_per_s: float = float("nan"),
    observation_time_s: float = float("nan"),
):
    domain, raster, flow = straight_channel_case()
    x_regions = ((0, 3, 0), (4, 7, 1), (8, 10, 2), (11, 13, 4), (14, 16, 5))
    for start, end, vessel_id in x_regions:
        region = np.zeros(domain.shape, dtype=bool)
        region[start : end + 1] = raster.lumen_mask[start : end + 1]
        raster.vessel_id[region] = vessel_id
    final_region = np.zeros(domain.shape, dtype=bool)
    final_region[17:] = raster.lumen_mask[17:]
    raster.vessel_id[final_region] = 3
    fields = _straight_hydrodynamic_fields(domain, raster, flow)
    return build_molecular_target_candidates(
        domain,
        raster,
        flow,
        fields,
        _branching_vessels(),
        injection_rate_per_s=injection_rate_per_s,
        observation_time_s=observation_time_s,
    )


def _branching_vessels() -> list[Vessel]:
    return [
        _vessel(0, -1, [1], 0.0, 1.0, 100.0, 5.0),
        _vessel(1, 0, [2, 3], 1.0, 2.0, 100.0, 5.0),
        _vessel(2, 1, [4, 5], 2.0, 3.0, 60.0, 4.0),
        _vessel(3, 1, [], 2.0, 3.0, 40.0, 3.0),
        _vessel(4, 2, [], 3.0, 4.0, 30.0, 2.0),
        _vessel(5, 2, [], 3.0, 4.0, 30.0, 2.0),
    ]


def _straight_hydrodynamic_fields(domain, raster, flow):
    vessel = Vessel(
        vid=0,
        parent_id=-1,
        children=[],
        x_p=np.asarray([0.0, 0.0, 3.0]),
        x_d=np.asarray([19.0, 0.0, 3.0]),
        radius=3.0,
        flow_rate=50.0,
    )
    geometry = build_continuous_vessel_geometry([vessel], domain)
    return build_particle_hydrodynamic_fields(
        domain, raster, flow, continuous_geometry=geometry
    )


def _vessel(
    vessel_id: int,
    parent_id: int,
    children: list[int],
    start_x: float,
    end_x: float,
    flow_rate: float,
    radius: float,
) -> Vessel:
    return Vessel(
        vid=vessel_id,
        parent_id=parent_id,
        children=children,
        x_p=np.asarray([start_x, 0.0, 0.0]),
        x_d=np.asarray([end_x, 0.0, 0.0]),
        radius=radius,
        flow_rate=flow_rate,
    )


if __name__ == "__main__":
    unittest.main()
