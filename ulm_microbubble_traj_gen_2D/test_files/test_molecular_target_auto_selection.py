from __future__ import annotations

import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_auto_selection import (
    select_automatic_influence_anchor,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_candidate_io import (
    save_spatially_heterogeneous_target_mask,
    save_spatially_heterogeneous_target_report,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_candidates import (
    MolecularTargetCandidate,
    MolecularTargetCandidateCatalog,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_spatial_heterogeneity import (
    build_spatially_heterogeneous_target,
    evaluate_random_fourier_field,
    generate_random_fourier_coefficients,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.vessel_bed_topology import VesselBedTopology


class MolecularTargetAutoSelectionTests(unittest.TestCase):
    def test_anchor_excludes_single_units_and_inaccessible_subtrees(self) -> None:
        catalog = _catalog(
            _candidate("segment:1", 0.10, visits=100.0, kind="single_vessel_unit"),
            _candidate("subtree:2", 0.11, visits=0.99),
            _candidate("subtree:3", 0.13, visits=1.0),
        )

        result = select_automatic_influence_anchor(catalog, 0.10)

        self.assertEqual(result.anchor_candidate_id, "subtree:3")
        reasons = {
            evaluation.candidate_id: evaluation.rejection_reason
            for evaluation in result.evaluations
        }
        self.assertEqual(reasons["segment:1"], "not_a_downstream_subtree")
        self.assertEqual(reasons["subtree:2"], "expected_bubble_visits_below_one")

    def test_anchor_ties_use_depth_compactness_then_stable_id(self) -> None:
        catalog = _catalog(
            _candidate("subtree:z", 0.10, visits=10.0, depth=2, radius=0.1),
            _candidate("subtree:c", 0.10, visits=10.0, depth=3, radius=2.0),
            _candidate("subtree:b", 0.10, visits=10.0, depth=3, radius=1.0),
            _candidate("subtree:a", 0.10, visits=10.0, depth=3, radius=1.0),
        )

        result = select_automatic_influence_anchor(catalog, 0.10)

        self.assertEqual(result.anchor_candidate_id, "subtree:a")

    def test_spatial_target_is_repeatable_area_weighted_and_wall_limited(self) -> None:
        catalog = _catalog(_candidate("subtree:1", 0.40, visits=10.0))
        anchor = select_automatic_influence_anchor(catalog, 0.40)

        first = _spatial_target(catalog, anchor, seed=42)
        second = _spatial_target(catalog, anchor, seed=42)
        different = _spatial_target(catalog, anchor, seed=43)

        np.testing.assert_array_equal(
            first.target_positive_wall_mask,
            second.target_positive_wall_mask,
        )
        np.testing.assert_allclose(
            first.influence_wall_random_field,
            second.influence_wall_random_field,
        )
        self.assertFalse(
            np.array_equal(
                first.target_positive_wall_mask,
                different.target_positive_wall_mask,
            )
        )
        self.assertFalse(
            np.any(first.target_positive_wall_mask & ~first.influence_wall_mask)
        )
        self.assertTrue(np.any(first.influence_wall_mask[:, 7]))
        self.assertTrue(np.any(first.influence_wall_mask[:, 13]))
        self.assertAlmostEqual(
            first.achieved_positive_wall_fraction_within_influence,
            0.50,
            delta=0.08,
        )
        self.assertGreaterEqual(first.patch_count, 1)

    def test_random_field_is_defined_by_physical_coordinates_not_grid_indices(self) -> None:
        wavevectors, phases = generate_random_fourier_coefficients(5.0, 17, 256)
        shared_coordinates = np.asarray([[1.0, 2.0], [4.0, 8.0], [9.5, 3.5]])

        direct = evaluate_random_fourier_field(
            shared_coordinates,
            np.asarray([0.0, 0.0]),
            wavevectors,
            phases,
        )
        refined_grid_values = evaluate_random_fourier_field(
            np.vstack((shared_coordinates, [[1.5, 2.5], [4.5, 8.5]])),
            np.asarray([0.0, 0.0]),
            wavevectors,
            phases,
        )

        np.testing.assert_allclose(direct, refined_grid_values[:3], rtol=0.0, atol=0.0)

    def test_speed_shear_and_residence_descriptors_do_not_change_spatial_target(self) -> None:
        catalog = _catalog(_candidate("subtree:1", 0.40, visits=10.0))
        changed_candidate = replace(
            catalog.candidates[0],
            inlet_flow_um3_s=999.0,
            residence_time_s=999.0,
            mean_wall_shear_pa=999.0,
        )
        changed_catalog = replace(catalog, candidates=(changed_candidate,))

        first_anchor = select_automatic_influence_anchor(catalog, 0.40)
        changed_anchor = select_automatic_influence_anchor(changed_catalog, 0.40)
        first = _spatial_target(catalog, first_anchor, seed=8)
        changed = _spatial_target(changed_catalog, changed_anchor, seed=8)

        np.testing.assert_array_equal(
            first.target_positive_wall_mask,
            changed.target_positive_wall_mask,
        )

    def test_spatial_npz_and_report_preserve_reproduction_parameters(self) -> None:
        catalog = _catalog(_candidate("subtree:1", 0.40, visits=10.0))
        anchor = select_automatic_influence_anchor(catalog, 0.40)
        result = _spatial_target(catalog, anchor, seed=19)

        with tempfile.TemporaryDirectory() as directory:
            target_path = Path(directory) / "target.npz"
            report_path = Path(directory) / "target.json"
            save_spatially_heterogeneous_target_mask(target_path, catalog, result)
            save_spatially_heterogeneous_target_report(
                report_path,
                catalog,
                anchor,
                result,
            )
            with np.load(target_path, allow_pickle=False) as data:
                self.assertEqual(
                    str(data["selection_mode"].item()),
                    "automatic_spatial_heterogeneity",
                )
                self.assertEqual(int(data["random_seed"].item()), 19)
                np.testing.assert_array_equal(data["target_mask"], result.target_positive_wall_mask)
                np.testing.assert_allclose(
                    data["random_wavevectors_um_inv"],
                    result.random_wavevectors_um_inv,
                )
            report = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["schema_version"], "v2")
        self.assertEqual(report["random_field"]["seed"], 19)
        self.assertIn("wall shear stress", report["default_probability_excludes"])

    def test_legacy_catalog_cannot_drive_v10_automatic_generation(self) -> None:
        catalog = replace(_catalog(_candidate("subtree:1", 0.10, visits=10.0)), schema_version="v2")

        with self.assertRaisesRegex(ValueError, "v3 candidate catalog"):
            select_automatic_influence_anchor(catalog, 0.10)


def _spatial_target(catalog, anchor, *, seed: int):
    return build_spatially_heterogeneous_target(
        catalog,
        anchor,
        influence_wall_area_fraction=0.40,
        positive_wall_fraction_within_influence=0.50,
        correlation_length_um=2.0,
        random_seed=seed,
        random_field_modes=256,
    )


def _candidate(
    candidate_id: str,
    area_fraction: float,
    *,
    visits: float,
    kind: str = "downstream_subtree",
    depth: int = 1,
    radius: float = 1.0,
) -> MolecularTargetCandidate:
    return MolecularTargetCandidate(
        candidate_id=candidate_id,
        kind=kind,
        label=candidate_id,
        root_unit_id=0,
        member_unit_ids=(0,),
        parent_candidate_id=None,
        depth=depth,
        topology_depth=depth,
        inlet_flow_um3_s=1.0,
        root_flow_fraction=0.1,
        network_flow_fraction=0.1,
        volume_um3=1.0,
        residence_time_s=1.0,
        mean_wall_shear_pa=1.0,
        endothelial_wall_area_um2=100.0 * area_fraction,
        endothelial_wall_area_fraction=area_fraction,
        wall_area_centroid_x_um=10.0,
        wall_area_centroid_z_um=10.0,
        radius_of_gyration_um=radius,
        expected_bubble_visits=visits,
    )


def _catalog(*candidates: MolecularTargetCandidate) -> MolecularTargetCandidateCatalog:
    shape = (21, 21)
    x = np.arange(shape[0], dtype=np.float64)
    z = np.arange(shape[1], dtype=np.float64)
    solid_wall = np.zeros(shape, dtype=bool)
    solid_wall[3:18, 7] = True
    solid_wall[3:18, 13] = True
    weights = np.zeros(shape, dtype=np.float64)
    weights[solid_wall] = 1.0
    accessible = solid_wall.copy()
    unit_grid = np.zeros(shape, dtype=np.int32)
    return MolecularTargetCandidateCatalog(
        x_coordinates_um=x,
        z_coordinates_um=z,
        unit_id_grid=unit_grid,
        candidate_support_mask=solid_wall.copy(),
        lumen_mask=np.zeros(shape, dtype=bool),
        solid_wall_mask=solid_wall,
        open_boundary_mask=np.zeros(shape, dtype=bool),
        topology=VesselBedTopology(units=(), segment_id_to_unit_id={}, root_unit_ids=()),
        candidates=tuple(candidates),
        unresolved_junction_cells=0,
        network_endothelial_wall_area_um2=float(np.sum(weights)),
        network_inlet_flow_um3_s=10.0,
        injection_rate_per_s=10.0,
        observation_time_s=1.0,
        wall_area_weight_um2=weights,
        wall_segment_id_grid=np.where(solid_wall, 1, -1).astype(np.int32),
        accessible_wall_mask=accessible,
        expected_bubble_visits_by_unit=np.asarray([10.0]),
        mapped_endothelial_wall_area_um2=float(np.sum(weights)),
        unmapped_endothelial_wall_area_um2=0.0,
    )


if __name__ == "__main__":
    unittest.main()
