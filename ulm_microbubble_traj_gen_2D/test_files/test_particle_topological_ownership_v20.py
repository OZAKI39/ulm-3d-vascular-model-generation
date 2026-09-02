from __future__ import annotations

from types import SimpleNamespace
import unittest

import numpy as np
from shapely.geometry import LineString

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    _continuous_lumen_polygon,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_topological_ownership import (
    build_topological_commitment_catalog,
    inspect_topological_crossings,
)


class ParticleTopologicalOwnershipV20Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.vessels = _branching_tree()
        polygon = _continuous_lumen_polygon(self.vessels, quad_segs=32)
        self.catalog = build_topological_commitment_catalog(
            self.vessels,
            SimpleNamespace(lumen_polygon=polygon),
        )

    def test_catalog_has_one_grid_independent_section_per_directed_edge(self) -> None:
        self.assertEqual(self.catalog.section_count, 2)
        np.testing.assert_array_equal(self.catalog.root_vessel_id, [1])
        np.testing.assert_array_equal(self.catalog.parent_vessel_id, [1, 1])
        np.testing.assert_array_equal(self.catalog.child_vessel_id, [2, 3])
        self.assertTrue(
            np.all(
                self.catalog.commitment_distance_um
                > self.catalog.transition_end_distance_um
            )
        )
        lines = [
            LineString(
                (
                    self.catalog.point_xz_um[index]
                    - self.catalog.half_width_um[index]
                    * self.catalog.tangent_xz[index],
                    self.catalog.point_xz_um[index]
                    + self.catalog.half_width_um[index]
                    * self.catalog.tangent_xz[index],
                )
            )
            for index in range(self.catalog.section_count)
        ]
        self.assertFalse(lines[0].intersects(lines[1]))

    def test_forward_and_reverse_crossings_change_only_adjacent_owners(self) -> None:
        section = 0
        point = self.catalog.point_xz_um[section]
        normal = self.catalog.downstream_normal_xz[section]
        forward = inspect_topological_crossings(
            (point - normal)[None, :],
            (point + normal)[None, :],
            np.asarray([1], dtype=np.int32),
            np.asarray([True]),
            self.catalog,
            use_numba=False,
        )
        reverse = inspect_topological_crossings(
            (point + normal)[None, :],
            (point - normal)[None, :],
            np.asarray([2], dtype=np.int32),
            np.asarray([True]),
            self.catalog,
            use_numba=True,
        )
        self.assertAlmostEqual(float(forward.fraction[0]), 0.5)
        self.assertEqual(int(forward.new_vessel_id[0]), 2)
        self.assertAlmostEqual(float(reverse.fraction[0]), 0.5)
        self.assertEqual(int(reverse.new_vessel_id[0]), 1)

    def test_sibling_transfer_must_pass_through_parent_state(self) -> None:
        first_point = self.catalog.point_xz_um[0]
        first_normal = self.catalog.downstream_normal_xz[0]
        second_point = self.catalog.point_xz_um[1]
        second_normal = self.catalog.downstream_normal_xz[1]
        child_to_junction = inspect_topological_crossings(
            (first_point + first_normal)[None, :],
            (first_point - first_normal)[None, :],
            np.asarray([2], dtype=np.int32),
            np.asarray([True]),
            self.catalog,
            use_numba=False,
        )
        parent_to_sibling = inspect_topological_crossings(
            (second_point - second_normal)[None, :],
            (second_point + second_normal)[None, :],
            child_to_junction.new_vessel_id,
            np.asarray([True]),
            self.catalog,
            use_numba=False,
        )
        self.assertEqual(int(child_to_junction.new_vessel_id[0]), 1)
        self.assertEqual(int(parent_to_sibling.new_vessel_id[0]), 3)

    def test_direction_resolves_a_start_exactly_on_the_section(self) -> None:
        point = self.catalog.point_xz_um[0]
        normal = self.catalog.downstream_normal_xz[0]
        forward = inspect_topological_crossings(
            point[None, :],
            (point + normal)[None, :],
            np.asarray([1], dtype=np.int32),
            np.asarray([True]),
            self.catalog,
            use_numba=False,
        )
        tangent = inspect_topological_crossings(
            point[None, :],
            (point + self.catalog.tangent_xz[0])[None, :],
            np.asarray([1], dtype=np.int32),
            np.asarray([True]),
            self.catalog,
            use_numba=False,
        )
        self.assertEqual(float(forward.fraction[0]), 0.0)
        self.assertEqual(int(forward.new_vessel_id[0]), 2)
        self.assertTrue(np.isnan(tangent.fraction[0]))
        self.assertEqual(int(tangent.new_vessel_id[0]), 1)

    def test_numba_and_python_crossing_kernels_match(self) -> None:
        starts = self.catalog.point_xz_um - self.catalog.downstream_normal_xz
        ends = self.catalog.point_xz_um + self.catalog.downstream_normal_xz
        owners = self.catalog.parent_vessel_id.copy()
        live = np.ones(owners.size, dtype=bool)
        reference = inspect_topological_crossings(
            starts, ends, owners, live, self.catalog, use_numba=False
        )
        accelerated = inspect_topological_crossings(
            starts, ends, owners, live, self.catalog, use_numba=True
        )
        np.testing.assert_array_equal(accelerated.fraction, reference.fraction)
        np.testing.assert_array_equal(
            accelerated.new_vessel_id, reference.new_vessel_id
        )
        np.testing.assert_array_equal(
            accelerated.section_index, reference.section_index
        )


def _branching_tree() -> list[Vessel]:
    junction = np.asarray([20.0, 0.0, 0.0])
    return [
        Vessel(
            vid=0,
            parent_id=-1,
            children=[1, 2],
            x_p=np.asarray([0.0, 0.0, 0.0]),
            x_d=junction,
            radius=4.0,
            flow_rate=2.0,
        ),
        Vessel(
            vid=1,
            parent_id=0,
            children=[],
            x_p=junction,
            x_d=np.asarray([45.0, 0.0, 12.0]),
            radius=2.0,
            flow_rate=1.0,
        ),
        Vessel(
            vid=2,
            parent_id=0,
            children=[],
            x_p=junction,
            x_d=np.asarray([45.0, 0.0, -12.0]),
            radius=2.0,
            flow_rate=1.0,
        ),
    ]


if __name__ == "__main__":
    unittest.main()
