import json
import numpy as np
from particle_3d.rbc import RBCGeometry
from particle_3d.rbc_distribution import stratified_indices
from particle_3d.particle2_cases import simple_shear_case,SHEAR_ANGLE_STEPS


def test_all_64_distribution_shapes_with_same_orientation_solver(population,p2_data):
    indices=stratified_indices(population.samples,2026092064)
    geometries=[RBCGeometry.from_population(population,i) for i in indices]
    _,metrics=simple_shear_case(geometries,SHEAR_ANGLE_STEPS[-1])
    saved=json.loads((p2_data/"07_distribution_wide_metrics.json").read_text())
    assert len(metrics)==64 and len(set(m["rbc_id"] for m in metrics))==64
    for m,old,g in zip(metrics,saved,geometries):
        assert m["passed"] and m["rbc_id"]==g.provenance.rbc_id
        assert m["measured_period_s"]==old["measured_period_s"]
        assert m["r"]==g.r and m["D_um"]==g.provenance.D_um
        assert m["relative_error"]<=m["period_relative_error_bound"]
