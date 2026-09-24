import numpy as np
from particle_3d.validation_boundary import ValidationBoundaryClassifier


def test_all_centers_finite_inside_and_outlet_classified(real_cases,center_samples,p2_audited):
    classifier=ValidationBoundaryClassifier(p2_audited[3])
    for di,samples in enumerate(center_samples):
        assert samples.inside_lumen.all()
        rows,summary=real_cases[di]
        events=[classifier.first_event([a[f"{s}_m"] for s in "xyz"],[b[f"{s}_m"] for s in "xyz"]) for a,b in zip(rows[:-1],rows[1:])]
        assert all(e is None for e in events[:-1])
        assert events[-1].role==summary["exit_boundary"]=="OUTLET_02"
        assert summary["checked_segment_count"]==len(events)
    for rows,summary in real_cases:
        xyz=np.array([[r[f"{a}_m"] for a in "xyz"] for r in rows])
        assert np.isfinite(xyz).all() and all(r["inside_lumen"] for r in rows)
        np.testing.assert_array_equal(xyz,[[r[f"{a}_m"] for a in "xyz"] for r in real_cases[summary["dt_index"]][0]])
        assert not summary["wall_crossing"] and not summary["inlet_crossing"]
        assert summary["finite_size_wall_clearance"]=="NOT_VALIDATED_PARTICLE3"
