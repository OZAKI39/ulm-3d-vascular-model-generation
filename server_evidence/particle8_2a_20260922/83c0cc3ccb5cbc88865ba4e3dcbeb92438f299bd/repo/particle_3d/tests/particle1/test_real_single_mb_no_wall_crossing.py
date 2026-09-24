import numpy as np


def test_every_real_segment_reclassified_without_wall_or_inlet_contact(trajectory_cases,classified_segments):
    for (rows,summary),events in zip(trajectory_cases,classified_segments):
        assert len(events)==summary['checked_segment_count']==len(rows)-1
        assert all(event is None for event in events[:-1])
        assert events[-1] is not None
        assert events[-1].role.startswith('OUTLET_')
        assert not summary['wall_crossing'] and not summary['inlet_crossing']
        assert all(r['boundary_event']=='ACTIVE' for r in rows[:-1])


def test_finite_records_and_unmodified_explicit_euler_segments(trajectory_cases):
    for rows,summary in trajectory_cases:
        dt=summary['validation_dt_s']
        for row in rows:
            assert all(np.isfinite(float(v)) for k,v in row.items() if k not in ['boundary_event','timestep_role'])
        for a,b in zip(rows[:-2],rows[1:-1]):
            position=np.array([float(a[f'{axis}_m']) for axis in 'xyz'])
            velocity=np.array([float(a[f'V_{axis}_m_s']) for axis in 'xyz'])
            actual=np.array([float(b[f'{axis}_m']) for axis in 'xyz'])
            np.testing.assert_array_equal(actual,position+dt*velocity)
        # No arbitrary displacement acceptance threshold: full trial maximum is reported.
        assert summary['visual_step_jump_review']=='PENDING_USER_REVIEW'
