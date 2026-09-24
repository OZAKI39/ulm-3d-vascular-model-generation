import numpy as np


def test_first_event_position_and_time_match_original_outlet_triangle(trajectory_cases,classified_segments,p1_audited):
    for (rows,summary),events in zip(trajectory_cases,classified_segments):
        event=summary['event']
        assert summary['exit_boundary'] in ['OUTLET_01','OUTLET_02','OUTLET_03']
        assert events[-1].role==summary['exit_boundary']==rows[-1]['boundary_event']
        start=np.array(event['segment_start_m']); end=np.array(event['unmodified_trial_endpoint_m'])
        point=np.array(event['position_m'])
        np.testing.assert_allclose(point,start+event['segment_fraction']*(end-start),rtol=0,atol=8*np.finfo(float).eps*np.max(np.abs(point)))
        assert summary['exit_time_s']==(event['segment_index']+event['segment_fraction'])*summary['validation_dt_s']
        surface=p1_audited[3][summary['exit_boundary']]
        nodes=surface.faces.reshape(-1,4)[event['role_triangle_id'],1:]
        triangle=np.asarray(surface.points)[nodes]
        # Independent least-squares face coordinates, without the locator/plane code.
        basis=(triangle[1:]-triangle[0]).T
        weights=np.linalg.lstsq(basis,point-triangle[0],rcond=None)[0]
        np.testing.assert_allclose(basis@weights,point-triangle[0],atol=32*np.finfo(float).eps*np.max(np.abs(point)),rtol=0)
        assert min(*weights,1-weights.sum())>=-256*np.finfo(float).eps
