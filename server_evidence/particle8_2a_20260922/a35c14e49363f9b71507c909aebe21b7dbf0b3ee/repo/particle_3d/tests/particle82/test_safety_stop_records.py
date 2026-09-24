import json,math

def test_all_720_stops_have_explicit_unchanged_evidence(evidence_root):
    s=json.loads((evidence_root/'admission/stop_audit/SAFETY_STOP_AUDIT.json').read_text())
    assert s['count']==len(s['rows'])==len({r['stable_id'] for r in s['rows']})==720
    assert sum(s['by_basin'].values())==720
    for r in s['rows']:
        for k in ['x_m','y_m','z_m','physical_time_s','elapsed_time_s','nearest_wall_gap_m','g_nf_m','h_lower_m','speed_m_s']:
            assert math.isfinite(r[k])
        assert r['subdivision_depth'] is None
        assert r['maximum_recorded_subdivision_depth']>=0 and r['provider_calls']>0 and r['tetra_id']>=0
        assert r['point_tracer_basin'] and 'NOT_TERMINAL_TRIAL_DEPTH' in r['subdivision_depth_availability']
