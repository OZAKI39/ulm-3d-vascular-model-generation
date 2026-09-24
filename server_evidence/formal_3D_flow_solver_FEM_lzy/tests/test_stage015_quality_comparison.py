from stage015_helpers import REPORT,read


def test_no_uncomputed_tetra_quality_is_reported_as_improvement():
    result=read(REPORT/'quality_comparison.json')
    assert result['baseline']['total_lt_0_1']==153 and result['baseline']['cap_adjacent_lt_0_1']==129
    for candidate in result['candidates'].values():
        if candidate['geometry_status']=='FAIL':
            assert candidate['volume_qc'] is None
            assert candidate['volume_status']=='NOT RUN: geometry gate failed'
