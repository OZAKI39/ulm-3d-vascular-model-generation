from copy import deepcopy
from fem3d.cap_selection import select_candidate
from stage015_helpers import REPORT,read


def test_fewer_tetrahedra_never_overrides_worse_cap_tail():
    policy=read(REPORT/'acceptance_policy.json')
    q={'cap_adjacent_lt_0_1':10,'total_lt_0_1':30,'P1':.4,'P5':.5,'median':.75,'tetrahedra':160000,'invalid_tetra':0}
    good={'candidate':'candidate_A','geometry_status':'PASS','volume_qc':q}
    bad=deepcopy(good);bad['candidate']='candidate_B';bad['volume_qc'].update(cap_adjacent_lt_0_1=11,tetrahedra=100000)
    assert select_candidate([bad,good],policy)['selected_candidate']=='candidate_A'


def test_geometry_failure_disqualifies_numerically_better_candidate():
    policy=read(REPORT/'acceptance_policy.json')
    row={'candidate':'candidate_A','geometry_status':'FAIL','volume_qc':{'cap_adjacent_lt_0_1':0}}
    assert select_candidate([row],policy)['selected_candidate'] is None


def test_recorded_selection_obeys_eligibility():
    result=read(REPORT/'candidate_selection.json')
    if result['selected_candidate'] is None:
        assert result['status']=='FAIL' and not result['eligible_candidates']
    else:
        assert result['selected_candidate'] in result['eligible_candidates']
