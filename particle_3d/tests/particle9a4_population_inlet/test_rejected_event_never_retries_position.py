from .conftest import RejectChecker

def test_position_draw_once(simple_source):
 simple_source.checker=RejectChecker()
 row=simple_source.proposal(1)
 assert row['status']=='REJECTED_SOURCE_EVENT'
 assert simple_source.sampler.calls==simple_source.checker.calls==1
