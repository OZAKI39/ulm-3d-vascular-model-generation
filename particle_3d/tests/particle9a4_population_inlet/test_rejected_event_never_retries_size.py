from .conftest import RejectChecker

def test_size_draw_once(simple_source):
 simple_source.checker=RejectChecker()
 for i in range(1,18):simple_source.proposal(i)
 assert simple_source.distribution.calls==simple_source.sampler.calls==simple_source.checker.calls==17
