import multiprocessing as mp
from particle_3d.particle82a_admission import common_event,method_a,method_b

_CTX=None


def evaluate(i):
    event=common_event(i,ctx=_CTX)
    return method_a(event,_CTX,guard=16),method_b(event,_CTX,guard=16)


def test_parallel_results_match_serial(simple_context):
    global _CTX
    _CTX=simple_context
    serial=list(map(evaluate,range(1,9)))
    with mp.get_context('fork').Pool(2) as pool:parallel=pool.map(evaluate,range(1,9))
    assert serial==parallel
