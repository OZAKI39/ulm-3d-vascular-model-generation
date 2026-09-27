import pytest
import numpy as np
@pytest.mark.parametrize('index',[0,1])
def test_independent_100_steps(parity_cases,index):
    d=parity_cases[index]
    assert d['steps']==100 and d['independently_evolved']
    assert len(d['errors'])==100 and d['accepted_substeps_bridge']==d['accepted_substeps_standalone']
    assert all(e['exact_equal'] and e['neighbor_mismatch_count']==0 for e in d['errors'])
    assert any(not np.array_equal(a['position'],b['position']) for a,b in zip(d['bridge'][0]['particles'],d['bridge'][-1]['particles']))
