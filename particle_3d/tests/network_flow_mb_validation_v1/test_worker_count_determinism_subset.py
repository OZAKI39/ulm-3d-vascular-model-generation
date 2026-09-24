import gzip,json
import numpy as np
def test_worker_count_determinism_subset(report):
    for pid in [3,15,22]:
        ps=[report/'outputs'/name/'trajectories'/f'mb_{pid:06d}' for name in ['determinism_w1','determinism_w3','NEW']]
        data=[np.load(p.with_suffix('.npz'))['samples'] for p in ps]
        assert all(a.shape==data[0].shape and a.tobytes()==data[0].tobytes() for a in data)
        logs=[gzip.decompress(p.with_suffix('.audit.jsonl.gz').read_bytes()) for p in ps]
        assert logs[0]==logs[1]==logs[2]
        meta=[json.loads(p.with_suffix('.json').read_text()) for p in ps]
        for k in ['end_reason','failure_detail','exit_outlet','residence_time_s']:assert meta[0][k]==meta[1][k]==meta[2][k]
