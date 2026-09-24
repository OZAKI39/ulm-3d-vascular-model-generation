import json
import numpy as np

def test_all_64_ids_and_every_worker_level(evidence_root):
    root=evidence_root/'scaling';audit=json.loads((root/'SERVER_SCALING_BENCHMARK.json').read_text())
    assert audit['all_exact'] and len(audit['fixed_ids'])==64
    assert [r['workers'] for r in audit['worker_scaling_results']]==[1,2,4,8,16]
    reference=json.loads((root/'workers_1/catalog.json').read_text())['records']
    for group in audit['parallel_vs_serial']:
        assert group['all_exact'] and len(group['ids'])==64
        folder=root/f"workers_{group['workers']}"
        assert json.loads((folder/'catalog.json').read_text())['records']==reference
        for pid in audit['fixed_ids']:
            name=f'mb_{pid:06d}.npz'
            x=np.load(root/'workers_1/trajectories'/name)['samples'];y=np.load(folder/'trajectories'/name)['samples']
            assert x.shape==y.shape and x.dtype==y.dtype and x.tobytes()==y.tobytes()
