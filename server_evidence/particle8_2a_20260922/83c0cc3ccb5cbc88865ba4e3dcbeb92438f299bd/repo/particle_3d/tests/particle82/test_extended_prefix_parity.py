import json,os
from pathlib import Path
import numpy as np

def test_every_extended_accepted_prefix(evidence_root):
    original=Path(os.environ.get('PARTICLE82_PREVIOUS_ROOT','/root/particle8_2_runs/readonly_inputs/p81/particle_3d/outputs/particle8_1'))
    s=json.loads((evidence_root/'continuation/EXTENDED_GUARD_AUDIT.json').read_text())
    assert [r['guard_factor'] for r in s['factors']]==[4,16]
    for trial in s['factors']:
        assert trial['all_original_prefixes_exact'] and len(trial['rows'])==720
        assert trial['completed_from_stops']+trial['still_stopped']==720
        for r in trial['rows']:
            pid=r['particle_id'];x=np.load(original/'trajectories'/f'mb_{pid:06d}.npz')['samples']
            y=np.load(evidence_root/f"continuation/guard_x{trial['guard_factor']}/trajectories/mb_{pid:06d}.npz")['samples']
            assert len(y)>=len(x) and y[:len(x)].tobytes()==x.tobytes()
            meta=json.loads((evidence_root/f"continuation/guard_x{trial['guard_factor']}/trajectories/mb_{pid:06d}.json").read_text())
            assert meta['continuation_checkpoint']['saved_prefix_count']==len(x)
            assert meta['continuation_checkpoint']['historical_states_not_reintegrated']
            if len(y)>len(x):
                suffix=y[len(x)-1:]
                residual=np.diff(suffix[:,1:4],axis=0)-np.diff(suffix[:,0])[:,None]*suffix[1:,4:7]
                assert np.max(np.abs(residual))<1e-18
                assert np.all(np.diff(suffix[:,0])>0)
