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
