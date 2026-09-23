from pathlib import Path
import json,hashlib,numpy as np
from audit_math import align_saved
ROOT=Path(__file__).resolve().parents[1]

def test_all_input_trajectories_hashes_unchanged():
    for p,h in json.loads((ROOT/'data/readonly_baseline.json').read_text()).items():
        if '/trajectories/' in p or '/imported_P82A/' in p:
            assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p

def test_velocity_timestamp_alignment_and_no_mutation():
    inv=json.loads((ROOT/'data/trajectory_inventory.json').read_text())
    r=next(r for r in inv['records'] if r['sample_count']>10)
    with np.load(r['local_path']) as f:s=f['samples'].copy()
    before=s.copy();x=align_saved(s)
    np.testing.assert_array_equal(s,before)
    np.testing.assert_array_equal(x[1:],s[:-1,1:4])
