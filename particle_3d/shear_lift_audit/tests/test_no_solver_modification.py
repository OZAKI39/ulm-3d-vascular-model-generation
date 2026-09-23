from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]

def test_all_frozen_sources_and_flows_unchanged():
    for p,h in json.loads((ROOT/'data/readonly_baseline.json').read_text()).items():
        if '/src/' in p or '/frozen_reference/' in p or '/frozen_flow/' in p:
            assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
