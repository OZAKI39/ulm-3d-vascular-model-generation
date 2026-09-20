import json
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
def artifact(name):
    path=ROOT/'reports/sv1_3'/(name+'.json')
    assert path.exists(), 'Missing SV1.3 evidence: '+name
    return json.loads(path.read_text())
def gpu_evidence(name):
    env=artifact('gpu_environment')
    path=ROOT/'reports/sv1_3'/(name+'.json')
    if not path.exists():
        decision=artifact('gpu_decision')
        assert decision['status'] in ('GPU_BLOCKED','GPU_BACKEND_NOT_REACHABLE')
        pytest.skip('Hardware/backend unavailable, recorded evidence: '+decision['reason'])
    data=json.loads(path.read_text())
    if data.get('status') in ('NOT_RUN','BLOCKED'):
        assert data.get('reason')
        pytest.skip(data['reason'])
    return data
