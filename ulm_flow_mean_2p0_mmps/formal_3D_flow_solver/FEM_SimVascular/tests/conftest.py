import json
import gzip
from pathlib import Path
import pytest
from sv_validation.postprocess import SolutionMeasurements
ROOT = Path(__file__).resolve().parents[1]
@pytest.fixture(scope="session")
def root(): return ROOT
@pytest.fixture(scope="session")
def report():
    def load(name):
        p=ROOT / "reports/sv1" / (name + ".json")
        return json.loads(p.read_text() if p.exists() else gzip.decompress(p.with_suffix(".json.gz").read_bytes()))
    return load
@pytest.fixture(scope="session")
def measured(report):
    q = report("flow_qc")
    p = json.loads((ROOT / "configs/time_policy.json").read_text())
    m = SolutionMeasurements(ROOT / "outputs/sv1/SV_MESH/mesh_arrays.npz", q["Q_target_m3_s"], p["Umean_m_s"])
    u, pressure = m.read(ROOT / q["path"])
    return m, u, pressure, m.measure(u, pressure)
