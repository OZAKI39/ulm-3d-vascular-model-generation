import numpy as np
import pytest

from mevo_mb.acoustics.kwave_backend import KWaveCW, extract_complex
from mevo_mb.utils.config import load_config


def test_absolute_complex_phase():
    f = 2e6
    t = (np.arange(128) + 1013) / (32*f)
    phasor = np.array([2+3j, -4+1j])
    signals = np.real(phasor[:, None]*np.exp(2j*np.pi*f*t))
    assert np.allclose(extract_complex(signals, t, f), phasor, atol=1e-6)


@pytest.mark.acoustic
def test_cpu_kwave(tmp_path):
    c = load_config("configs/smoke.yaml")
    c["acoustics"].update(backend="python", device="cpu")
    adapter = KWaveCW(c, tmp_path)
    p, g, meta = adapter.run(np.array([1000.]), "cpu")
    assert np.isfinite(p).all() and np.abs(p).max() > 1
    assert np.isfinite(g).all() and meta["device"] == "cpu"


@pytest.mark.gpu
def test_gpu_kwave(tmp_path):
    c = load_config("configs/smoke.yaml")
    adapter = KWaveCW(c, tmp_path)
    p, _, meta = adapter.run(np.array([1000.]), "gpu")
    assert np.abs(p).max() > 1 and meta["device"] == "gpu"
