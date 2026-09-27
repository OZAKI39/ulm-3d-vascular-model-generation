"""Load the project-local binary runtime before importing the transfer tools."""
import ctypes
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def prepare_runtime():
    library = ROOT/'outputs/topbrain_brava_transfer/runtime/lib/libusb-1.0.so.0'
    if library.exists():
        ctypes.CDLL(str(library),mode=ctypes.RTLD_GLOBAL)


def require_versions():
    import importlib.metadata as m
    for name,version in [('open3d','0.20.0')]:
        if m.version(name)!=version:raise RuntimeError(f'Run tools/setup_topbrain_brava_transfer.py: expected {name} {version}')
    prepare_runtime()
