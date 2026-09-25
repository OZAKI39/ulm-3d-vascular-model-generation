#!/usr/bin/env python3
"""Offline verification of the existing Open3D runtime; installs nothing."""
from pathlib import Path
import sys
import json
import importlib.metadata as metadata
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from vascular_processing.transfer_runtime import require_versions

def main():
    require_versions()
    import open3d as o3d
    print(json.dumps(dict(python=sys.executable,open3d=metadata.version('open3d'),
        similarity_scaling=o3d.pipelines.registration.TransformationEstimationPointToPoint(with_scaling=True).with_scaling,
        packages_installed=False,network_used=False),indent=2))
    return 0

if __name__=='__main__':raise SystemExit(main())
