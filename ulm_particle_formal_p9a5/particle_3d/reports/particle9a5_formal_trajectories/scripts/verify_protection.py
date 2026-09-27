"""Verify retained current science, excluding explicitly deleted history."""
from pathlib import Path
import subprocess, sys
ROOT = Path(__file__).resolve().parents[4]
R = ROOT/'particle_3d/reports/particle9a5_formal_trajectories'
subprocess.run([sys.executable,str(ROOT/'scripts/verify_current_data.py'),
                '--output',str(R/'data/protected_verification.json')],check=True)
