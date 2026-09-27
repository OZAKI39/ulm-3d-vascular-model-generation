"""Full current workflow regression, including trajectory parity re-integration."""
from pathlib import Path
import json, shutil, subprocess, sys
ROOT = Path(__file__).resolve().parents[4]
R = ROOT/'particle_3d/reports/particle9a5_formal_trajectories'
out = R/'logs/final_regression'
subprocess.run([sys.executable,str(ROOT/'scripts/check_current.py'),
                '--output',str(out),'--full-simulation'],check=True)
receipt=json.loads((out/'summary.json').read_text())
shutil.copy2(out/'tests.xml',R/'logs/final_tests.xml')
(R/'logs/final_tests.txt').write_text(json.dumps(receipt,indent=2)+'\n'+(out/'tests.log').read_text())
