import sys,subprocess
from pathlib import Path
r=Path(__file__).resolve().parent
raise SystemExit(subprocess.call([sys.executable,'-m','pytest',str(r/'tests'),'-q',*sys.argv[1:]]))
