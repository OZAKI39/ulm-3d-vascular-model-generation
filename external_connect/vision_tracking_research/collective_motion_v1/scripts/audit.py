import sys,subprocess,importlib.metadata,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
paths=[VIDEO,ROOT.parent/'scripts/replay_benchmark.py',ROOT.parent/'smart_online_tracking_v2/src/live_bridge.py',Path('/mnt/f/CHRONOS_control/chronos_common.py'),Path('/home/lzy/projects/external_connect/CHRONOS_control/chronos_common.py')]
for folder in (ROOT.parent/'smart_online_tracking_v2/results/expert_freezes').iterdir():
    if folder.is_dir(): paths.extend(p for p in folder.rglob('*') if p.is_file())
paths += [SIM/'results/straight_tube_exp1/f10s/a_navigation_v1/runs/X_PAIR_3MS_DT1/trajectory.h5', SIM/'results/straight_tube_exp1/f10s/analysis/FINAL_GATE.json']
dump(ROOT/'evidence/protected_inputs_before.json',{str(p):sha(p) for p in paths if p.exists()})
packages=['openpiv','pycpd','numpy','scipy','av','opencv-python-headless','norfair','h5py']
dump(ROOT/'evidence/initial_audit.json',{'execution':identity(),'versions':{k:importlib.metadata.version(k) for k in packages},'vision_baseline':'35755720a01de6ac6d6f14ac38e2a3ca84260236','simulation_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=SIM,text=True).strip(),'gpu':subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,memory.used,driver_version','--format=csv'],text=True),'cleanup':{'deleted_files':[],'reason':'No active temporary entry proven unused; historical untracked artifacts preserved.'}})
(ROOT/'evidence/git_status_before.txt').write_text(subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True))
(ROOT/'requirements-lock.txt').write_text(subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True))
print('audit complete')
