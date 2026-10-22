"""Freeze selection using development/validation only, then open fresh procedural test."""
import sys,json,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
p=OUT/'results/CONFIGURATION_FREEZE.json';assert not p.exists(),'REFUSE_REFREEZE'
record=dict(identity=identity(),backend='B4',reference_backend='B2',preprocess='causal',gap=3,window=48,overlap=24,background_alpha=.08,snr_min=1.2,texture_std_min=2.,full_channel_strategy='REGULAR_GPU_BATCH_EVERY_FRAME_WITH_SUPPORT_GATED_OUTPUT',schedule=1,selected=False,reason='Development/validation causal gap3 had 2 false accepts per22 versus local gap2 6 per22. Regular CuPy pilot faster than masked sparse after actual successful computation; no claim same coverage speedup from skipping. Multipass not chosen: no reliable aggregate false-accept improvement.',known_unresolved=['symmetric rotation bias exceeds vorticity limit','opposed local motion can pass affine fit','no independent real motion truth','3-frame span adds49.95ms midpoint age'],selection_inputs={f:sha(OUT/'results'/f) for f in ['accuracy_development.json','accuracy_validation.json']},test_index_sha256=sha(OUT/'data/final_test/index.json'),training_ready=False)
dump(p,record)
subprocess.run([sys.executable,str(ROOT/'scripts/accuracy_study.py'),'final_test'],check=True)
