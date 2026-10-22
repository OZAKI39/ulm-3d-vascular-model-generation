"""Post-hoc failure description; never used to tune or choose a model."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from vcommon import OUT,DATA,dump
r=json.loads((OUT/'results/analysis.json').read_text());label=np.load(DATA/'final_test/label.npy');truth=np.load(DATA/'final_test/truth.npy');mask=label==0;diagnostic={}
for key in [r['champions'][n] for n in ['direct_gray','direct_fused','direct_heat','direct_fused_robust','geo_fused']]:
 p=np.load(OUT/'results'/f'{key}_final_test.npz')['mean'];diagnostic[key]=dict(prediction_std=np.std(p[mask],0).tolist(),truth_std=np.std(truth[mask],0).tolist(),std_ratio=(np.std(p[mask],0)/np.std(truth[mask],0)).tolist(),bias=np.mean(p[mask]-truth[mask],0).tolist())
dump(OUT/'results/posthoc_collapse_diagnostic.json',dict(scope='POST_HOC_FAILURE_DIAGNOSTIC_ONLY; no model selection/retraining/test tuning',models=diagnostic,zero_motion_MAE=np.mean(abs(truth[mask]),0).tolist()))
