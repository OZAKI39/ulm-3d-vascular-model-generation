import sys,time,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from common import *
from distribution import CandidateDetector
from methods import piv
from affine import parameters
z=np.load(ROOT/'results/synthetic/test/7201_step.npz');ims=z['frames'];H=z['H'];det=CandidateDetector();oo=[det(i)[0] for i in ims];rows=[];mode=json.loads((ROOT/'results/selection.json').read_text())['mode']
for window in P['windows']:
    for end in range(window-1,len(ims)):
        begin=end-window+1;dt=(window-1)*P['synthetic_dt_s'];HH=H[end]@np.linalg.inv(H[begin]);truth=parameters(HH[:2,:2],HH[:2,2],dt,np.array([128.,128.]))
        r=piv(ims[begin],ims[end],dt,'local' if mode=='multipass' else mode,oo[begin],oo[end],multipass=mode=='multipass');r.pop('vectors',None)
        rows.append({'window':window,'end':end,'end_t':end*P['synthetic_dt_s'],'interval_truth':truth,'instantaneous_vx_after_step':30 if end>7 else 0,**r})
dump(ROOT/'results/window_response.json',{'onset_t':7*P['synthetic_dt_s'],'scope':'causal endpoint method; interval-average response cannot jump instantaneously','rows':rows})
