"""Offline three-frame CLI. No device connection; rejected fields stay null."""
import sys,argparse,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import cv2
from vcommon import ROOT,OUT,dump
from inference import Predictor,estimate_from_prediction
from heatmaps import HeatmapDetector
from legacy import distribution
from state_v2 import new_state,to_v1
import jsonschema
p=argparse.ArgumentParser();p.add_argument('--images',nargs=3,required=True);p.add_argument('--pts',nargs=3,type=float,required=True);p.add_argument('--model',default='geo_fused');p.add_argument('--device',default='cuda');p.add_argument('--output',type=Path,required=True);p.add_argument('--legacy',action='store_true');a=p.parse_args()
f=json.loads((OUT/'results/model_freeze.json').read_text());key=f['champions_by_validation_only'][a.model];cal=f['models'][key];m=Predictor(OUT/cal['checkpoint'],device=a.device);det=HeatmapDetector();received=time.monotonic();gray=[];heat=[]
for path in a.images:
 image=cv2.imread(path)
 if image is None:raise ValueError('UNREADABLE_IMAGE:'+path)
 objects,g,h=det(image);gray.append(g);heat.append(h)
if len(set(g.shape for g in gray))!=1:raise ValueError('INCONSISTENT_ROI_SHAPES')
pred,timing=m.frames(gray,heat,a.pts);estimate=estimate_from_prediction(pred,cal,a.model,m.sha);dist,_=distribution(objects,gray[0].shape)
s=new_state('OFFLINE_THREE_FRAME_CLI','f10v_cli',2,a.pts[-1],received,[0,0,gray[0].shape[1],gray[0].shape[0]],dist,estimate,{},dict(frame_ids=[0,1,2],pts_s=a.pts,history_span_s=a.pts[-1]-a.pts[0],reference='interval-equivalent'),model_metadata=dict(name=a.model,weight_sha256=m.sha,calibration_version='F10V-synthetic-groups-v1'))
jsonschema.validate(s,json.loads((ROOT/'configs/collective_state_v2.schema.json').read_text()));dump(a.output,to_v1(s) if a.legacy else s);print(a.output)
