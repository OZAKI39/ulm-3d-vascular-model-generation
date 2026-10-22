"""Development/validation before freeze; final is inaccessible without freeze."""
import sys,time,argparse,traceback,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
from backends import Backend
from studies import multipass,score
import cv2
def main():
 p=argparse.ArgumentParser();p.add_argument('split',choices=['development','validation','final_test']);a=p.parse_args();cv2.setNumThreads(2)
 protocol=json.loads((ROOT/'configs/research_protocol.json').read_text());index=json.loads((OUT/'data'/a.split/'index.json').read_text());final=a.split=='final_test'
 if final:
  freeze=json.loads((OUT/'results/CONFIGURATION_FREEZE.json').read_text());assert not (OUT/'results/final_test_opened.json').exists(),'FINAL_ALREADY_OPENED'
  dump(OUT/'results/final_test_opened.json',dict(identity=identity(),freeze_sha256=sha(OUT/'results/CONFIGURATION_FREEZE.json')))
 engines={n:Backend(n) for n in (['B0','B1','B2','B3','B4'] if final else ['B2','B3','B4'])};rows=[];fails=[];mon=ResourceMonitor()
 for meta in index:
  data=np.load(OUT/meta['file']);f=data['frames'];background=f[0].astype(np.float32);backgrounds=[]
  for image in f:backgrounds.append(background.copy());cv2.accumulateWeighted(image,background,.08)
  configs=[]
  if final:
   configs=[(n,freeze['preprocess'],freeze['gap'],48,False) for n in engines]
   configs += [('B4','local',2,48,False),('B4',freeze['preprocess'],freeze['gap'],32,True)]
  else:
   configs=[('B4',mode,gap,48,False) for mode in protocol['modes'] for gap in protocol['gaps']]
   configs += [(n,'local',2,32,True) for n in ['B2','B3','B4']]
   configs += [('B4','local',2,ws,False) for ws in [32,64]]
  for name,mode,gap,ws,multi in configs:
   start,end=15-gap,15;engine=engines[name];engine.ws=ws;engine.ov=ws//2
   try:
    bg0=backgrounds[start] if mode=='causal' else backgrounds[0];bg1=backgrounds[end] if mode=='causal' else backgrounds[0]
    kw=dict(mode=mode,background_a=bg0,background_b=bg1)
    r=multipass(engine,f[start],f[end],float(data['times'][end]-data['times'][start]),**kw) if multi else engine.run(f[start],f[end],float(data['times'][end]-data['times'][start]),actual_coordinates=True,**kw)
    result=score(r,data,meta,start,end);rows.append(dict(split=a.split,case=meta['case'],seed=meta['seed'],backend=name,mode=mode,gap=gap,ws=ws,multipass=multi,**result))
   except Exception as e:
    fails.append(dict(case=meta['case'],seed=meta['seed'],backend=name,mode=mode,gap=gap,ws=ws,multipass=multi,error=repr(e),traceback=traceback.format_exc()))
  print(a.split,meta['case'],meta['seed'],'rows',len(rows),'failures',len(fails),flush=True)
 dump(OUT/f'results/accuracy_{a.split}.json',dict(identity=identity(),rows=rows,failures=fails,resources=mon.close(),scope=protocol['independence_scope']))
if __name__=='__main__':main()
