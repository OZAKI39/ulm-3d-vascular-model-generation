"""One new-process cold initialization + first real cached-triple inference."""
import time
start=time.perf_counter()
import sys,json,argparse
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from vcommon import OUT,DATA,dump,identity
from inference import Predictor
import torch,numpy as np
p=argparse.ArgumentParser();p.add_argument('--model',required=True);a=p.parse_args();torch.set_num_threads(4)
imports_s=time.perf_counter()-start;f=json.loads((OUT/'results/model_freeze.json').read_text());k=f['champions_by_validation_only'][a.model];t=time.perf_counter();model=Predictor(OUT/f['models'][k]['checkpoint']);load_s=time.perf_counter()-t
x=np.load(DATA/'validation/gray.npy',mmap_mode='r')[0].astype(np.float32)/255;h=np.load(DATA/'validation/heat.npy',mmap_mode='r')[0].astype(np.float32)/255;d=np.load(DATA/'validation/dt.npy',mmap_mode='r')[0];images=torch.tensor(np.stack([(x-.65)/.25,h],1)[None]);t=time.perf_counter();out=model.batch(images,torch.tensor(d[None]));pred=model.outputs(out,d[None]);torch.cuda.synchronize();first_s=time.perf_counter()-t
record=dict(model=a.model,imports_s=imports_s,checkpoint_and_cuda_initialization_s=load_s,first_inference_with_transfer_s=first_s,total_since_script_start_s=time.perf_counter()-start,scope='new process, first inference; excludes interpreter spawn/SSH and video decode',identity=identity());dump(OUT/'results'/f'cold_process_{a.model}.json',record);print(json.dumps({k:record[k] for k in ['model','total_since_script_start_s']}))
