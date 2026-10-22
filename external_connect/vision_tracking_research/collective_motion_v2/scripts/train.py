"""Bounded formal training. Never opens calibration/final_test datasets."""
import sys,os,argparse,json,time,random
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np,torch
from torch.utils.data import DataLoader
from networks import build_model
from learning import CachedTriples,loss_fn
from vcommon import P,ROOT,OUT,DATA,dump,sha,identity,ResourceMonitor,SCALES,stats

def seed_all(seed):random.seed(seed);np.random.seed(seed);torch.manual_seed(seed);torch.cuda.manual_seed_all(seed)
def worker_seed(worker):np.random.seed(torch.initial_seed()%2**32)
def loader(split,batch,shuffle=False,corrupt=None):
 return DataLoader(CachedTriples(split,corrupt),batch_size=batch,shuffle=shuffle,num_workers=4,pin_memory=True,persistent_workers=True,prefetch_factor=2,worker_init_fn=worker_seed)
def move(batch):return {k:v.cuda(non_blocking=True) if torch.is_tensor(v) else v for k,v in batch.items()}

@torch.no_grad()
def validate(model,dl,amp):
 model.eval();pred=[];truth=[];labels=[];logits=[]
 for b in dl:
  b=move(b)
  with torch.autocast('cuda',dtype=torch.bfloat16,enabled=amp):out=model(b['images'],b['dt'])
  pred.append(out['mean'].float().cpu());truth.append(b['truth'].cpu());labels.append(b['label'].cpu());logits.append(out['logits'].float().cpu())
 pred=torch.cat(pred);truth=torch.cat(truth);labels=torch.cat(labels);logits=torch.cat(logits);mask=labels==0
 mae=(pred[mask]-truth[mask]).abs().mean();ce=torch.nn.functional.cross_entropy(logits,labels)
 return dict(score=float(mae+.2*ce),normalized_mae=float(mae),classification_loss=float(ce),classification_accuracy=float((logits.argmax(1)==labels).float().mean()),mae_parameters=((pred[mask]-truth[mask]).abs().mean(0)*torch.tensor(SCALES)).tolist())

def train(name,seed,batch,amp,epochs=None,resume=False,per_model=None):
 seed_all(seed);torch.set_num_threads(4);torch.backends.cuda.matmul.allow_tf32=True;torch.backends.cudnn.allow_tf32=True
 model=build_model(name).cuda();optimizer=torch.optim.AdamW(model.parameters(),lr=P['training']['lr'],weight_decay=P['training']['weight_decay'])
 epochs=epochs or P['training']['epochs_max'];scheduler=torch.optim.lr_scheduler.CosineAnnealingLR(optimizer,epochs)
 folder=OUT/'training'/f'{name}_seed{seed}';folder.mkdir(parents=True,exist_ok=True)
 train_dl=loader('train',batch,True,'train' if name=='direct_fused_robust' else None);val_dl=loader('validation',batch)
 start_epoch=0;best=float('inf');best_epoch=0;history=[]
 if resume and (folder/'last.pt').exists():
  c=torch.load(folder/'last.pt',weights_only=False,map_location='cpu');model.load_state_dict(c['model']);optimizer.load_state_dict(c['optimizer']);scheduler.load_state_dict(c['scheduler']);start_epoch=c['epoch']+1;best=c['best'];best_epoch=c['best_epoch'];history=c['history']
  torch.set_rng_state(c['torch_rng']);torch.cuda.set_rng_state_all(c['cuda_rng']);np.random.set_state(c['numpy_rng']);random.setstate(c['python_rng'])
 monitor=ResourceMonitor();torch.cuda.reset_peak_memory_stats();steps=0;stop='MAX_EPOCHS';wait_times=[];step_times=[];validation_times=[]
 config=dict(name=name,seed=seed,batch_size=batch,amp='BF16' if amp else 'FP32',epochs_max=epochs,initialization='random; no pretrained weights',data={s:sha(DATA/s/'manifest.json') for s in ['train','validation']},identity=identity(),parameters=sum(p.numel() for p in model.parameters()))
 dump(folder/'config.json',config)
 wall_limit=per_model or P['training']['per_model_wall_s']
 for epoch in range(start_epoch,epochs):
  model.train();totals=[];components=[];last=time.perf_counter()
  for raw in train_dl:
   now=time.perf_counter();wait_times.append(now-last);b=move(raw);optimizer.zero_grad(set_to_none=True)
   with torch.autocast('cuda',dtype=torch.bfloat16,enabled=amp):out=model(b['images'],b['dt']);loss,parts=loss_fn(out,b,local='no_local' not in name)
   if not torch.isfinite(loss):raise FloatingPointError(f'nonfinite loss epoch{epoch} step{steps}')
   loss.backward();grad=torch.nn.utils.clip_grad_norm_(model.parameters(),10);optimizer.step()
   totals.append(float(loss.detach()));components.append({k:float(v) for k,v in parts.items()});steps+=1
   # float(loss) synchronizes this measured full step; CUDA events not substituted for wall time.
   last=time.perf_counter();step_times.append(last-now)
  scheduler.step();t=time.perf_counter();metrics=validate(model,val_dl,amp);validation_times.append(time.perf_counter()-t)
  row=dict(epoch=epoch+1,train_loss=float(np.mean(totals)),components={k:float(np.mean([r[k] for r in components])) for k in components[0]},validation=metrics,elapsed_s=time.perf_counter()-monitor.start,learning_rate=optimizer.param_groups[0]['lr']);history.append(row)
  improved=metrics['score']<best-1e-5
  if improved:best=metrics['score'];best_epoch=epoch
  checkpoint=dict(model=model.state_dict(),optimizer=optimizer.state_dict(),scheduler=scheduler.state_dict(),epoch=epoch,best=best,best_epoch=best_epoch,history=history,config=config,torch_rng=torch.get_rng_state(),cuda_rng=torch.cuda.get_rng_state_all(),numpy_rng=np.random.get_state(),python_rng=random.getstate())
  torch.save(checkpoint,folder/'last.pt')
  if improved:torch.save(checkpoint,folder/'best.pt')
  if (epoch+1)%P['training']['checkpoint_every_epochs']==0:torch.save(checkpoint,folder/f'epoch_{epoch+1:03}.pt')
  dump(folder/'history.json',history);print(json.dumps(dict(run=folder.name,**row)),flush=True)
  if epoch+1>=P['training']['min_epochs'] and epoch-best_epoch>=P['training']['patience']:stop='VALIDATION_EARLY_STOP';break
  if time.perf_counter()-monitor.start>=wall_limit:stop='PRE_REGISTERED_WALL_BUDGET';break
 resource=monitor.close();summary=dict(config=config,epochs_completed=len(history),optimizer_steps=steps,stop_reason=stop,best_epoch=best_epoch+1,best_validation=best,resource=resource,peak_cuda_allocated_mb=torch.cuda.max_memory_allocated()/2**20,data_wait_s=stats(wait_times),step_s=stats(step_times),validation_s=validation_times,checkpoint_sha256={p.name:sha(p) for p in folder.glob('*.pt')})
 dump(folder/'summary.json',summary);return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--name');p.add_argument('--seed',type=int,default=1729);p.add_argument('--batch',type=int,default=96);p.add_argument('--epochs',type=int);p.add_argument('--all',action='store_true');p.add_argument('--resume',action='store_true');p.add_argument('--fp32',action='store_true');a=p.parse_args()
 if not torch.cuda.is_available():raise RuntimeError('CUDA_REQUIRED_FOR_FORMAL_TRAINING')
 if a.all:
  start=time.perf_counter();results=[];runs=[(n,s) for s in P['training_seeds'] for n in P['models']]+[('direct_fused_2frame',1729),('geo_no_local_supervision',1729)]
  for n,s in runs:
   remaining=P['training']['max_total_gpu_wall_s']-(time.perf_counter()-start)
   if remaining<180:
    results.append(dict(name=n,seed=s,status='NOT_RUN',reason='pre_registered_total_wall_budget'));continue
   folder=OUT/'training'/f'{n}_seed{s}'
   if a.resume and (folder/'summary.json').exists():results.append(json.loads((folder/'summary.json').read_text()));continue
   result=train(n,s,a.batch,not a.fp32,a.epochs,a.resume,min(P['training']['per_model_wall_s'],remaining));results.append(result)
   dump(OUT/'training/suite.json',dict(elapsed_s=time.perf_counter()-start,runs=results,identity=identity()))
 else:train(a.name,a.seed,a.batch,not a.fp32,a.epochs,a.resume)
