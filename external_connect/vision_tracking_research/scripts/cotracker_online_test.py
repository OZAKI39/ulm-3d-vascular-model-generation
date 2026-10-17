"""Official CoTracker3 Online, only arrived frames, max 16 input images, mock admission."""
from pathlib import Path
import sys,time,json,threading,collections,subprocess,hashlib,traceback
ROOT=Path(__file__).resolve().parents[1];sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'third_party/co-tracker'))
from replay_benchmark import ReplayStream,CFG,diststats,dump,DT
import torch,numpy as np,cv2
from cotracker.predictor import CoTrackerOnlinePredictor
OUT=ROOT/'results/cotracker_online';OUT.mkdir(exist_ok=False)
records=[];windows=collections.deque(maxlen=16);lock=threading.Lock();stop=threading.Event();sample_error=[]
try:
 torch.set_num_threads(2);torch.cuda.reset_peak_memory_stats();model=CoTrackerOnlinePredictor(checkpoint=str(ROOT/'third_party/checkpoints/scaled_online.pth')).cuda().eval()
 ann=json.loads((ROOT/'evidence/visual_annotations.json').read_text())['positive_only_keyframes']['120'];query=np.asarray(ann)-[650,350];query=torch.tensor(np.c_[np.zeros(len(query)),query],dtype=torch.float32,device='cuda')[None]
 stream=ReplayStream(120,160).start()
 def sample():
  last=None
  try:
   while not stop.is_set():
    p=stream.snapshot()
    if p and p[0]!=last:
     fid,rx,utc,bgr=p; rgb=cv2.cvtColor(bgr[350:700,650:830],cv2.COLOR_BGR2RGB)
     with lock:windows.append((fid,rx,rgb))
     last=fid
    if stream.done.is_set():break
    time.sleep(.0005)
  except Exception as e:sample_error.append(repr(e))
 th=threading.Thread(target=sample,daemon=True);th.start();start=time.monotonic();last_end=None;session=0
 while time.monotonic()-start<300:
  with lock:batch=list(windows)
  if len(batch)<16 or (last_end is not None and batch[-1][0]<last_end+8):
   if stream.done.is_set():break
   time.sleep(.001);continue
  ids=[a[0] for a in batch]
  if ids[-1]-ids[0]!=15:
   records.append({'error':'NONCONTIGUOUS_RECEIVED_WINDOW','frames':ids});break
  reset=last_end is None or ids[-1]-last_end!=8
  rgb=torch.from_numpy(np.stack([a[2] for a in batch])).permute(0,3,1,2)[None].float().cuda()
  if reset:
   # Initial queries belong to frame120. After an overrun re-seed current-window
   # image features; do not silently keep an identity across the gap.
   if ids[0]!=120:
    pts=cv2.goodFeaturesToTrack(cv2.cvtColor(batch[0][2],cv2.COLOR_RGB2GRAY),maxCorners=9,qualityLevel=.05,minDistance=8)
    query=torch.tensor(np.c_[np.zeros(len(pts)),pts[:,0]],device='cuda',dtype=torch.float32)[None]
   model(rgb[:,:1],is_first_step=True,queries=query);session+=1
  torch.cuda.synchronize();t=time.monotonic()
  with torch.inference_mode():tracks,visible=model(video_chunk=rgb,is_first_step=False)
  torch.cuda.synchronize();done=time.monotonic();idx=tracks.shape[1]-1;positions=tracks[0,idx].cpu().numpy()+[650,350];visibility=visible[0,idx].cpu().numpy()
  age=(done-batch[-1][1])*1000
  rec={'frame_id':ids[-1],'input_frame_ids':ids,'session':session,'reset':reset,'inference_ms':(done-t)*1000,'latest_frame_age_ms':age,'oldest_frame_age_ms':(done-batch[0][1])*1000,'time_since_replay_start_ms':(done-start)*1000,'returned_history_frames':int(tracks.shape[1]),'visible_points':int(visibility.sum()),'n_queries':len(visibility),'xy_latest_px':positions.tolist(),'visible_latest':visibility.tolist(),'fresh':age<=100,'bubble_control_admitted':False,'max_input_images':16,'output_semantics':'released at actual inference completion; no retroactive real-time observations'}
  records.append(rec);print(json.dumps({k:rec[k] for k in ['frame_id','inference_ms','latest_frame_age_ms','visible_points','reset']}),flush=True);last_end=ids[-1]
 stop.set();th.join(3);stream.stop();measured=[r for r in records if 'inference_ms' in r]
 summary={'status':'RUN','test_kind':'CAUSAL_REPLAY_TEST','online_only':True,'repo_head':subprocess.check_output(['git','-C',str(ROOT/'third_party/co-tracker'),'rev-parse','HEAD'],text=True).strip(),'execution_git_head':subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'torch_version':torch.__version__,'gpu':torch.cuda.get_device_name(0),'roi_xywh':[650,350,180,350],'model_interpolation_shape':list(model.interp_shape),'window':16,'step':model.step,'input_requested_frames':[120,160],'decoded_frames':len(stream.records),'total_wall_s':time.monotonic()-start,'emissions':len(measured),'fresh_emissions':sum(r['fresh'] for r in measured),'input_window_span_ms':15*DT*1000,'nominal_update_spacing_ms':8*DT*1000,'inference_ms':diststats([r['inference_ms'] for r in measured]),'latest_frame_age_ms':diststats([r['latest_frame_age_ms'] for r in measured]),'peak_cuda_allocated_mb':torch.cuda.max_memory_allocated()/1e6,'peak_cuda_reserved_mb':torch.cuda.max_memory_reserved()/1e6,'sample_errors':sample_error,'CAMERA_TO_CONTROLLER_END_TO_END_LATENCY':'NOT_MEASURED','caveat':'Tiny bounded contrast; n-window quantiles descriptive only, not a broad accuracy benchmark. Queries are image landmarks, not certified bubbles.'}
 dump(OUT/'summary.json',summary);dump(OUT/'windows.json',records);dump(OUT/'decoder.json',stream.records)
except Exception as e:
 dump(OUT/'failure.json',{'error':repr(e),'traceback':traceback.format_exc(),'partial_windows':records});raise
