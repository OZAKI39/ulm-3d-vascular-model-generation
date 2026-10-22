"""Two routes only: reused OpenCV components or Cellpose cyto3, then Norfair.

All current outputs are UNREVIEWED image regions. Semantic bubble/group admission
requires a later frozen expert evaluation; motion is never an identity proof.
"""
from pathlib import Path
import sys,os,time,collections,uuid
import numpy as np
import cv2
from scipy.ndimage import find_objects
from norfair import Detection,Tracker
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'scripts'))
import replay_benchmark as legacy
cv2.setNumThreads(1)

def stats(values):
    a=np.asarray(values,float)
    return {k:float(f(a)) for k,f in [('mean',np.mean),('p50',lambda a:np.percentile(a,50)),
                                    ('p95',lambda a:np.percentile(a,95)),('maximum',np.max)]} if len(a) else None

def region_shape(mask,origin=(0,0)):
    m=cv2.moments(np.uint8(mask),binaryImage=True);area=m['m00']
    if area<1:return None
    q=np.array([[m['mu20'],m['mu11']],[m['mu11'],m['mu02']]])/area
    eig,v=np.linalg.eigh(q);xy=np.array([m['m10'],m['m01']])/area+origin
    contours,_=cv2.findContours(np.uint8(mask),cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
    c=max(contours,key=cv2.contourArea).reshape(-1,2)+origin if contours else np.empty((0,2))
    return {'centroid_xy_px':xy.tolist(),'area_px2':area,'Rg_px':float(np.sqrt(np.trace(q))),
            'Q_px2':q.tolist(),'aspect_ratio':float(np.sqrt(max(eig[1],0)/max(eig[0],.1))),
            'major_axis_deg':float(np.rad2deg(np.arctan2(v[1,-1],v[0,-1]))%180),
            'contour_xy_px':c.tolist(),'semantics':'UNREVIEWED_IMAGE_REGION; not admitted group shape'}

class VisionEngine:
    def __init__(self,route,origin=(0,0),exclusions=(),detector_config=None):
        if route not in ('opencv','cellpose'):raise ValueError(route)
        self.route=route;self.origin=np.asarray(origin);self.exclusions=list(exclusions)
        self.detector_config=dict(detector_config or {})
        self.detector=legacy.Detector('cc') if route=='opencv' else None
        self.model=None;self.torch=None
        if route=='cellpose':
            os.environ['CELLPOSE_LOCAL_MODELS_PATH']=str(ROOT/'models/cellpose')
            import torch
            from cellpose import models
            self.torch=torch;torch.set_num_threads(2)
            weights=self.detector_config.get('weights')
            self.model=(models.CellposeModel(gpu=True,pretrained_model=weights) if weights
                        else models.CellposeModel(gpu=True,model_type='cyto3'))
            torch.cuda.reset_peak_memory_stats()
        self.previous_gray=None;self.previous_time=None;self.previous_seq=None;self.session=None
        self.history={};self.reset_events=[];self.tracker=None;self.generation=0
        self.reset('initial')
    def reset(self,reason):
        self.tracker=Tracker(distance_function='euclidean',distance_threshold=15,
                             initialization_delay=1,hit_counter_max=3,past_detections_length=5)
        self.history={};self.previous_gray=None;self.previous_time=None;self.previous_seq=None
        self.generation+=1;self.reset_events.append({'reason':reason,'generation':self.generation})
    def detect(self,bgr):
        h,w=bgr.shape[:2]
        if self.route=='opencv':
            legacy.CFG['roi_xywh']=[0,0,w,h]
            legacy.CFG['fixed_structure_exclusion_xyxy']=self.exclusions
            previous=legacy.CFG['contrast_threshold']
            legacy.CFG['contrast_threshold']=self.detector_config.get('contrast_threshold',previous)
            try:allowed,gray,labels,_=self.detector(bgr)
            finally:legacy.CFG['contrast_threshold']=previous
            objects=[{'label':o['label'],'xy_local':o['xy'],'area':o['area'],
                      'score':min(o['peak']/30,1),'score_semantics':'uncalibrated local contrast score'} for o in allowed]
            # Preserve large connected regions as unreviewed aggregate candidates.
            ids,counts=np.unique(labels,return_counts=True)
            for label,area in zip(ids,counts):
                if label and 450<area<=min(50000,.4*h*w):
                    yy,xx=np.where(labels==label)
                    objects.append({'label':int(label),'xy_local':[float(xx.mean()),float(yy.mean())],
                                    'area':int(area),'score':None,'score_semantics':'no calibrated aggregate confidence'})
        else:
            gray=cv2.cvtColor(bgr,cv2.COLOR_BGR2GRAY)
            labels,flows,_=self.model.eval(bgr,channels=[0,0],diameter=self.detector_config.get('diameter',30),augment=False,resample=True,
                                           flow_threshold=.4,cellprob_threshold=0,min_size=self.detector_config.get('min_size',15))
            self.torch.cuda.synchronize()
            ids,counts=np.unique(labels,return_counts=True);objects=[]
            probability=flows[2]
            for label,area in zip(ids,counts):
                if not label:continue
                region=labels==label;yy,xx=np.where(region)
                raw=float(np.median(probability[region]))
                objects.append({'label':int(label),'xy_local':[float(xx.mean()),float(yy.mean())],
                                'area':int(area),'score':float(1/(1+np.exp(-np.clip(raw,-50,50)))),
                                'score_semantics':'sigmoid Cellpose cell probability; not calibrated bubble probability'})
        return objects,gray,labels
    def process(self,bgr,header):
        t0=time.monotonic();seq=header['seq'];source_time=header['source_time_s'];session=header['session']
        gap=None if self.previous_seq is None else seq-self.previous_seq
        if self.session is not None and session!=self.session:self.reset('source_session_changed')
        elif self.previous_time is not None and (source_time<=self.previous_time or source_time-self.previous_time>.15 or gap>5):
            self.reset('time_or_frame_gap')
        self.session=session
        objects,gray,labels=self.detect(bgr);t1=time.monotonic()
        det=[]
        for obj in objects:
            obj['xy_px']=(np.asarray(obj['xy_local'])+self.origin).tolist()
            obj['candidate_class']='SINGLE_BUBBLE_CANDIDATE' if obj['area']<=450 else 'AGGREGATE_CANDIDATE'
            obj['frame_id']=seq
            det.append(Detection(np.array([obj['xy_px']]),data=obj))
        # Norfair's original association/filter; no fixed direction, speed or F10S prior.
        tracked=self.tracker.update(detections=det);t2=time.monotonic()
        targets=[]
        for tr in tracked:
            obj=tr.last_detection.data;observed=obj['frame_id']==seq
            key=f'{session}:{self.generation}:{tr.id}'
            hist=self.history.setdefault(key,collections.deque(maxlen=5));velocity=None
            if observed:
                xy=np.array(obj['xy_px']);
                if hist and source_time>hist[-1][0]:velocity=((xy-hist[-1][1])/(source_time-hist[-1][0])).tolist()
                hist.append((source_time,xy.copy()))
            targets.append({'id':key,'status':'OBSERVED_UNREVIEWED' if observed else 'LOST',
                            'xy_px':obj['xy_px'] if observed else None,'velocity_px_per_source_s':velocity if observed else None,
                            'candidate_class':obj['candidate_class'],'score':obj['score'],'score_semantics':obj['score_semantics'],
                            'short_history':[[float(t),*p.tolist()] for t,p in hist] if observed else [],
                            'bubble_identity_verified':False,'physical_identity_guaranteed':False})
        shapes=[]
        boxes=find_objects(labels)
        for obj in objects:
            box=boxes[obj['label']-1]
            if box is None:continue
            ys,xs=box
            sh=region_shape(labels[box]==obj['label'],self.origin+[xs.start,ys.start])
            shapes.append({'label':obj['label'],'candidate_class':obj['candidate_class'],**sh})
        t3=time.monotonic();flow=None
        if self.previous_gray is not None and self.previous_gray.shape==gray.shape and seq%3==0:
            pts=cv2.goodFeaturesToTrack(self.previous_gray,maxCorners=60,qualityLevel=.03,minDistance=8)
            if pts is not None:
                dst,ok,err=cv2.calcOpticalFlowPyrLK(self.previous_gray,gray,pts,None,winSize=(15,15),maxLevel=2)
                back,bok,_=cv2.calcOpticalFlowPyrLK(gray,self.previous_gray,dst,None,winSize=(15,15),maxLevel=2)
                keep=(ok[:,0]>0)&(bok[:,0]>0)&(np.linalg.norm(back[:,0]-pts[:,0],axis=1)<1)&(err[:,0]<20)
                vectors=dst[:,0][keep]-pts[:,0][keep]
                flow={'points':int(keep.sum()),'median_displacement_px':np.median(vectors,axis=0).tolist() if len(vectors) else None,
                      'semantics':'local image motion, not liquid velocity'}
        t4=time.monotonic()
        self.previous_gray=gray;self.previous_time=source_time;self.previous_seq=seq
        self.history={k:v for k,v in self.history.items() if k in {t['id'] for t in targets}}
        message={'schema':2,'source':header['source'],'session':session,'frame_id':seq,'source_time_s':source_time,
                 'source_time_semantics':header['source_time_semantics'],'exposure_timestamp':None,
                 'timestamp_status':'EXPOSURE_TIMESTAMP_UNKNOWN','arrival_mono':header['arrival_mono'],
                 'age_upper_bound_at_arrival_ms':header['age_upper_bound_at_arrival_ms'],
                 'route':self.route,'roi_xywh':header['roi_xywh'],'candidates':objects,'targets':targets,
                 'unreviewed_region_shapes':shapes,'group_shape':None,'group_shape_status':'GROUP_SHAPE_NOT_ADMITTED',
                 'group_neighborhood':None,'apparent_incorporation_events':[],
                 'flow':flow,'status':'NO_VALID_TARGET','reason':'EXPERT_GROUND_TRUTH_PENDING',
                 'visual_candidate_available':bool(objects),'biological_valid':False,
                 'control_permission':False,'CAMERA_TO_CONTROLLER_END_TO_END_LATENCY':'NOT_MEASURED'}
        timing={'detect_ms':(t1-t0)*1000,'track_ms':(t2-t1)*1000,'morphology_ms':(t3-t2)*1000,
                'flow_ms':(t4-t3)*1000,'engine_ms':(t4-t0)*1000}
        return message,timing,labels

def neighborhood_for_admitted_masks(single_xy,group_mask,radius_px=60):
    """Geometry helper; caller must pass an independently admitted group mask."""
    shape=region_shape(group_mask)
    if shape is None:return None
    center=np.array(shape['centroid_xy_px']);xy=np.asarray(single_xy,float).reshape(-1,2)
    rel=xy-center;keep=np.linalg.norm(rel,axis=1)<=radius_px
    return {'radius_px':radius_px,'single_count':int(keep.sum()),'relative_xy_px':rel[keep].tolist()}
