"""Reuse the frozen OpenCV detector. No tracker is instantiated; no persistent IDs."""
import sys
import cv2
import numpy as np
from scipy.ndimage import gaussian_filter
from common import ROOT
sys.path.insert(0,str(ROOT.parent/'scripts'))
import replay_benchmark as legacy

class CandidateDetector:
    def __init__(self): self.detector=legacy.Detector('cc')
    def __call__(self, image):
        if image.ndim==2: image=cv2.cvtColor(image,cv2.COLOR_GRAY2BGR)
        h,w=image.shape[:2]
        # Detector configuration is scoped and restored; no old workspace files change.
        old=legacy.CFG.copy()
        try:
            legacy.CFG.update(roi_xywh=[0,0,w,h],fixed_structure_exclusion_xyxy=[],contrast_threshold=10)
            objects,gray,labels,mask=self.detector(image)
        finally: legacy.CFG.clear(); legacy.CFG.update(old)
        for o in objects:
            x,y=o['xy']
            o['near_roi_boundary']=min(x,y,w-1-x,h-1-y)<16
            o['quality_score']=min(1.,o['peak']/40.)  # contrast heuristic, NOT calibrated probability
        return objects, gray

def distribution(objects, shape):
    h,w=shape[:2]; xy=np.array([o['xy'] for o in objects],float).reshape(-1,2)
    density=np.zeros((h,w),np.float32)
    for x,y in xy: density[min(h-1,max(0,round(y))),min(w-1,max(0,round(x)))]+=1
    density=gaussian_filter(density,12,mode='reflect')
    out={'scope':'VISIBLE_2D_CANDIDATE_DISTRIBUTION_NOT_PHYSICAL_CLUSTER',
         'visible_count':len(xy),'valid_area_fraction':1.0,'valid_area_scope':'geometric ROI; optical focus validity separate',
         'density_summary':{'mean_per_px2':len(xy)/(w*h),'peak_per_px2':float(density.max()),'bandwidth_px':12},
         'ingress_egress_indicators':{'boundary_candidate_fraction':np.mean([o['near_roi_boundary'] for o in objects]) if objects else None,'count_change_is_not_material_flux':True},
         'group_center_xy':None,'covariance_2x2':None,'gyration_radius':None,'distribution_axis_ratio':None,
         'distribution_orientation':None,'principal_axis':None,'eigenvalues':None,
         'field_validity':{'group_center_xy':bool(len(xy)),'distribution_orientation':False},'sparse':len(xy)<8}
    if len(xy): out['group_center_xy']=xy.mean(0).tolist()
    if len(xy)>=3:
        cov=np.cov(xy.T,bias=True); eig,vec=np.linalg.eigh(cov)
        ratio=np.sqrt(eig[-1]/max(eig[0],1e-12)); orient=ratio>=1.15
        out.update(covariance_2x2=cov.tolist(),gyration_radius=float(np.sqrt(np.trace(cov))),distribution_axis_ratio=float(ratio),eigenvalues=eig.tolist(),
                   distribution_orientation=float(np.arctan2(vec[1,-1],vec[0,-1])%np.pi) if orient else None,
                   principal_axis=vec[:,-1].tolist() if orient else None)
        out['field_validity']['distribution_orientation']=bool(orient)
    return out,density
