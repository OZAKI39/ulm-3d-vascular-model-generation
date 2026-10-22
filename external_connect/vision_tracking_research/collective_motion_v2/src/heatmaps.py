"""Actual frozen OpenCV centers; Gaussian rasterization is not ground truth."""
import cv2
import numpy as np
from legacy import CandidateDetector
from vcommon import P

def centers_to_heatmap(centers,shape,sigma=3.):
    h,w=shape;impulse=np.zeros(shape,np.float32)
    for x,y in np.asarray(centers).reshape(-1,2):
        if 0<=x<w and 0<=y<h:
            impulse[min(h-1,round(y)),min(w-1,round(x))]+=1
    heat=cv2.GaussianBlur(impulse,(0,0),sigma,borderType=cv2.BORDER_CONSTANT)*(2*np.pi*sigma*sigma)
    return np.clip(heat,0,1)

class HeatmapDetector:
    def __init__(self):self.detector=CandidateDetector()
    def __call__(self,image):
        objects,gray=self.detector(image)
        centers=np.array([x['xy'] for x in objects],float).reshape(-1,2)
        heat=centers_to_heatmap(centers,gray.shape,P['heatmap']['sigma_raw_px'])
        return objects,gray,heat
