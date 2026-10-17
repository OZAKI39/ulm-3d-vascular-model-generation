from pathlib import Path
import cv2,numpy as np,json,time
ROOT=Path(__file__).resolve().parents[1]; cv2.setNumThreads(1)
im=cv2.imread(str(ROOT/'evidence/frame_00000.png')); x,y,w,h=640,80,620,960; crop=im[y:y+h,x:x+w]; gray=cv2.cvtColor(crop,cv2.COLOR_BGR2GRAY)
bg=cv2.GaussianBlur(gray,(0,0),4); contrast=cv2.subtract(bg,gray); mask=np.uint8(contrast>=6)*255
n,lab,stats,centers=cv2.connectedComponentsWithStats(mask,8)
objects=[]; drawn=crop.copy()
for k in range(1,n):
 xx,yy,ww,hh,area=stats[k]
 if 6<=area<=800 and max(ww/hh,hh/ww)<=3.5 and min(xx,yy)>2 and xx+ww<w-2 and yy+hh<h-2:
  objects.append({'center':centers[k].tolist(),'area':int(area),'box':[int(v) for v in stats[k][:4]]})
  cv2.rectangle(drawn,(xx,yy),(xx+ww,yy+hh),(0,190,255),1)
cv2.imwrite(str(ROOT/'evidence/calibration_detections.png'),drawn)
(ROOT/'evidence/calibration.json').write_text(json.dumps({'frame':0,'roi':[x,y,w,h],'sigma':4,'threshold':6,'minarea':6,'maxarea':800,'n':len(objects),'objects':objects},indent=2)); print(len(objects))
