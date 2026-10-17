"""Non-destructive assembly of separately predicted SAM instances."""
import numpy as np


def assemble_candidates(local_masks, targets, choices, included, shape, occupied=None):
    occupancy=np.zeros(shape,np.uint16);items={};rejected={}
    for i,(masks,target,choice,keep) in enumerate(zip(local_masks,targets,choices,included)):
        if not keep:continue
        info=target['proposals'][choice]
        if not info['allowed']:
            rejected[i]='；'.join(info['reasons']);continue
        mask=np.asarray(masks[choice],bool);yy,xx=np.where(mask)
        y0,x0,_,_=target['crop_yxyx'];indices=(yy+y0,xx+x0)
        if not len(yy):rejected[i]='空候选';continue
        if occupied is not None and np.any(occupied[indices]):
            rejected[i]='与已有标签重叠；未覆盖';continue
        items[i]=indices;occupancy[indices]+=1
    # Reject both intersecting instances: never assign overlapping pixels by order.
    for i,indices in list(items.items()):
        if np.any(occupancy[indices]>1):
            rejected[i]='与另一单泡候选重叠；请切换候选或只保留一个';del items[i]
    labels=np.zeros(shape,np.uint32)
    for i,indices in items.items():labels[indices]=i+1
    return labels,items,rejected
