"""Explicit image-size ranges and whole-mask merging for unreviewed proposals."""
import numpy as np


def diameter_range(values):
    a=np.asarray(values,dtype=float)
    if a.shape!=(2,) or not np.isfinite(a).all() or not (4<=a[0]<a[1]<=100):
        raise ValueError('直径范围需满足 4 ≤ 最小直径 < 最大直径 ≤ 100 px；输入的是直径，不是半径')
    return a.tolist()


def size_plan(request):
    """Old scalar requests retain their original one-scale/filter semantics."""
    model=request['model']
    if model not in ('cyto3','cellpose4','microsam_lm'):raise ValueError('Unknown proposal model')
    if 'diameter_range_px' in request:
        if 'diameter_px' in request:raise ValueError('Do not mix scalar and range diameter requests')
        lo,hi=diameter_range(request['diameter_range_px'])
        # Missing policy means the old strict behaviour for historical requests.
        policy=request.get('diameter_filter','STRICT_RANGE')
        if policy not in ('STRICT_RANGE','REFERENCE_ONLY'):raise ValueError('Unknown diameter filter policy')
        strict=policy=='STRICT_RANGE'
        scales=[lo,float(np.sqrt(lo*hi)),hi] if model=='cyto3' else [None]
        return {'mode':'DIAMETER_RANGE','range_px':[lo,hi],'inference_diameters_px':scales,
                'diameter_filter':policy,
                'range_measure':'area_equivalent_diameter_px=2*sqrt(area_px/pi)',
                'range_applies_to':('inference_scales_and_filter' if model=='cyto3' else 'candidate_filter_only') if strict
                                   else ('inference_scales_only' if model=='cyto3' else 'advisory_only'),
                'max_extent_px':2.5*hi if strict else None,
                'max_roi_area_fraction':None if strict else .5}
    if 'diameter_filter' in request:raise ValueError('Diameter filter policy requires a diameter range')
    d=float(request.get('diameter_px',16))
    if not np.isfinite(d) or not 4<=d<=100:raise ValueError('Invalid expected diameter')
    return {'mode':'LEGACY_NOMINAL_DIAMETER','nominal_diameter_px':d,
            'inference_diameters_px':[d] if model=='cyto3' else [None],
            'range_applies_to':'legacy_nominal_envelope','max_extent_px':2.5*d}


def collect_candidates(raw_by_scale,relative_roi,global_roi,plan):
    """Keep original whole masks; overlapping alternatives are logged, never clipped.

    Ranking prefers the scale nearest the candidate's area-equivalent diameter.
    This is a deterministic editing heuristic, not a calibrated confidence score.
    """
    y0,x0,y1,x1=relative_roi;gy0,gx0,_,_=global_roi
    scales=plan['inference_diameters_px']
    if len(raw_by_scale)!=len(scales):raise ValueError('Inference scale count mismatch')
    labels=np.zeros((y1-y0,x1-x0),np.uint32);pending=[];rejected=[];objects=[]
    for scale_index,(raw,scale) in enumerate(zip(raw_by_scale,scales)):
        if raw.ndim!=2 or raw.dtype.kind not in 'ui' or not (0<=y0<y1<=raw.shape[0] and 0<=x0<x1<=raw.shape[1]):
            raise ValueError('Invalid raw candidate image or ROI')
        local=raw[y0:y1,x0:x1]
        for ident in np.unique(local):
            if not ident:continue
            yy,xx=np.where(local==ident);area=len(yy);equiv=float(2*np.sqrt(area/np.pi));extent=int(max(np.ptp(yy)+1,np.ptp(xx)+1))
            reasons=[]
            if area<5:reasons.append('too_small')
            if plan['mode']=='DIAMETER_RANGE':
                lo,hi=plan['range_px']
                if plan.get('diameter_filter','STRICT_RANGE')=='STRICT_RANGE':
                    if equiv<lo:reasons.append('below_diameter_range')
                    if equiv>hi:reasons.append('above_diameter_range')
                    if extent>plan['max_extent_px']:reasons.append('excessive_extent_for_single_bubble')
                elif area>plan['max_roi_area_fraction']*labels.size:
                    reasons.append('dominates_roi_requires_review')
            else:
                d=plan['nominal_diameter_px']
                if area>np.pi*(d*1.8/2)**2 or extent>plan['max_extent_px']:reasons.append('too_large_for_single_bubble')
            if yy.min()==0 or yy.max()==labels.shape[0]-1 or xx.min()==0 or xx.max()==labels.shape[1]-1:
                reasons.append('roi_border_requires_review')
            item={'raw_id':int(ident),'scale_index':scale_index,'inference_diameter_px':scale,
                  'area_px':area,'equivalent_diameter_px':equiv,'extent_px':extent,
                  'centroid_yx':[float(yy.mean()+gy0),float(xx.mean()+gx0)],'reasons':reasons}
            if plan['mode']=='DIAMETER_RANGE':item['outside_reference_range']=not lo<=equiv<=hi
            if reasons:rejected.append(item)
            else:
                rank=abs(float(np.log(equiv/scale))) if scale is not None else 0.
                pending.append((rank,scale_index,int(ident),yy,xx,item))
    if plan['mode']=='DIAMETER_RANGE':pending.sort(key=lambda v:v[:3])
    for _,_,_,yy,xx,item in pending:
        overlaps,counts=np.unique(labels[yy,xx],return_counts=True)
        hits=[(int(i),int(n)) for i,n in zip(overlaps,counts) if i]
        if hits:
            duplicate=False
            if len(hits)==1:
                kept,intersection=hits[0];other=objects[kept-1]['area_px']
                iou=intersection/(item['area_px']+other-intersection)
                contained=intersection/min(item['area_px'],other)
                duplicate=iou>=.5 or contained>=.8
            item['reasons']=['duplicate_across_scales' if duplicate else 'cross_scale_overlap_requires_review']
            item['overlaps_kept_ids']=[i for i,_ in hits];rejected.append(item);continue
        ident=len(objects)+1;labels[yy,xx]=ident;objects.append(dict(item,id=ident))
    return labels,objects,rejected
