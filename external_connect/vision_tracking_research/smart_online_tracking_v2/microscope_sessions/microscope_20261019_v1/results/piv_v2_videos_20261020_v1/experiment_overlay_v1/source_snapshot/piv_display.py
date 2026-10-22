"""Shared PIV arrow geometry for live display and recorded overlays.

Display geometry and optional smoothing are separate from the original
observations. Neither renderer nor display filter mutates measured velocities.
"""
import numpy as np

GAIN_SECONDS = .8
MAX_LENGTH_PX = 40.
SPACING_PX = 48.
REFERENCE_SPEED = 25.
ARROW_RGB = (255, 226, 59)
OUTLINE_RGB = (18, 22, 28)
STROKE_PX = 1.35
HALO_PX = 2.5
SMOOTH_SPACING_PX = 72.
SMOOTH_TAU_S = .12


def arrow_shape(start, velocity, gain=GAIN_SECONDS):
    start = np.asarray(start, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    speed = float(np.linalg.norm(velocity))
    if not np.isfinite(speed) or speed < 1e-9:return None
    length = min(speed*gain, MAX_LENGTH_PX)
    if length < 3.:return None  # Do not turn near-zero motion into tiny arrowhead dots.
    direction = velocity/speed
    end = start + direction*length
    normal = np.array([-direction[1], direction[0]])
    head_length = min(7., length*.38)
    half_width = head_length*.50
    head = np.array([end, end-direction*head_length+normal*half_width,
                     end-direction*head_length-normal*half_width])
    return start,end,head


def piv_arrows(vectors, offset=(0.,0.), sample=True):
    vectors = np.asarray(vectors, dtype=float).reshape(-1,4)
    vectors = vectors[np.isfinite(vectors).all(axis=1)]
    if not len(vectors):return []
    rows=vectors
    if sample:
        cells = np.floor(vectors[:,:2]/SPACING_PX).astype(np.int64)
        distance = np.sum((vectors[:,:2]-(cells+.5)*SPACING_PX)**2,axis=1)
        # Raw view selects actual measurements, with stable tie breaking.
        order = np.lexsort((np.arange(len(vectors)),distance,cells[:,1],cells[:,0]))
        first = np.r_[True,np.any(np.diff(cells[order],axis=0),axis=1)]
        starts = np.flatnonzero(first)
        chosen = order[starts]
        cell_order = np.argsort(np.minimum.reduceat(order,starts))
        rows = vectors[chosen[cell_order]]
    speed = np.linalg.norm(rows[:,2:],axis=1)
    keep = speed*GAIN_SECONDS >= 3.
    rows,speed = rows[keep],speed[keep]
    if not len(rows):return []
    length = np.minimum(speed*GAIN_SECONDS,MAX_LENGTH_PX)
    direction = rows[:,2:]/speed[:,None]
    begin = rows[:,:2]+offset
    end = begin+direction*length[:,None]
    normal = np.column_stack((-direction[:,1],direction[:,0]))
    head_length = np.minimum(7.,length*.38)[:,None]
    base = end-direction*head_length
    head = np.stack((end,base+normal*head_length*.50,base-normal*head_length*.50),axis=1)
    return list(zip(begin,end,head))


def draw_cv_arrows(image, arrows):
    import cv2
    if not arrows:return
    shafts = np.rint([[start,end] for start,end,_ in arrows]).astype(np.int32)
    heads = np.rint([head[[1,0,2]] for _,_,head in arrows]).astype(np.int32)
    # Original-resolution exports are usually downscaled by the player. Keep
    # a 2 px bright core so the narrow halo cannot swallow it on resize.
    cv2.polylines(image,shafts,False,OUTLINE_RGB[::-1],3,cv2.LINE_AA)
    cv2.polylines(image,heads,False,OUTLINE_RGB[::-1],3,cv2.LINE_AA)
    cv2.polylines(image,shafts,False,ARROW_RGB[::-1],2,cv2.LINE_AA)
    cv2.polylines(image,heads,False,ARROW_RGB[::-1],2,cv2.LINE_AA)


class PIVDisplayFilter:
    """Causal, display-only coarse field. Never used by scientific metrics.

    Every returned cell has current valid measurements. Missing cells are
    removed immediately; no invented particles, extrapolation or stale hold.
    """
    def __init__(self):self.reset()

    def reset(self):
        self.context=None;self.time=None;self.previous={};self.mask=None

    def update(self,state,ttl_ms=250.):
        if not state.get('valid'):
            self.reset();return np.empty((0,4))
        t=float(state['arrival_monotonic'])
        geometry=state.get('geometry')
        context=(state['session_id'],state.get('geometry_revision'),
                 state.get('motion_reference'),repr(geometry))
        if context!=self.context or self.time is None or not 0<t-self.time<=ttl_ms/1000:
            self.reset();self.context=context
            if geometry:
                import cv2
                w,h=geometry['image_wh'];self.mask=np.zeros((h,w),np.uint8)
                cv2.fillPoly(self.mask,[np.rint(geometry['channel_polygon']).astype(np.int32)],1)
                for poly in geometry.get('excluded_polygons',[]):
                    cv2.fillPoly(self.mask,[np.rint(poly).astype(np.int32)],0)
        v=np.asarray(state.get('vectors',[]),float).reshape(-1,4)
        v=v[np.isfinite(v).all(axis=1)]
        if not len(v):self.previous={};self.time=t;return np.empty((0,4))
        cells,inverse=np.unique(np.floor(v[:,:2]/SMOOTH_SPACING_PX).astype(np.int64),axis=0,return_inverse=True)
        count=np.bincount(inverse)
        velocity=np.column_stack([np.bincount(inverse,weights=v[:,j])/count for j in [2,3]])
        anchors=(cells+.5)*SMOOTH_SPACING_PX
        if self.mask is not None:
            pos=np.rint(anchors+geometry['crop_xywh'][:2]).astype(int)
            h,w=self.mask.shape;inside=(pos[:,0]>=0)&(pos[:,0]<w)&(pos[:,1]>=0)&(pos[:,1]<h)
            inside[inside]=self.mask[pos[inside,1],pos[inside,0]]>0
            cells,anchors,velocity=cells[inside],anchors[inside],velocity[inside]
        alpha=1. if self.time is None else -np.expm1(-(t-self.time)/SMOOTH_TAU_S)
        current={}
        for i,cell in enumerate(cells):
            key=tuple(cell);old=self.previous.get(key);new=velocity[i]
            # Rest and reversals are shown immediately, not hidden by filter lag.
            if old is not None and np.linalg.norm(new)*GAIN_SECONDS>=3 and np.dot(old,new)>0:
                velocity[i]=old+alpha*(new-old)
            current[key]=velocity[i].copy()
        self.previous=current;self.time=t
        return np.column_stack((anchors,velocity))


def display_vectors(state,smooth=True):
    return state.get('piv_display_vectors',state.get('vectors',[])) if smooth else state.get('vectors',[])


def display_is_smoothed(state,layers):
    return layers.get('Smooth PIV',True) and 'piv_display_vectors' in state


DISPLAY_SPEC=dict(velocity_units='px/s',gain_seconds=GAIN_SECONDS,
    maximum_arrow_length_px=MAX_LENGTH_PX,sampling_cell_px=SPACING_PX,
    zero_motion='No directional arrow below 3 display pixels',reference_speed_px_s=REFERENCE_SPEED,
    glyph='Thin open arrowhead; no filled triangle',stroke_screen_px=STROKE_PX,
    smooth_grid_px=SMOOTH_SPACING_PX,smooth_time_constant_s=SMOOTH_TAU_S,
    smoothing='Display only: local spatial mean then causal EMA; no missing-cell persistence; raw data unchanged')
