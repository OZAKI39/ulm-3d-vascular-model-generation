"""Shared PIV arrow geometry for live display and recorded overlays.

Only display sampling/scaling changes; returned glyphs never modify velocities.
"""
import numpy as np

GAIN_SECONDS = .8
MAX_LENGTH_PX = 48.
SPACING_PX = 48.
REFERENCE_SPEED = 25.
ARROW_RGB = (255, 226, 59)
OUTLINE_RGB = (18, 22, 28)


def arrow_shape(start, velocity, gain=GAIN_SECONDS):
    start = np.asarray(start, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    speed = float(np.linalg.norm(velocity))
    if not np.isfinite(speed) or speed < 1e-9:return None
    length = min(speed*gain, MAX_LENGTH_PX)
    if length < 1.:return None  # Do not inflate near-zero motion into a direction.
    direction = velocity/speed
    end = start + direction*length
    normal = np.array([-direction[1], direction[0]])
    head_length = min(8., length*.45)
    half_width = head_length*.55
    head = np.array([end, end-direction*head_length+normal*half_width,
                     end-direction*head_length-normal*half_width])
    return start,end,head


def piv_arrows(vectors, offset=(0.,0.)):
    vectors = np.asarray(vectors, dtype=float).reshape(-1,4)
    vectors = vectors[np.isfinite(vectors).all(axis=1)]
    if not len(vectors):return []
    cells = np.floor(vectors[:,:2]/SPACING_PX).astype(np.int64)
    distance = np.sum((vectors[:,:2]-(cells+.5)*SPACING_PX)**2,axis=1)
    # Select a measured vector, never a cell average. Stable ties retain the
    # first input vector, and cell order remains the original acquisition order.
    order = np.lexsort((np.arange(len(vectors)),distance,cells[:,1],cells[:,0]))
    first = np.r_[True,np.any(np.diff(cells[order],axis=0),axis=1)]
    starts = np.flatnonzero(first)
    chosen = order[starts]
    cell_order = np.argsort(np.minimum.reduceat(order,starts))
    rows = vectors[chosen[cell_order]]
    speed = np.linalg.norm(rows[:,2:],axis=1)
    keep = speed*GAIN_SECONDS >= 1.
    rows,speed = rows[keep],speed[keep]
    if not len(rows):return []
    length = np.minimum(speed*GAIN_SECONDS,MAX_LENGTH_PX)
    direction = rows[:,2:]/speed[:,None]
    begin = rows[:,:2]+offset
    end = begin+direction*length[:,None]
    normal = np.column_stack((-direction[:,1],direction[:,0]))
    head_length = np.minimum(8.,length*.45)[:,None]
    base = end-direction*head_length
    head = np.stack((end,base+normal*head_length*.55,base-normal*head_length*.55),axis=1)
    return list(zip(begin,end,head))


def draw_cv_arrows(image, arrows):
    import cv2
    if not arrows:return
    shafts = np.rint([[start,end] for start,end,_ in arrows]).astype(np.int32)
    heads = np.rint([head for _,_,head in arrows]).astype(np.int32)
    cv2.polylines(image,shafts,False,OUTLINE_RGB[::-1],5,cv2.LINE_AA)
    cv2.polylines(image,heads,True,OUTLINE_RGB[::-1],3,cv2.LINE_AA)
    cv2.polylines(image,shafts,False,ARROW_RGB[::-1],2,cv2.LINE_AA)
    # Separate convex fills avoid even/odd holes where arrowheads overlap.
    for triangle in heads:cv2.fillConvexPoly(image,triangle,ARROW_RGB[::-1],cv2.LINE_AA)


DISPLAY_SPEC=dict(velocity_units='px/s',gain_seconds=GAIN_SECONDS,
    maximum_arrow_length_px=MAX_LENGTH_PX,sampling_cell_px=SPACING_PX,
    zero_motion='No directional arrow',reference_speed_px_s=REFERENCE_SPEED)
