"""Finite convex surface to actual WALL gap; center lumen membership stays in P0."""
from dataclasses import dataclass
import numpy as np
from .convex_triangle import triangle_gap,triangle_closest_many,capsule_triangle_many,TriangleGap,feature
from .particle_shapes import roundoff_length,Sphere,Capsule


@dataclass(frozen=True)
class WallGapResult:
    gap_m: float
    state: str
    inside_lumen: bool | None
    wall_triangle_id: int
    wall_feature: str
    wall_point_m: np.ndarray
    particle_point_m: np.ndarray
    normal_inward: np.ndarray
    distance_m: float
    shape_mode: str
    candidate_count: int
    roundoff_m: float


def wall_gap(shape,wall,*,inside_lumen=None):
    nearest,_=wall.nearest_center_triangle(shape.center_m)
    first=triangle_gap(shape,wall.triangles[nearest]);rho=shape.bounding_radius_m
    # A certified current clearance upper bound expands a distance-query search,
    # not a physical coating: omitted triangles obey d(center,T)-rho>best_gap.
    radius=rho+max(0.,first.gap_m)
    ids=wall.candidates(shape.center_m,radius)
    if nearest not in ids:ids=np.sort(np.r_[ids,nearest])
    center_points,_=triangle_closest_many(shape.center_m,wall.triangles[ids])
    lower=np.linalg.norm(center_points-shape.center_m,axis=1)-rho
    best=first;best_id=nearest;pad=roundoff_length(shape.center_m,rho,wall.triangles)
    if isinstance(shape,(Sphere,Capsule)):
        capsule=shape if isinstance(shape,Capsule) else Capsule(shape.center_m,[0,0,1],shape.radius_m,0.)
        # Each finite triangle lies in its centroid ball. Enlarging that triangle
        # to the ball can only LOWER its signed convex distance. Capsule/ball
        # distance is point-to-axis-segment distance minus both radii, including
        # overlaps. This certified lower bound rejects distant BVH candidates
        # before the exact prism query; final gaps still use finite triangles.
        delta=wall.triangle_centers[ids]-capsule.center_m
        along=np.clip(delta@capsule.axis_world,-capsule.cylindrical_length_m/2,capsule.cylindrical_length_m/2)
        ball_lower=np.linalg.norm(delta-along[:,None]*capsule.axis_world,axis=1)-capsule.radius_m-wall.triangle_bounding_radii[ids]
        exact_ids=ids[ball_lower<=first.gap_m+pad]
        if nearest not in exact_ids:exact_ids=np.sort(np.r_[exact_ids,nearest])
        gaps,wp,pp,ns,bary=capsule_triangle_many(capsule,wall.triangles[exact_ids])
        k=int(np.argmin(gaps));best_id=int(exact_ids[k]);best=TriangleGap(float(gaps[k]),wp[k],pp[k],ns[k],feature(bary[k]),pad)
        lower=np.full(len(ids),np.inf)  # all candidates already evaluated by the same finite-feature mathematics
    for k in np.argsort(lower,kind='stable'):
        i=int(ids[k])
        if lower[k]>best.gap_m+pad:break
        if i==nearest:continue
        result=triangle_gap(shape,wall.triangles[i])
        if result.gap_m<best.gap_m or (result.gap_m==best.gap_m and i<best_id):best,best_id=result,i
    g=best.gap_m
    state='SEPARATED' if g>pad else 'PENETRATING' if g< -pad else 'TOUCHING'
    return WallGapResult(g,state,inside_lumen,best_id,best.wall_feature,best.wall_point_m,best.particle_point_m,
                         best.normal_inward,abs(g),shape.mode,len(ids),pad)


def touching_contacts(shape,wall):
    result=[];pad=roundoff_length(shape.center_m,shape.bounding_radius_m,wall.triangles)
    ids=wall.candidates(shape.center_m,shape.bounding_radius_m)
    if isinstance(shape,(Sphere,Capsule)):
        capsule=shape if isinstance(shape,Capsule) else Capsule(shape.center_m,[0,0,1],shape.radius_m,0.)
        gs,wp,pp,ns,bary=capsule_triangle_many(capsule,wall.triangles[ids])
        return [(int(ids[k]),TriangleGap(float(gs[k]),wp[k],pp[k],ns[k],feature(bary[k]),pad)) for k in np.flatnonzero(np.abs(gs)<=pad)]
    for i in ids:
        g=triangle_gap(shape,wall.triangles[i])
        if abs(g.gap_m)<=pad:result.append((int(i),g))
    return result
