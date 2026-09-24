"""One-way volume/area-constrained capsule geometry; no membrane mechanics.

Radius feasibility need not be monotone in curved walls. A deterministic
largest-radius-first branch-and-bound covers the entire area-feasible interval.
Pruning uses a Hausdorff Lipschitz bound, never a coarse grid or tube-radius guess.
"""
from dataclasses import dataclass
import heapq
import numpy as np
from .particle_shapes import Capsule,Ellipsoid,EPS,roundoff_length,unit
from .wall_gap import wall_gap
from .convex_triangle import triangle_gap,triangle_closest_many,capsule_triangle_many,TriangleGap,feature

AREA_BUDGET_ROLE='MODEL_DERIVED_FROM_PARTICLE2_OBLATE_GEOMETRY'
SEARCH_NODE_LIMIT=20000  # implementation guard, not a geometric/physical parameter


def oblate_area(a_m,c_m):
    a=float(a_m);c=float(c_m)
    if not (np.isfinite([a,c]).all() and 0<c<=a):raise ValueError('oblate a>=c>0 required')
    if a==c:return float(4*np.pi*a**2)
    e2=max(0.,1-(c/a)**2)
    if e2<np.sqrt(EPS):
        # atanh(e)/e=1+e²/3+e⁴/5+..., including the exact sphere limit.
        ratio=1+e2/3+e2**2/5+e2**3/7
    else:
        e=np.sqrt(e2);ratio=np.arctanh(e)/e
    return float(2*np.pi*a*a*(1+(1-e2)*ratio))


def capsule_length(volume_m3,radius_m):
    length=volume_m3/(np.pi*radius_m**2)-4*radius_m/3
    if length<0 and abs(length)<=roundoff_length(radius_m):return 0.
    if length<0:raise ValueError('capsule radius exceeds equal-volume sphere')
    return float(length)


def capsule_area(volume_m3,radius_m):
    return float(2*volume_m3/radius_m+(4/3)*np.pi*radius_m**2)


def area_feasible_interval(geometry):
    volume=float(geometry.volume_m3);budget=oblate_area(geometry.a_m,geometry.c_m)
    upper=float(np.cbrt(3*volume/(4*np.pi)))
    low=upper/2
    while capsule_area(volume,low)<=budget:low/=2
    high=upper;tol=roundoff_length(upper)
    for _ in range(64):
        if high-low<=tol:break
        middle=(low+high)/2
        if capsule_area(volume,middle)<=budget:high=middle
        else:low=middle
    return high,upper,budget


@dataclass(frozen=True)
class DeformationDecision:
    status: str
    shape: object | None
    gap: object | None
    original_oblate_gap_m: float
    area_budget_m2: float
    area_budget_role: str
    radius_interval_m: tuple
    search_evaluations: int
    radius_optimality_bound_m: float
    reason: str


@dataclass(frozen=True)
class _FeasibilityBound:
    """Search-only proof; replaced by an exact full-WALL gap before acceptance."""
    gap_m: float = 0.


def select_rbc_shape(center,geometry,q,free_velocity,wall,*,inside_lumen=True,velocity_scale=None,node_limit=SEARCH_NODE_LIMIT,use_radius_cache=True):
    oblate=Ellipsoid.from_rbc(center,geometry,q)
    # A negative actual triangle gap certifies a conflict without finding the
    # most deeply overlapped triangle. FREE acceptance still queries the full WALL.
    nearest,center_distance=wall.nearest_center_triangle(center)
    quick=triangle_gap(oblate,wall.triangles[nearest])
    if quick.gap_m>=-quick.roundoff_m:
        # An affine map takes the ellipsoid to the unit ball. Exact triangle
        # closest-point geometry there cheaply identifies a true intersected
        # feature; then the signed Euclidean kernel certifies the conflict.
        # No point cloud, effective radius or altered accepted shape is used.
        ids=wall.candidates(center,oblate.bounding_radius_m)
        if len(ids):
            mapped=((wall.triangles[ids]-np.asarray(center))@oblate.rotation)/oblate.axes_m
            points,_=triangle_closest_many(np.zeros(3),mapped)
            distances=np.linalg.norm(points,axis=1);k=int(np.argmin(distances))
            if distances[k]<1-512*EPS:quick=triangle_gap(oblate,wall.triangles[int(ids[k])])
    original=quick if quick.gap_m < -quick.roundoff_m else wall_gap(oblate,wall,inside_lumen=inside_lumen)
    lower,upper,budget=area_feasible_interval(geometry)
    area_interval=(lower,upper)
    def decision(status,shape=None,gap=None,count=0,bound=0.,reason=''):
        if isinstance(shape,Capsule) and use_radius_cache:
            # Feasibility bounds accelerate only the search. Accepted states
            # always get the exact full-WALL minimum and its real witness.
            gap=wall_gap(shape,wall,inside_lumen=True)
            if gap.gap_m < -geometry_tol:
                return DeformationDecision('DEFORMATION_SEARCH_UNRESOLVED',None,None,original.gap_m,budget,AREA_BUDGET_ROLE,area_interval,count,bound,'Cached feasibility certificate disagrees with final exact WALL query')
        return DeformationDecision(status,shape,gap,original.gap_m,budget,AREA_BUDGET_ROLE,area_interval,count,bound,reason)
    if original.gap_m>=-original.roundoff_m:return decision('FREE_OBLATE',oblate,original)
    if not inside_lumen:return decision('CENTER_OUTSIDE_LUMEN',reason='P0 center membership fails; not redefined as wall gap')
    velocity=np.asarray(free_velocity,dtype=float);speed=np.linalg.norm(velocity)
    scale=speed if velocity_scale is None else float(velocity_scale)
    if not np.isfinite(speed) or speed<=256*EPS*max(scale,np.finfo(float).tiny):
        return decision('DEFORMATION_AXIS_UNRESOLVED',reason='Flow direction is below float64 resolution relative to supplied velocity scale')
    axis=unit(velocity);cache={};best=None;active_triangle=None
    # Reserve the coordinate-scale WALL budget for independent witnesses and
    # continuous-sweep checks. Radius search uses the stricter shape-scale
    # roundoff only; it must not spend the entire WALL tolerance on a negative gap.
    geometry_tol=roundoff_length(upper)
    radius_tol=roundoff_length(upper)
    history_gap=np.full(len(wall.triangles),np.nan);history_r=np.zeros(len(wall.triangles));history_h=np.zeros(len(wall.triangles))
    def cached_feasibility(cap):
        # Signed convex distance is 1-Lipschitz under a support/Hausdorff
        # perturbation. For fixed center/axis, |dR|+|d(L/2)| is such a bound.
        # Thus previous exact triangle gaps can safely exclude most distant
        # features during radius refinement, without assuming radius monotonicity.
        ids=wall.candidates(center,cap.bounding_radius_m);half=cap.cylindrical_length_m/2
        delta=wall.triangle_centers[ids]-np.asarray(center);z=np.clip(delta@axis,-half,half)
        lower_bound=np.linalg.norm(delta-z[:,None]*axis,axis=1)-cap.radius_m-wall.triangle_bounding_radii[ids]-wall.roundoff_m
        known=np.isfinite(history_gap[ids])
        lower_bound[known]=np.maximum(lower_bound[known],history_gap[ids[known]]-abs(cap.radius_m-history_r[ids[known]])-abs(half-history_h[ids[known]])-wall.roundoff_m)
        needed=ids[lower_bound < -geometry_tol]
        # Chunking permits an exact negative witness to stop the search query.
        # A feasible result is returned only after EVERY candidate is certified.
        for begin in range(0,len(needed),256):
            batch=needed[begin:begin+256];gs,wp,pp,ns,bary=capsule_triangle_many(cap,wall.triangles[batch])
            history_gap[batch]=gs;history_r[batch]=cap.radius_m;history_h[batch]=half
            k=int(np.argmin(gs))
            if gs[k] < -geometry_tol:
                return TriangleGap(float(gs[k]),wp[k],pp[k],ns[k],feature(bary[k]),wall.roundoff_m)
        return _FeasibilityBound()
    def evaluate(radius):
        nonlocal active_triangle
        key=float(radius)
        if key not in cache:
            cap=Capsule(center,axis,key,capsule_length(geometry.volume_m3,key))
            if use_radius_cache:
                cache[key]=(cap,cached_feasibility(cap))
                return cache[key]
            # One truly intersected triangle already certifies infeasibility.
            # Reuse that constraint to prune radius bands, but a candidate is
            # declared feasible ONLY after a complete WALL distance query.
            local=None if active_triangle is None else triangle_gap(cap,wall.triangles[active_triangle])
            if local is not None and local.gap_m < -geometry_tol:
                result=local
            else:
                result=wall_gap(cap,wall,inside_lumen=True)
                active_triangle=result.wall_triangle_id
            cache[key]=(cap,result)
        return cache[key]
    # Every capsule contains a sphere of radius R at its center. This necessary
    # WALL bound cannot skip a feasible radius, even in a curved vessel.
    if center_distance<lower-geometry_tol:
        return decision('DEFORMATION_SURROGATE_INFEASIBLE',reason='Even the smallest area-feasible center sphere intersects WALL')
    upper=min(upper,center_distance)
    cap,gap=evaluate(upper)
    if gap.gap_m>=-geometry_tol:
        return decision('CAPILLARY_DEFORMED',cap,gap,len(cache),geometry_tol,'Largest radius reaches the certified volume/center-sphere upper bound')
    # |dR| + |d(L/2)| bounds support displacement for fixed center/axis.
    # L decreases with R, so a whole interval has bound 5/3+V/(pi*R_low³).
    queue=[(-upper,lower)];unresolved=[]
    while queue:
        negative_hi,lo=heapq.heappop(queue);hi=-negative_hi
        if best is not None and hi<=best[0]+radius_tol:continue
        if len(cache)>=node_limit:
            return decision('DEFORMATION_SEARCH_UNRESOLVED',count=len(cache),reason='Deterministic branch-and-bound evaluation guard reached')
        top_shape,top_gap=evaluate(hi)
        if top_gap.gap_m>=-geometry_tol:
            best=(hi,top_shape,top_gap)
            continue
        # For R in [lo,hi], L(R)>=L(hi). Thus Capsule(lo,L(hi)) is contained
        # in every candidate shape. Shrinking R at fixed shaft raises its signed
        # gap by exactly hi-lo. This excludes an upper radius band without any
        # assumption that the original wall-feasibility function is monotone.
        certified_hi=hi+top_gap.gap_m+geometry_tol
        if certified_hi<lo:continue
        if certified_hi<hi-radius_tol:
            heapq.heappush(queue,(-certified_hi,lo))
            continue
        mid=.5*(lo+hi);cap,gap=evaluate(mid)
        if gap.gap_m>=-geometry_tol and (best is None or mid>best[0]):best=(mid,cap,gap)
        lipschitz=5/3+geometry.volume_m3/(np.pi*lo**3)
        if gap.gap_m+lipschitz*(hi-lo)/2< -geometry_tol:continue
        # A small interval can still lie more than radius_tol ABOVE the best
        # feasible point. Do not abandon it merely because its own width is
        # small: continue excluding/refining until the optimality bound closes.
        # Only a true representable-radius floor leaves an unresolved interval.
        if mid<=lo or mid>=hi:
            for r in [lo,hi]:
                candidate,g=evaluate(r)
                if g.gap_m>=-geometry_tol and (best is None or r>best[0]):best=(r,candidate,g)
            if best is None or hi>best[0]+radius_tol:unresolved.append((lo,hi))
            continue
        heapq.heappush(queue,(-hi,mid));heapq.heappush(queue,(-mid,lo))
    if unresolved and (best is None or max(hi for _,hi in unresolved)>best[0]+radius_tol):
        return decision('DEFORMATION_SEARCH_UNRESOLVED',count=len(cache),reason=f'Potential feasible interval remains unresolved at geometry roundoff; best_radius={None if best is None else best[0]!r}; intervals={unresolved!r}; radius_tol={radius_tol!r}')
    if best is None:
        return decision('DEFORMATION_SURROGATE_INFEASIBLE',count=len(cache),reason='Entire area-feasible radius interval excluded by certified gap bounds')
    return decision('CAPILLARY_DEFORMED',best[1],best[2],len(cache),radius_tol,'Largest wall-feasible radius certified to geometry roundoff over the full area-feasible interval')
