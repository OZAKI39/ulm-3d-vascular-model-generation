"""Exact spheres and conservative swept-ball broadphase (no GJK, no skins)."""
import numpy as np
from scipy.spatial import cKDTree
from ..particle_shapes import Sphere, roundoff_length
from ..pair_geometry import ParticlePairGapResult
from ..nearfield_regularization import contract as inherited_contract


def sphere_gap(a, b, particle_i_id=0, particle_j_id=1):
    if type(a) is not Sphere or type(b) is not Sphere:
        raise TypeError('MB_MB_INTERACTION_V1 supports exact spheres only')
    i, j = int(particle_i_id), int(particle_j_id)
    if i == j: raise ValueError('Distinct stable IDs required')
    if i > j:
        r = sphere_gap(b, a, j, i)
        return ParticlePairGapResult(i,j,a.mode,b.mode,r.gap_m,r.state,r.point_j_m,
            r.point_i_m,-r.normal_j_to_i,r.distance_m,0,r.roundoff_budget_m,(j,i))
    delta=a.center_m-b.center_m; length=float(np.linalg.norm(delta))
    # A coincident-center state is invalid, but has an unambiguous overlap flag.
    n=delta/length if length else np.array([1.,0.,0.])
    gap=float(length-a.radius_m-b.radius_m)
    pad=roundoff_length(a.center_m,b.center_m,a.radius_m,b.radius_m)
    state='SEPARATED' if gap>pad else 'PENETRATING' if gap < -pad else 'TOUCHING'
    return ParticlePairGapResult(i,j,a.mode,b.mode,gap,state,a.center_m-a.radius_m*n,
        b.center_m+b.radius_m*n,n,abs(gap),0,pad,(i,j))


def candidates(shapes, *, velocities=None, dt=0.):
    """Enclose each entire straight path, including the continuum handoff shell.

    Broadphase overestimation is harmless. Exact geometry chooses activation.
    Bounding spheres include half path length; therefore crossing paths cannot
    disappear when both endpoints are clear. No historical neighbor list reused.
    """
    ids=sorted(shapes)
    if len(ids)<2: return []
    x=np.array([shapes[i].center_m for i in ids]);r=np.array([shapes[i].radius_m for i in ids])
    delta=np.zeros_like(x) if velocities is None else dt*np.array([velocities[i] for i in ids])
    center=x+.5*delta
    pad=roundoff_length(x,x+delta,r)
    bound=r+.5*np.linalg.norm(delta,axis=1)+np.maximum(.05*r,2e-9)+pad
    pairs=cKDTree(center).query_pairs(2*float(bound.max()),output_type='ndarray')
    if not len(pairs):return []
    d=center[pairs[:,0]]-center[pairs[:,1]]
    keep=np.linalg.norm(d,axis=1)<=bound[pairs[:,0]]+bound[pairs[:,1]]
    return sorted((ids[a],ids[b]) for a,b in pairs[keep])


def components(ids, pairs):
    parent={int(i):int(i) for i in ids}
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in pairs:
        a,b=find(i),find(j)
        if a!=b:parent[max(a,b)]=min(a,b)
    groups={}
    for i in sorted(parent):groups.setdefault(find(i),[]).append(i)
    return [tuple(v) for _,v in sorted(groups.items())]


def contract():
    return dict(contract_name='MB_MB_INTERACTION_V1',version=1,
        parameter_origin=inherited_contract(),
        formalization='SPHERICAL_NORMAL_PAIR_MODEL_ONLY',
        broadphase='FRESH_CONSERVATIVE_SWEPT_BALL_KDTREE_CPU_OR_GPU_CELL_LIST',
        neighbor_skin_m=0.,gap='norm(x_i-x_j)-a_i-a_j',
        narrow_phase='EXACT_ANALYTICAL_SPHERE_NO_GJK',
        temporal_certificate='EXACT_MINIMUM_GAP_ON_HELD_VELOCITY_STRAIGHT_PATH',
        coupled_contact='R_METRIC_SIMULTANEOUS_NONPENETRATION',
        nominal_dt_s=.0005,adaptive_subdivision='INHERITED_P3_REAL_PHYSICAL_TIME',
        no_penetration_budget='512*float64_eps*max_absolute_coordinate_or_radius',
        contact_meaning='CONTINUUM_HANDOFF_NOT_MOLECULAR_CONTACT',
        scope_limit='NO_PAIR_TANGENTIAL_LUBRICATION_DEFORMATION_SHELL_ACOUSTICS_OR_TWO_WAY_CFD')
