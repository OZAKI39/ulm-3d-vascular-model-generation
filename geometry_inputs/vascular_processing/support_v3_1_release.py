"""Fixed-plate supports versus a single rigidly lifted mold; bounded FCL sweeps."""
from __future__ import annotations
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from .support_removal import swept_path
from .support_inputs import transform


def release_poses(mold,support,axis=None,angle=0.):
    # Local moving coordinates reduce the Lipschitz rotational travel bound.
    pivot=mold.bounds.mean(0); local=mold.copy();local.apply_translation(-pivot)
    start=np.eye(4);start[:3,3]=pivot
    rise=float(support.bounds[1,2]-mold.bounds[0,2]+5.)
    if axis is None:
        end=start.copy();end[2,3]+=rise
        return local,[start,end]
    # Track the lowest mold vertex while tilting. Rotating about the centre
    # after a tiny lift alone can put a box corner THROUGH the build plate.
    # <=1-degree keyframes plus 0.05-mm clearance bound that curved motion.
    lifted=start.copy();lifted[2,3]+=.05
    poses=[start,lifted]
    for theta in np.linspace(0,angle,int(abs(angle))+1)[1:]:
        tilted=start.copy();tilted[:3,:3]=Rotation.from_euler(axis,theta,degrees=True).as_matrix()
        minimum=transform(local.vertices,tilted)[:,2].min()
        tilted[2,3]+=max(0.,.05-minimum);poses.append(tilted)
    end=tilted.copy();end[2,3]+=rise+float(np.ptp(local.vertices[:,2]))
    return local,poses+[end]


def bed_clearance(local,poses):
    radius=float(np.linalg.norm(local.vertices,axis=1).max());lower=np.inf
    for a,b in zip(poses[:-1],poses[1:]):
        theta=Rotation.from_matrix(a[:3,:3].T@b[:3,:3]).magnitude()
        # Arc-versus-chord deviation bound applies to every moving vertex.
        endpoint_min=min(transform(local.vertices,a)[:,2].min(),transform(local.vertices,b)[:,2].min())
        lower=min(lower,float(endpoint_min-radius*(1-np.cos(theta/2))))
    return dict(passed=lower>=-1e-5,minimum_swept_lower_bound_mm=lower,method='Rigid rotational arc/chord vertex bound including build-plate plane')


def plate_release(mold, support, inner_casting=None, print_transform=None):
    if not len(support.faces):return dict(passed=True,status='NO_PLATE_ROOTED_SUPPORT',attempts=[],level=0)
    manager=trimesh.collision.CollisionManager();manager.add_object('fixed_plate_support',support)
    attempts=[]
    choices=[(None,0.)]+[(a,float(s*d)) for d in (5,10,15) for a in 'xy' for s in (1,-1)]
    for axis,angle in choices:
        local,poses=release_poses(mold,support,axis,angle)
        result=swept_path(local,[],poses,static_manager=manager)
        bed=bed_clearance(local,poses);result['buildplate_clearance']=bed
        if not bed['passed']:result.update(passed=False,status='MOLD_CROSSES_BUILD_PLATE')
        end=transform(local.vertices,poses[-1]); full_exit=bool(end[:,2].min()>support.bounds[1,2])
        result.update(axis=axis,signed_angle_deg=angle,level=1 if axis is None else 2,full_mold_exit=full_exit,
                      transforms=[p.tolist() for p in poses],moving_local_pivot=mold.bounds.mean(0).tolist())
        # For a closed intact five-wall mold, collision-free full separation
        # of a cavity-intersecting fixed component must use its sole opening.
        result['opening_exit_basis']='Intact five-wall frozen box, continuous collision-free path, complete final Z separation'
        attempts.append(result)
        if result['passed'] and full_exit:
            return dict(passed=True,status='PLATE_RELEASE_PASS',level=result['level'],path=result,attempts=attempts,
                        physical_release_validated=False)
    return dict(passed=False,status='PLATE_RELEASE_FAILED',level=None,attempts=attempts,physical_release_validated=False)
