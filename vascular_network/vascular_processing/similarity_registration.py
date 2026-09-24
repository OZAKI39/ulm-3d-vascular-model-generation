"""Deterministic, proper similarity ICP; transforms only the donor."""
import numpy as np
from scipy.spatial.transform import Rotation
from .topbrain_qc import TopBrainError

REGISTRATION_RULES=dict(twists_degrees=[0,90,180,270],minimum_fitness=.5,maximum_normalized_rmse=.12,
                        relative_scale_range=[.5,2.],maximum_root_error_fraction=.25,
                        initialization='ICA/MCA root translation + radial 90th percentile size + proximal direction')


def transform_points(points, matrix):
    return np.asarray(points)@matrix[:3,:3].T+matrix[:3,3]


def proper_similarity(matrix):
    matrix=np.asarray(matrix)
    if matrix.shape!=(4,4) or not np.isfinite(matrix).all():
        raise TopBrainError('INVALID_REGISTRATION','Nonfinite or malformed ICP transform')
    linear=matrix[:3,:3]
    if np.linalg.det(linear)<=1e-12:
        raise TopBrainError('INVALID_REGISTRATION','Reflection or collapsed ICP transform')
    s=np.linalg.svd(linear,compute_uv=False)
    if np.linalg.det(linear)<=0 or not np.allclose(s,s.mean(),rtol=1e-5):
        raise TopBrainError('INVALID_REGISTRATION','Reflection or nonuniform scale')
    return float(s.mean())


def align_direction(source,target):
    a=np.asarray(source)/np.linalg.norm(source);b=np.asarray(target)/np.linalg.norm(target)
    cross=np.cross(a,b);dot=np.clip(a@b,-1,1)
    if np.linalg.norm(cross)<1e-10:
        if dot>0:return np.eye(3)
        axis=np.cross(a,np.eye(3)[np.argmin(abs(a))]);axis/=np.linalg.norm(axis)
        return Rotation.from_rotvec(np.pi*axis).as_matrix()
    return Rotation.from_rotvec(np.arccos(dot)*cross/np.linalg.norm(cross)).as_matrix()


def register(source,target):
    from .transfer_runtime import require_versions
    require_versions()
    import open3d as o3d
    if source.side!=target.side:raise TopBrainError('SIDE_MISMATCH','No contralateral reflection')
    size=lambda g:float(np.percentile(np.linalg.norm(g.points-g.root,axis=1),90))
    a,b=size(source),size(target)
    if min(a,b)<1e-8:raise TopBrainError('INVALID_REGISTRATION','Degenerate geometry')
    initial_scale=b/a;rotation=align_direction(source.direction,target.direction)
    def cloud(points):return o3d.geometry.PointCloud(o3d.utility.Vector3dVector(points))
    src,dst=cloud(source.points),cloud(target.points);candidates=[]
    for twist in REGISTRATION_RULES['twists_degrees']:
        r=Rotation.from_rotvec(np.deg2rad(twist)*target.direction).as_matrix()@rotation
        initial=np.eye(4);initial[:3,:3]=initial_scale*r;initial[:3,3]=target.root-initial[:3,:3]@source.root
        result=None
        for threshold in [.2,.12,.08]:
            result=o3d.pipelines.registration.registration_icp(src,dst,threshold*b,initial,
                o3d.pipelines.registration.TransformationEstimationPointToPoint(with_scaling=True),
                o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=60))
            initial=np.asarray(result.transformation).copy()
            if not np.isfinite(initial).all() or np.linalg.det(initial[:3,:3])<=1e-12:break
        try:scale=proper_similarity(initial)
        except TopBrainError:continue
        root_error=float(np.linalg.norm(transform_points(source.root[None],initial)[0]-target.root))
        valid=(result.fitness>=.5 and result.inlier_rmse/b<=.12 and .5<=scale/initial_scale<=2 and root_error/b<=.25)
        candidates.append(dict(transform=initial,scale=scale,fitness=float(result.fitness),inlier_rmse=float(result.inlier_rmse),
                               normalized_rmse=float(result.inlier_rmse/b),root_error=root_error,
                               initial_scale=initial_scale,twist_degrees=twist,status='VALID' if valid else 'INVALID_REGISTRATION'))
    valid=[c for c in candidates if c['status']=='VALID']
    result=sorted(valid or candidates,key=lambda c:(-c['fitness'],c['normalized_rmse'],c['twist_degrees']))[0] if candidates else dict(status='INVALID_REGISTRATION')
    return dict(**result,source_case=source.case_id,target_case=target.case_id,side=source.side,rules=REGISTRATION_RULES)
