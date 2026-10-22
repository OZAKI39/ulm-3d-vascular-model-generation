"""Deterministic rigid orientation search; the casting frame is immutable."""
from __future__ import annotations
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from .support_detection import section_polygon,contiguous_runs

def rotations(count=128):
    if count<80 or count>160:raise ValueError('COARSE_COUNT_MUST_BE_80_TO_160')
    # Principal directions guarantee that broad flat box faces can meet the bed.
    directions=[np.array(v,float) for v in [(0,0,1),(0,0,-1),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0)]]
    n=count//4-len(directions);gold=np.pi*(3-np.sqrt(5))
    for i in range(n):
        z=1-2*(i+.5)/n;r=np.sqrt(1-z*z)
        directions.append(np.array([r*np.cos(i*gold),r*np.sin(i*gold),z]))
    result=[]
    for direction in directions:
        alignment=trimesh.geometry.align_vectors(direction,[0,0,1])[:3,:3]
        for roll in [0,90,180,270]:
            result.append(Rotation.from_euler('z',roll,degrees=True).as_matrix()@alignment)
    return result

def bed_transform(vertices,rotation,build_bounds):
    points=np.asarray(vertices)@rotation.T;lo=points.min(0);hi=points.max(0)
    bounds=np.asarray(build_bounds);matrix=np.eye(4);matrix[:3,:3]=rotation
    matrix[:2,3]=(bounds[0,:2]+bounds[1,:2]-lo[:2]-hi[:2])/2
    matrix[2,3]=-lo[2]
    return matrix

def rigid(matrix):
    m=np.asarray(matrix);return bool(np.allclose(m[3],[0,0,0,1]) and np.allclose(m[:3,:3].T@m[:3,:3],np.eye(3),atol=1e-8) and np.isclose(np.linalg.det(m[:3,:3]),1.))

def fits(vertices,matrix,build_bounds,margin=2.):
    q=np.asarray(vertices)@matrix[:3,:3].T+matrix[:3,3];bounds=np.asarray(build_bounds)
    return bool(np.all(q.min(0)>=bounds[0]+[margin,margin,-1e-4]) and np.all(q.max(0)<=bounds[1]-[margin,margin,0]))

def coarse(mold,vascular,box,branches,profile,count=128):
    rows=[];bounds=np.asarray(profile['build_bounds_mm']);points=vascular.triangles_center
    branch_points=np.vstack([b['xyz'] for b in branches]);port_flags=np.concatenate([np.full(len(b['xyz']),b['is_port']) for b in branches])
    from scipy.spatial import cKDTree
    is_port=port_flags[cKDTree(branch_points).query(points)[1]]
    for i,r in enumerate(rotations(count)):
        m=bed_transform(mold.vertices,r,bounds);normal=vascular.face_normals@r.T;down=normal[:,2]<-np.sqrt(.5)
        bnormal=box.face_normals@r.T;bc=box.triangles_center@r.T+m[:3,3]
        boxarea=float(box.area_faces[(bnormal[:,2]<-np.sqrt(.5))&(bc[:,2]>.2)].sum())
        vc=points@r.T+m[:3,3];area=float(vascular.area_faces[down&(vc[:,2]>.2)].sum())
        portarea=float(vascular.area_faces[down&is_port&(vc[:,2]>.2)].sum())
        critical=0.;longest=0.
        for b in branches:
            t=np.asarray(b['tangent'])@r.T;mask=np.abs(t[:,2])<.5
            for a,z in contiguous_runs(mask):
                length=float(b['arc'][z]-b['arc'][a]);critical+=length;longest=max(longest,length)
        transformed=mold.copy();transformed.apply_transform(m);lo,hi=transformed.bounds
        # Six real coarse sections, explicitly a surrogate rather than the
        # 0.10 mm production detector used for every retained orientation.
        levels=np.linspace(lo[2]+.15,hi[2]-.15,6);islands=0
        for path in transformed.section_multiplane([0,0,0],[0,0,1],levels):
            poly=section_polygon(path);islands+=len(poly.geoms) if hasattr(poly,'geoms') else int(not poly.is_empty)
        ok=fits(mold.vertices,m,bounds)
        # A visible surface estimate only. Exact roots are compared later.
        reach=float(np.mean(vc[down,2])) if down.any() else 0.
        row=dict(orientation_id=f'O{i:03}',rotation=r.tolist(),print_transform=m.tolist(),casting_restore_transform=np.linalg.inv(m).tolist(),
            gravity_in_casting=r[2].tolist(),vascular_unsupported_area_mm2=area,critical_interval_length_mm=critical,
            long_unsupported_span_mm=longest,coarse_layer_islands=islands,box_overhang_burden_mm2=boxarea,port_overhang_burden_mm2=portarea,
            estimated_support_reach_mm=reach,z_height_mm=float(hi[2]-lo[2]),footprint_mm2=float(np.prod((hi-lo)[:2])),build_volume_fit=ok,
            rejection_reason='' if ok else 'BUILD_VOLUME_EXCEEDED',surrogate_not_full_fdm=True)
        row['coarse_rank_value']=area+boxarea+2*critical+.5*islands+.05*reach
        rows.append(row)
    return rows

def diverse_top(rows,n=12):
    selected=[]
    for row in sorted(rows,key=lambda r:(not r['build_volume_fit'],round(r['coarse_rank_value'],6),r['orientation_id'])):
        if not row['build_volume_fit']:continue
        direction=np.array(row['gravity_in_casting'])
        if any(np.dot(direction,r['gravity_in_casting'])>np.cos(np.radians(5)) for r in selected):continue
        selected.append(row)
        if len(selected)==n:break
    return selected
