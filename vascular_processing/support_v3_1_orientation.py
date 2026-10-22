"""Casting-local inverted orientations and explicitly coarse visibility metrics."""
from __future__ import annotations
import numpy as np
import trimesh
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree
from .support_inputs import transform
from .support_v2_orientation import bed_transform, fits, rigid
from .support_detection import section_polygon


def rotation_for(name):
    if name == 'UPRIGHT_O000':
        return np.eye(3)
    inverted = Rotation.from_euler('x', 180, degrees=True).as_matrix()
    if name == 'INV_0':
        return inverted
    axis, sign, angle = name[4].lower(), name[5], int(name[6:])
    if axis not in 'xy' or sign not in 'PM' or angle not in (15, 30, 45):
        raise ValueError('UNKNOWN_CASTING_LOCAL_ORIENTATION')
    # Local-axis rotation after inversion: (Ri Rlocal Ri^-1) Ri.
    return inverted @ Rotation.from_euler(axis, angle if sign == 'P' else -angle, degrees=True).as_matrix()


def orientation_audit(name, matrix):
    r = np.asarray(matrix)[:3, :3]
    normal = r @ [0., 0., 1.]
    angle = float(np.degrees(np.arccos(np.clip(-normal[2], -1, 1))))
    expected = 180. if name == 'UPRIGHT_O000' else 0. if name == 'INV_0' else float(name[6:])
    error = float(np.max(abs(r - rotation_for(name))))
    good = rigid(matrix) and error < 1e-10 and abs(angle-expected) < 1e-8
    return dict(opening_normal_print=normal.tolist(), bottom_normal_print=(-normal).tolist(),
                angle_to_print_down_deg=angle, expected_tilt_deg=expected, angle_error_deg=abs(angle-expected),
                matrix_label_error=error, label_matrix_matches=bool(good),
                orientation_status='ORIENTATION_VERIFIED' if good else 'ORIENTATION_LABEL_MATRIX_MISMATCH',
                axis_convention='Casting-local X/Y; R = Rx(180) @ Rlocal(signed tilt)', scale_factor=1.)


def orientations(mold, profile, upright=None):
    names = ['UPRIGHT_O000', 'INV_0'] + [f'INV_{a}{s}{d}' for a in 'XY' for d in (15,30,45) for s in 'PM']
    rows = []
    for name in names:
        m = np.asarray(upright['print_transform']) if name == 'UPRIGHT_O000' and upright else bed_transform(mold.vertices, rotation_for(name), profile['build_bounds_mm'])
        rows.append(dict(orientation_id=name, print_transform=m.tolist(), casting_restore_transform=np.linalg.inv(m).tolist(),
                         rotation=m[:3,:3].tolist(), build_volume_fit=fits(mold.vertices,m,profile['build_bounds_mm']),
                         **orientation_audit(name,m)))
    return rows


def vertical_clear(mesh, targets, clearance=.2):
    """Ray-only ranking proxy, excluding the target's intentional air gap."""
    targets = np.asarray(targets,dtype=float)
    if not len(targets): return np.zeros(0, bool)
    origins = targets.copy(); origins[:,2] -= clearance
    locations, ray, _ = mesh.ray.intersects_location(origins, np.tile([0,0,-1.], (len(origins),1)), multiple_hits=True)
    clear = origins[:,2] > 0
    if len(ray):
        blocked = ray[(locations[:,2] >= -.0001) & (locations[:,2] < origins[ray,2]-.00001)]
        clear[np.unique(blocked)] = False
    clear[mesh.contains(origins)] = False
    return clear


def cavity_visibility(model_print, inner_casting, matrix, samples=15):
    """Plate-to-cavity rays must intersect the REAL casting opening rectangle."""
    lo, hi = np.asarray(inner_casting)
    x,y = np.meshgrid(np.linspace(lo[0]+.5, hi[0]-.5, samples),np.linspace(lo[1]+.5,hi[1]-.5,samples))
    cast = np.c_[x.ravel(),y.ravel(),np.full(x.size,(lo[2]+hi[2])/2)]
    targets = transform(cast,matrix)
    start = targets.copy(); start[:,2] = 0.
    cast_start = transform(start,np.linalg.inv(matrix))
    delta = cast-cast_start
    t = np.divide(hi[2]-cast_start[:,2],delta[:,2],out=np.full(len(cast),np.nan),where=abs(delta[:,2])>1e-9)
    mouth = cast_start+delta*t[:,None]
    crosses = np.isfinite(t)&(t>=0)&(t<=1)&np.all(mouth[:,:2]>lo[:2],axis=1)&np.all(mouth[:,:2]<hi[:2],axis=1)
    valid = crosses & vertical_clear(model_print,targets,.0001)
    return dict(fraction=float(valid.mean()), samples=len(valid), visible=int(valid.sum()),
                points_print=targets.tolist(), visible_mask=valid.tolist(), actual_opening_crossing_mask=crosses.tolist(),
                scope='COARSE_RAYS_NOT_RELEASE_OR_ACTUAL_SUPPORT_PROOF')


def surrogate(row,mold,vascular,box,branches,inputs,profile):
    m=np.asarray(row['print_transform']); r=m[:3,:3]
    model=mold.copy(); model.apply_transform(m)
    b=box.copy(); b.apply_transform(m)
    vc=transform(vascular.triangles_center,m); bc=transform(box.triangles_center,m)
    vdown=(vascular.face_normals@r.T)[:,2]<-np.sqrt(.5)-1e-10
    bdown=(box.face_normals@r.T)[:,2]<-np.sqrt(.5)-1e-10
    vdown &= vc[:,2]>.2; bdown &= bc[:,2]>.2
    tree=cKDTree(np.vstack([q['xyz'] for q in branches])); port=np.concatenate([np.full(len(q['xyz']),q['is_port']) for q in branches])
    is_port=port[tree.query(vascular.triangles_center)[1]]
    # Both faces of the original bottom slab; walls above the inner floor are separate.
    bottom=box.triangles_center[:,2]<=inputs['inner_bounds'][0][2]+1e-4
    clear=vertical_clear(model,vc[vdown],profile['roof_gap_mm'])
    weights=vascular.area_faces[vdown]
    visibility=cavity_visibility(model,inputs['inner_bounds'],m)
    first=section_polygon(model.section_multiplane([0,0,0],[0,0,1],[profile['first_layer_height']/2])[0])
    return dict(**row,bounds_mm=model.bounds.tolist(),z_height_mm=float(model.extents[2]),footprint_mm2=float(np.prod(model.extents[:2])),
                first_layer_model_contact_area_mm2=float(first.area),first_layer_island_count=len(first.geoms) if hasattr(first,'geoms') else int(not first.is_empty),
                box_bottom_unsupported_burden_mm2=float(box.area_faces[bdown&bottom].sum()),
                box_wall_unsupported_burden_mm2=float(box.area_faces[bdown&~bottom].sum()),
                box_ceiling_status='BOX_BOTTOM_SUPPORT_REQUIRED' if np.any(bdown&bottom) else 'NO_BOX_BOTTOM_REQUIRED',
                vascular_unsupported_burden_mm2=float(vascular.area_faces[vdown&~is_port].sum()),
                port_unsupported_burden_mm2=float(vascular.area_faces[vdown&is_port].sum()),
                direct_buildplate_support_fraction=float(weights[clear].sum()/weights.sum()) if weights.sum() else 1.,
                direct_buildplate_fraction_method='AREA_WEIGHTED_DOWNWARD_FACE_PROXY; layer-required fraction recomputed for finalists',
                cavity_buildplate_visibility_fraction=visibility['fraction'],visibility=visibility,
                direct_targets_print=vc[vdown].tolist(),direct_clear_mask=clear.tolist(),surrogate_only=True)


def surrogate_key(row):
    return (row['box_bottom_unsupported_burden_mm2'],-row['direct_buildplate_support_fraction'],
            row['vascular_unsupported_burden_mm2']+row['port_unsupported_burden_mm2'],
            -row['first_layer_model_contact_area_mm2'],row['z_height_mm'],row['orientation_id'])
