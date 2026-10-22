"""Bounded local changes to the existing explicit s1-11 scaffold only."""
from pathlib import Path
import json,pickle,shutil
import numpy as np
import trimesh
from .core_support_geometry import union,boolean,overlap,bed_roots
from .core_support_closure import canonical_support,local_box,quantized_top
from .support_inputs import save,digest


def material_components(mesh):
    """Physical solids including their void shells, not raw STL face islands."""
    import manifold3d as md
    m=md.Manifold(md.Mesh64(np.array(mesh.vertices,dtype=np.float64,order='C'),np.array(mesh.faces,dtype=np.uint64,order='C')))
    if m.status()!=md.Error.NoError:raise ValueError('SUPPORT_MATERIAL_GRAPH_NOT_MANIFOLD')
    parts=[];voids=[]
    for c in m.decompose():
        if abs(c.volume())<=1e-12:continue
        v=c.to_mesh64();part=trimesh.Trimesh(v.vert_properties[:,:3],v.tri_verts,process=False)
        (parts if c.volume()>0 else voids).append(part)
    cavities={i:[] for i in range(len(parts))}
    for void in voids:
        owners=[i for i,p in enumerate(parts) if np.all(void.bounds[0]>=p.bounds[0]-1e-7) and np.all(void.bounds[1]<=p.bounds[1]+1e-7) and p.contains([void.centroid])[0]]
        if not owners:raise ValueError('UNASSIGNED_INTERNAL_VOID_SHELL')
        i=min(owners,key=lambda j:abs(parts[j].volume));cavities[i].append(void)
    parts=[trimesh.util.concatenate([p,*cavities[i]]) if cavities[i] else p for i,p in enumerate(parts)]
    return sorted(parts,key=lambda p:-p.volume)


def write_candidate(parent,folder,support,metadata):
    parent=Path(parent);folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    g=pickle.load((parent/'geometry.pkl').open('rb'));support,encoding=canonical_support(support)
    for name in ['vascular_core.stl','vascular_print_transform.json','inverse_transform.json','notches.json']:shutil.copy2(parent/name,folder/name)
    support.export(folder/'custom_support.stl',file_type=encoding)
    trimesh.util.concatenate([g['vascular'],support]).export(folder/'vascular_core_plus_support_review.stl',file_type=encoding)
    g['support']=support;shifted=support.copy();shifted.apply_translation([0,0,-3*g['profile']['layer_height']]);g['roofs']=boolean([support,shifted],'difference');g['curtains']=boolean([support,g['roofs']],'difference')
    g['practical_metadata']=metadata
    with (folder/'geometry.pkl').open('wb') as f:pickle.dump(g,f)
    r=json.loads((parent/'geometry.json').read_text());r.update(candidate_id=folder.name,folder=str(folder),practical_parent=str(parent),practical_change=metadata,support_volume_mm3=float(support.volume),intersection_mm3=overlap(support,g['vascular']),bed_roots=bed_roots(support),support_stl_encoding=encoding)
    save(folder/'geometry.json',r)
    assert digest(folder/'vascular_core.stl')==digest(parent/'vascular_core.stl')
    return r


def repair_true_missing(parent,folder,structural,round_number):
    if round_number not in [1,2,3]:raise ValueError('THREE_REPAIR_ROUNDS_MAXIMUM')
    parent=Path(parent);g=pickle.load((parent/'geometry.pkl').open('rb'));p=g['profile'];w=p['support_line_width'];h=p['layer_height'];v=g['vascular'];s=g['support'];branches={b['branch_id']:b for b in g['branches']};patches=[];changes=[];blocked=[]
    for row in structural['rows']:
        if not row['true_missing_support']:continue
        q=np.asarray(row['xyz']);b=branches[row['branch']];t=b['tangent'][row['sample_index'],:2];t=t/max(np.linalg.norm(t),1e-10);side=np.array([-t[1],t[0]])
        top=quantized_top(q[2],p['roof_gap_mm'],p);eps=4*np.spacing(np.float32(max(1.,np.abs(v.bounds).max())));mesh_top=top+h/2-eps
        width=(2 if round_number==1 else 3)*w;length=3*w
        shifts=[0.] if round_number<3 else [0.,w,-w]
        attempts=[]
        for shift in shifts:
            patch=local_box(q[:2]+shift*side,mesh_top,width,length,t,max(0.,mesh_top-(4 if round_number<3 else 8)*h))
            attempts.append((patch,shift,width,length,4 if round_number<3 else 8))
        if round_number==3 and row.get('nearest_support_path'):
            # At a nearby vascular crossing, a three-line rectangular pad can
            # intersect another branch. Extend only the actual missed bead lip
            # by one line, with one-layer depth, rather than moving the roof.
            near=np.asarray(row['nearest_support_path']['nearest_footprint_xy'])
            direction=q[:2]-near;distance=np.linalg.norm(direction)
            if distance>1e-10:
                direction/=distance
                lip=local_box((near+q[:2])/2,mesh_top,w,w+distance,direction,max(0.,mesh_top-h))
                attempts.append((lip,0.,w,w+distance,1))
        for patch,shift,patch_width,patch_length,depth in attempts:
            if overlap(patch,v)>=1e-7 or overlap(patch,s)<=1e-9:continue
            patches.append(patch);changes.append(dict(sample_id=row['sample_id'],new_class=row['new_class'],xyz=row['xyz'],width_lines=patch_width/w,length_lines=patch_length/w,shift_lines=shift/w,depth_layers=depth,roof_layer_top_mm=top,kind='LOCAL_TRUE_MISSING_ONLY'));break
        else:blocked.append(dict(sample_id=row['sample_id'],reason='NO_SHORT_ATTACHED_VESSEL_CLEAR_PATCH'))
    support=union([s,*patches]) if patches else s.copy()
    metadata=dict(stage='TRUE_MISSING_REPAIR',round=round_number,changes=changes,blocked=blocked,too_close_targets_modified=False,global_roof_regeneration=False)
    result=write_candidate(parent,folder,support,metadata);save(Path(folder)/'true_missing_repairs.json',metadata)
    return result
