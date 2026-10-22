"""Independent exported-geometry/root audits, without regenerating anatomy."""
from pathlib import Path
import json,pickle
import numpy as np
import trimesh
from .support_inputs import save,table,transform
from .support_v3_geometry import union
from .support_v2_mesh import support_stl


def geometry_audit(record,inputs):
    folder=Path(record['folder'])
    support=trimesh.load_mesh(folder/'custom_support_v3.stl');box=trimesh.load_mesh(folder/'box.stl')
    model=trimesh.load_mesh(folder/'mold_without_support.stl');assembly=trimesh.load_mesh(folder/'mold_with_custom_support_v3.stl')
    with (folder/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    repairs=[]
    for name,key in [('support_ribs.stl','ribs'),('support_roofs.stl','roofs'),('support_cut_necks.stl','necks')]:
        exported=trimesh.load_mesh(folder/name)
        if exported.is_volume:continue
        raw=union([r[key] for r in ribs])
        try:
            corrected=support_stl(raw);corrected.export(folder/name)
            repairs.append(dict(file=name,passed=True,cleanup_tolerance_mm=corrected.metadata.get('support_only_numerical_cleanup_tolerance_mm',0.),frozen_anatomy_modified=False))
        except ValueError as error:repairs.append(dict(file=name,passed=False,error=str(error)))
    # Exterior is a spatial domain, not a synonym for the PORT semantic.
    bounds=np.asarray(inputs['outer_bounds']);block=trimesh.creation.box(bounds[1]-bounds[0],transform=trimesh.transformations.translation_matrix(bounds.mean(0)))
    block.apply_transform(np.asarray(record['orientation']['print_transform']))
    outside=trimesh.boolean.difference([support,block],engine='manifold')
    if len(outside.faces):
        try:outside=support_stl(outside)
        except ValueError as error:repairs.append(dict(file='support_exterior.stl',passed=False,error=str(error)))
    outside.export(folder/'support_exterior.stl')
    allowed=union([r['feet'] for r in ribs]);overlap=trimesh.boolean.intersection([support,box],engine='manifold')
    extra=trimesh.boolean.difference([overlap,allowed],engine='manifold') if len(overlap.faces) else overlap
    lost=trimesh.boolean.difference([model,assembly],engine='manifold') if model.is_volume and assembly.is_volume else None
    rows=[];inner=np.array(inputs['inner_bounds']);restore=np.array(record['orientation']['casting_restore_transform'])
    for r in ribs:
        for j,part in enumerate(r['feet'].split()):
            if part.volume<1e-6:continue
            casting=transform([part.centroid],restore)[0]
            internal=bool(np.all(casting[:2]>inner[0,:2]) and np.all(casting[:2]<inner[1,:2]))
            root='BUILD_PLATE' if part.bounds[0,2]<.001 else 'SACRIFICIAL_BOX_SURFACE'
            rows.append(dict(root_id=f"{r['rib_id']}_ROOT{j}",rib_id=r['rib_id'],root_domain=root,
                support_domain='INTERIOR_CAVITY' if internal else 'EXTERIOR_PORT' if r['is_port'] else 'BOX_SURFACE',
                print_xyz=part.centroid.tolist(),casting_xyz=casting.tolist(),approved_root_volume_mm3=float(part.volume)))
    exports=[]
    for name in ['custom_support_v3.stl','mold_without_support.stl','mold_with_custom_support_v3.stl','support_ribs.stl','support_roofs.stl','support_cut_necks.stl','support_exterior.stl']:
        m=trimesh.load_mesh(folder/name);exports.append(dict(file=name,closed_positive_volume=bool(m.is_volume),faces=len(m.faces)))
    result=dict(illegal_box_intersection_mm3=abs(float(extra.volume)) if len(extra.faces) else 0.,semantic_export_repairs=repairs,
        frozen_model_removed_mm3=(abs(float(lost.volume)) if len(lost.faces) else 0.) if lost is not None else None,root_domains=rows,exports=exports,
        approved_root_rule='Only saved foot primitives may overlap the sacrificial box; no anatomical or port union',
        passed=all(x['closed_positive_volume'] for x in exports) and (not len(extra.faces) or abs(float(extra.volume))<1e-4) and lost is not None and (not len(lost.faces) or abs(float(lost.volume))<1e-4))
    save(folder/'exported_geometry_audit.json',result);table(folder/'root_domains.csv',rows)
    record['exported_geometry_audit']=result
    return result
