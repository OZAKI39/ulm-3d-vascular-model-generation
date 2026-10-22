"""Rigid orientation screening for the existing two-port ABS core, without slicing.

Scores are geometric overhang/contact proxies, NOT native support volume or an
assertion of a globally optimal FDM result. All original model files are frozen.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
import trimesh
import yaml
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from . import straight_core_print_quality as q
from .straight_channel_casting import design


def candidates(step, sphere_count):
    if not 0.25 <= step <= 5 or not 100 <= sphere_count <= 5000:
        raise ValueError('UNBOUNDED_ORIENTATION_SEARCH')
    result=[];seen=set()
    def add(name,r):
        # Yaw about print Z cannot change support geometry. Keep one per gravity.
        key=tuple(np.round(r[2],10))
        if key in seen:return
        seen.add(key);result.append((name,r))
    add('horizontal_original',np.eye(3))
    add('vertical_previous',Rotation.from_euler('y',-90,degrees=True).as_matrix())
    add('horizontal_clearance',Rotation.from_euler('x',90,degrees=True).as_matrix())
    add('tilt45',Rotation.from_euler('y',45,degrees=True).as_matrix())
    for angle in np.arange(0.,360.,step):
        add(f'in_plane_{angle:06.2f}',Rotation.from_euler('y',-angle,degrees=True).as_matrix())
    gold=np.pi*(3-np.sqrt(5))
    for i in range(sphere_count):
        z=1-2*(i+.5)/sphere_count;r=np.sqrt(1-z*z)
        direction=np.array([r*np.cos(i*gold),r*np.sin(i*gold),z])
        add(f'sphere_{i:04}',trimesh.geometry.align_vectors(direction,[0,0,1])[:3,:3])
    return result


def place(mesh,rotation,center):
    transform=np.eye(4);transform[:3,:3]=rotation
    moved=mesh.copy();moved.apply_transform(transform)
    shift=np.r_[center,0.]-np.r_[moved.bounds.mean(0)[:2],moved.bounds[0,2]]
    transform[:3,3]=shift;moved.apply_translation(shift)
    return moved,transform


def tube_clearance(mesh,rotation,cfg):
    d=design(cfg);axis=rotation[:,0];g=rotation[2]
    endpoints=np.array([[d['bend_x'],0,d['center_z']],
                        [d['outlet_bend_x'],0,d['center_z']]])
    z=endpoints@g-np.min(mesh.vertices@g)
    radial_vertical=cfg['diameter_mm']/2*np.sqrt(max(0.,1-axis[2]**2))
    return float(min(z)-radial_vertical),float(np.degrees(np.arcsin(abs(axis[2]))))


def support_metrics(mesh,rotation,study,cfg,main_faces):
    g=rotation[2];height=mesh.vertices@g;low=height.min()
    z=mesh.triangles_center@g-low;down=mesh.face_normals@g;areas=mesh.area_faces
    metrics={}
    for angle in study['overhang_angles_deg']:
        mask=(down < -np.cos(np.radians(angle))-1e-7)&(z>study['first_layer_height_mm'])
        projected=areas[mask]*-down[mask]
        metrics[str(angle)]=dict(overhang_surface_mm2=float(areas[mask].sum()),
            projected_contact_mm2=float(projected.sum()),
            height_weighted_projected_proxy_mm3=float(np.sum(projected*np.maximum(0,z[mask]-.24))),
            main_straight_overhang_mm2=float(areas[mask&main_faces].sum()))
    primary=metrics[str(study['ranking_angle_deg'])]
    score=(primary['overhang_surface_mm2']+
           primary['height_weighted_projected_proxy_mm3']/study['score_height_normalizer_mm']+
           study['main_straight_area_weight']*primary['main_straight_overhang_mm2'])
    clearance,angle=tube_clearance(mesh,rotation,cfg)
    return dict(score=score,metrics_by_angle=metrics,main_clearance_mm=clearance,
                main_angle_to_bed_deg=angle,print_height_mm=float(np.ptp(height)))


def first_layer(mesh,height):
    # A single geometric plane section for adhesion-area screening; no extrusion
    # paths, print layers, support generation or native slicer are invoked.
    section=mesh.section([0,0,1],[0,0,height])
    polygons=[] if section is None else list(section.to_2D()[0].polygons_full)
    return dict(area_mm2=float(sum(p.area for p in polygons)),components=len(polygons),
                z_mm=height,actual_adhesion_tested=False)


def pose_audit(before,after,transform,cfg):
    restored=after.copy();restored.apply_transform(np.linalg.inv(transform))
    err=float(abs(restored.vertices-before.vertices).max())
    if (err>1e-7 or not np.array_equal(before.faces,after.faces)
            or not after.is_volume or not np.isclose(before.volume,after.volume,rtol=1e-8)):
        raise ValueError('RIGID_CORE_GEOMETRY_CHECK_FAILED')
    if abs(after.bounds[0,2])>1e-8 or np.any(after.bounds[0,:2]<20) or np.any(after.bounds[1,:2]>236):
        raise ValueError('INVALID_BED_PLACEMENT')
    clearance,angle=tube_clearance(before,transform[:3,:3],cfg)
    return dict(passed=True,diameter_mm=cfg['diameter_mm'],maximum_inverse_error_mm=err,
        same_triangles=True,source_volume_mm3=before.volume,print_volume_mm3=after.volume,
        main_straight_clearance_mm=clearance,main_straight_angle_to_bed_deg=angle,
        print_height_mm=float(after.extents[2]),casting_to_print_transform=transform,
        print_to_casting_transform=np.linalg.inv(transform),source_geometry_changed=False,
        extra_feet_or_support_meshes_added=False)


def comparisons(rows,meshes,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    titles=['Previous: main tube vertical','Low-support candidate: about 45 deg','Horizontal: collars on bed, tube raised']
    for column,(row,mesh,title) in enumerate(zip(rows,meshes,titles)):
        v=mesh.vertices.copy();v[:,:2]-=mesh.bounds.mean(0)[:2]
        for j,indices in enumerate([[0,2],[0,1]]):
            ax=axes[j,column]
            ax.add_collection(PolyCollection(v[mesh.faces][:,:,indices],facecolors='#c45142',edgecolors='none'))
            ax.set_aspect('equal');ax.set_xlim(-22,22);ax.set_xlabel('Print X (mm)');ax.grid(alpha=.15)
            if j==0:
                ax.set_ylim(-2,36);ax.set_ylabel('Print Z (mm)');ax.axhline(0,color='#46546a',lw=2)
                ax.set_title(title+'\nTube clearance %.2f mm'%row['main_clearance_mm'])
            else:ax.set_ylim(-20,20);ax.set_ylabel('Print Y (mm)');ax.set_title('Top view')
    fig.suptitle('Frozen 1.60 mm ABS core: only orientation changes\nGeometric screening only; no sliced support or G-code')
    fig.tight_layout();fig.savefig(out/'orientation_comparison.png',dpi=170);plt.close(fig)


def package_study(out):
    out=Path(out);summary=json.loads((out/'orientation_summary.json').read_text())
    files=[out/n for n in ['orientation_summary.json','orientation_candidates.csv',
          'orientation_candidates.json','orientation_comparison.png','effective_study.yaml','使用说明.md']]
    for record in summary['deliveries']:
        folder=out/record['folder'];p=folder/record['project']['file']
        if q.digest(p)!=record['project']['sha256']:raise ValueError('STALE_ORIENTATION_PROJECT')
        files.extend(folder/n for n in [p.name,'print_oriented_core.stl','effective_settings.json','geometry_qc.json'])
        if record.get('gui_readback_passed'):
            files.extend(folder/n for n in ['gui_verification.json','bambu_gui_verified.png'])
    hashes={str(p.relative_to(out)):q.digest(p) for p in files}
    dest=q.ROOT/'deliverables/Straight_2Top_P1S04_ABS_orientation_options_unsliced.zip'
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,str(p.relative_to(out)))
        z.writestr('SHA256.json',json.dumps(hashes,indent=2))
    with zipfile.ZipFile(dest) as z:
        if z.testzip() is not None:raise ValueError('ORIENTATION_ZIP_CRC_FAILED')
        for n,h in hashes.items():
            if hashlib.sha256(z.read(n)).hexdigest()!=h:raise ValueError('ORIENTATION_ZIP_HASH_FAILED')
    q.save(out/'delivery_receipt.json',dict(zip=str(dest.relative_to(q.ROOT)),sha256=q.digest(dest),
           members_sha256=hashes,slicing_performed=False))
    return dest


def run_study(config_path):
    config_path=Path(config_path);study=yaml.safe_load(config_path.read_text())
    quality=yaml.safe_load((q.ROOT/study['quality_config']).read_text())
    cfg=q.load_config(q.ROOT/quality['source_config']);source_out=cfg['output']
    source=source_out/'Straight_2Top_P1S_04_core_ABS_project.3mf'
    assembly=source_out/'assembly/Straight_lumen_core.stl'
    if q.digest(source)!=quality['source_project_sha256'] or q.digest(assembly)!=quality['source_assembly_sha256']:
        raise ValueError('FROZEN_SOURCE_HASH_MISMATCH')
    before=q.read_objects(source)[0][0]['mesh']
    record=next(r for r in json.loads((source_out/'print_layout.json').read_text())['records'] if r['name']==q.NAME)
    mesh=before.copy();mesh.apply_transform(np.linalg.inv(record['casting_to_plate_transform']))
    out=(q.ROOT/study['output']).resolve()
    if out==source_out or source_out in out.parents:raise ValueError('SEPARATE_STUDY_OUTPUT_REQUIRED')
    out.mkdir(parents=True,exist_ok=True)
    d=design(cfg);p=mesh.triangles_center
    main=(p[:,0]>d['bend_x']+1e-4)&(p[:,0]<d['outlet_bend_x']-1e-4)&(p[:,2]<d['center_z']+cfg['diameter_mm']/2+1e-4)
    rows=[];rotations={}
    with threadpool_limits(limits=1):
        for name,r in candidates(study['in_plane_step_deg'],study['sphere_directions']):
            row=support_metrics(mesh,r,study,cfg,main);row['id']=name;row['rotation']=r.tolist()
            rows.append(row);rotations[name]=r
    ranked=sorted(rows,key=lambda r:(r['score'],r['id']));selected=None
    # Apply real first-layer footprint and critical tube clearance gates, rather
    # than choosing a low-overhang orientation that cannot stay on the plate.
    for row in ranked:
        if row['main_clearance_mm']<study['minimum_main_straight_clearance_mm']:continue
        moved,_=place(mesh,rotations[row['id']],study['plate_center_mm'])
        row['first_layer']=first_layer(moved,study['first_layer_height_mm'])
        if row['first_layer']['area_mm2']>=study['minimum_first_layer_area_mm2']:
            selected=row;break
    if selected is None:raise ValueError('NO_FEASIBLE_LOW_SUPPORT_POSE')
    horizontal=next(r for r in rows if r['id']=='horizontal_clearance')
    previous=next(r for r in rows if r['id']=='vertical_previous')
    deliveries=[]
    for role,row in [('01_low_support',selected),('02_horizontal_clearance',horizontal)]:
        folder=out/role;folder.mkdir(exist_ok=True)
        moved,T=place(mesh,rotations[row['id']],study['plate_center_mm']);audit=pose_audit(mesh,moved,T,cfg)
        if role=='02_horizontal_clearance':
            if abs(audit['main_straight_angle_to_bed_deg'])>1e-7 or not 1.15<audit['main_straight_clearance_mm']<1.25:
                raise ValueError('HORIZONTAL_TUBE_CLEARANCE_FAILED')
        settings=dict(quality['settings']);settings['support_threshold_angle']='45'
        settings['print_settings_id']='Straight 2Top '+role+' 0.08 - P1S'
        if role=='02_horizontal_clearance':settings.update(support_type='normal(auto)',support_style='snug')
        filename='Straight_2Top_P1S04_ABS_'+role.split('_',1)[1]+'_project.3mf'
        q.write_refined_project(source,folder/filename,moved,settings,
            plate_name='Straight 2Top - '+role.replace('_',' ')+' - UNSLICED')
        project=q.validate_project(folder/filename,{q.NAME:moved},'ABS',cfg)
        moved.export(folder/'print_oriented_core.stl');q.save(folder/'geometry_qc.json',audit)
        q.save(folder/'effective_settings.json',dict(settings=settings))
        q.save(folder/'orientation.json',row)
        for name in ['gui_verification.json','bambu_gui_verified.png']:(folder/name).unlink(missing_ok=True)
        row['first_layer']=first_layer(moved,study['first_layer_height_mm'])
        deliveries.append(dict(folder=role,candidate_id=row['id'],project=project,geometry=audit,
                               metrics=row,gui_readback_passed=False))
    compare=[previous,selected,horizontal]
    comparisons(compare,[place(mesh,rotations[r['id']],study['plate_center_mm'])[0] for r in compare],out)
    q.save(out/'orientation_candidates.json',rows)
    with (out/'orientation_candidates.csv').open('w',newline='') as f:
        fields=['id','score','main_clearance_mm','main_angle_to_bed_deg','print_height_mm','area45_mm2','height_weighted_proxy_mm3']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r in ranked:
            v={k:r[k] for k in fields[:5]};m=r['metrics_by_angle']['45']
            v.update(area45_mm2=m['overhang_surface_mm2'],height_weighted_proxy_mm3=m['height_weighted_projected_proxy_mm3']);writer.writerow(v)
    summary=dict(status='UNSLICED_ORIENTATION_OPTIONS_PHYSICAL_TEST_REQUIRED',tested_candidates=len(rows),
        selected=selected,horizontal=horizontal,previous=previous,deliveries=deliveries,
        metric='Surface angle screening plus height-weighted projected overhang; NOT slicer support volume',
        score_formula='A45 + projected_height_proxy/20 mm + 2*main_straight_A45',
        global_optimality_proven=False,generated_support_meshes=0,slicer_executed=False,
        source_project_sha256=q.digest(source),source_assembly_sha256=q.digest(assembly),
        old_model_files_preserved=True,physical_print_tested=False)
    q.save(out/'orientation_summary.json',summary)
    (out/'effective_study.yaml').write_text(config_path.read_text())
    (out/'使用说明.md').write_text((q.ROOT/'docs/Straight_2Top少支撑与水平离板姿态.md').read_text())
    dest=package_study(out)
    print('ORIENTATION_STUDY',len(rows),'BEST',selected['id'],'ZIP',dest,flush=True)
    return summary
