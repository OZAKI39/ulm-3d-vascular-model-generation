"""Change only the selected mold's wall/bottom thickness; export STL, no slicing."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import trimesh

from .support_aware_roi import DEFAULT_OUTPUT, ROOT, sha256, table, verify_snapshot
from .support_aware_port_design import configuration, box_from_bounds, load_companion, delivered_interface_qc
from .sacrificial_fixture_review import save_json
from .sacrificial_fixture import export_checked_stl
from . import abs_casting_mold as mold
from . import casting_mold_qc as qc
from .print_orientation import bed_transform


def resized_box(inner,thickness,cfg):
    """Keep the cavity and all vascular clearances; change external faces only."""
    thickness=float(thickness);inner=np.asarray(inner,dtype=float)
    if not np.isfinite(thickness) or thickness<=0:
        raise ValueError('BOX_THICKNESS_MUST_BE_POSITIVE_AND_FINITE')
    if inner.shape!=(2,3) or not np.isfinite(inner).all() or np.any(inner[1]<=inner[0]):
        raise ValueError('INVALID_INNER_BOX_BOUNDS')
    outer=inner.copy();outer[0]-=[thickness,thickness,thickness];outer[1]+=[thickness,thickness,0.]
    cfg['box'].update(policy='BOX_DIMENSIONS_REUSED',wall_thickness_mm=thickness,bottom_thickness_mm=thickness)
    result=box_from_bounds(inner,outer,cfg)
    result['dimensions_policy']='PRESERVE_SELECTED_INNER_BOUNDARIES_CHANGE_OUTER_THICKNESS_ONLY'
    return result


def measured_thickness(mesh):
    """Measure planes of the exported five-wall STL, independent of CAD metadata."""
    planes=[np.unique(mesh.vertices[:,axis]) for axis in range(3)]
    if [len(v) for v in planes]!=[4,4,3]:raise ValueError('UNEXPECTED_FIVE_WALL_STL_PLANES')
    x,y,z=planes
    return dict(negative_x_mm=float(x[1]-x[0]),positive_x_mm=float(x[-1]-x[-2]),
        negative_y_mm=float(y[1]-y[0]),positive_y_mm=float(y[-1]-y[-2]),bottom_mm=float(z[1]-z[0]))


def export_box_stl(root=DEFAULT_OUTPUT,thickness=2.):
    root=Path(root).resolve();cfg=configuration()
    summary_path=root/'roi_60x50_summary.json';summary=json.loads(summary_path.read_text())
    candidate=root/'candidates'/summary['selected_candidate'];source_dir=candidate/'manufacturing'
    geometry_path=source_dir/'geometry.json';geometry=json.loads(geometry_path.read_text())
    if geometry['failures']:raise ValueError('SELECTED_MANUFACTURING_GEOMETRY_HAS_FAILURES')
    box=resized_box(geometry['box']['inner_bounds_mm'],thickness,cfg)
    out=root/f'casting_mold_{thickness:g}mm';out.mkdir(exist_ok=True)
    # Snapshot all existing selected-model files, including the original 3 mm
    # slice results. Loading helpers that rewrite historical tables are avoided.
    protected=[summary_path,geometry_path,root/'roi_casting_transform.json',root/'resolved_evaluation_config.json']
    for directory in (source_dir,root/'selected',root/'casting_mold'):
        protected.extend(p for p in directory.rglob('*') if p.is_file())
    before={str(p):sha256(p) for p in protected};save_json(out/'protected_before.json',before)
    core=trimesh.load_mesh(source_dir/'core_with_ports.stl')
    ports=[]
    for original in geometry['ports']:
        port=dict(original)
        for key in ('target','end','direction'):port[key]=np.asarray(port[key],dtype=float)
        port['reference_mesh']=trimesh.load_mesh(source_dir/(port['port_id']+'.stl'))
        ports.append(port)
    print('Building box and union from the unchanged selected vascular core...',flush=True)
    box_path=out/f'intact_open_top_box_{thickness:g}mm.stl'
    union_path=out/f'final_abs_casting_mold_{thickness:g}mm.stl'
    oriented_path=out/f'final_abs_casting_mold_{thickness:g}mm_oriented.stl'
    export_checked_stl(box['mesh'],box_path,cfg)
    # QC uses reloaded delivered triangles, not just the CAD in memory.
    box['mesh']=trimesh.load_mesh(box_path)
    actual=measured_thickness(box['mesh'])
    if not all(abs(v-thickness)<2e-5 for v in actual.values()):raise ValueError('EXPORTED_STL_THICKNESS_MISMATCH')
    union,stages=mold.combine(core,box['mesh'],cfg);export_checked_stl(union,union_path,cfg)
    union=trimesh.load_mesh(union_path)
    box_qc=qc.intact_box_qc(box,cfg);audit=qc.union_audit(core,box['mesh'],union,cfg)
    overlap,_=qc.port_overlaps(dict(core=core,ports=ports),box,cfg)
    transform=np.asarray(json.loads((root/'roi_casting_transform.json').read_text())['transform_4x4'])
    companion=load_companion(candidate,transform,cfg,align=True)
    interfaces=delivered_interface_qc(union,companion,geometry['interfaces'],cfg)
    failures=audit['failures'].copy()
    if not all(box_qc[k] for k in ('four_side_walls_closed','bottom_closed','top_fully_open','no_assembly_clearance_holes')):failures.append('BOX_STRUCTURE_FAILED')
    if any(p['status']!='PASS' for p in overlap):failures.append('PORT_WALL_CONNECTION_FAILED')
    if not all(p['passed'] for p in interfaces):failures.append('PORT_INTERFACE_FAILED')
    rotation=np.asarray(summary['orientation']['transform_4x4'])[:3,:3]
    mat,_,_,_=bed_transform(union.vertices,rotation,cfg['printer']['provisional_build_volume_mm'],10.)
    placed=union.copy();placed.apply_transform(mat);export_checked_stl(placed,oriented_path,cfg)
    source_check=verify_snapshot(before)
    if not source_check['all_unchanged']:failures.append('SOURCE_MODEL_CHANGED')
    extension=[dict(port_id=p['port_id'],outside_stub_mm=float((p['end']-p['target'])@p['direction']-thickness)) for p in ports]
    result=dict(status='GEOMETRY_STL_CHECKS_PASSED' if not failures else 'GEOMETRY_STL_CHECKS_FAILED',failures=failures,
        selected_candidate=summary['selected_candidate'],source_summary=str(summary_path),source_core_sha256=sha256(source_dir/'core_with_ports.stl'),
        wall_thickness_mm=thickness,bottom_thickness_mm=thickness,measured_stl_thickness_mm=actual,
        inner_dimensions_mm=box['inner'][1]-box['inner'][0],outer_dimensions_mm=box['outer'][1]-box['outer'][0],
        inner_bounds_mm=box['inner'],outer_bounds_mm=box['outer'],dimensions_policy=box['dimensions_policy'],
        box_qc=box_qc,union_qc=audit,port_overlaps=overlap,interfaces=interfaces,source_protection=source_check,
        vascular_and_port_geometry_preserved=True,exposed_port_lengths=extension,print_transform_4x4=mat,
        inherited_roi_status=summary['status'],inherited_roi_failures=summary['failures'],
        slicing_status='NOT_RESLICED_FOR_NEW_THICKNESS',printer_job_sent=False,
        files={p.name:sha256(p) for p in (box_path,union_path,oriented_path)})
    save_json(out/'geometry_summary.json',result);save_json(out/'resolved_box_config.json',cfg['box'])
    save_json(out/'source_protection.json',source_check)
    for name,rows in [('actual_thickness', [actual]),('port_wall_overlap',overlap),('port_interfaces',interfaces)]:table(out/(name+'.csv'),rows)
    if failures:raise ValueError('STL_EXPORT_QC_FAILED: '+', '.join(failures))
    from .support_aware_review import scene
    scene([(box['mesh'],'#99b2c3',.5),(core,'#399d63',1)],out/'box_preview.png',f'自适应盒体：侧壁和底板 {thickness:g} mm',
        ["盒体外尺寸："+' × '.join(f'{v:.2f}' for v in result['outer_dimensions_mm'])+' mm','内腔及血管、端口位置保持；顶部开放'])
    (out/'README.md').write_text(f'''# 当前自适应盒体：{thickness:g} mm STL

侧壁及底板均为 **{thickness:g} mm**。内腔尺寸保持为 {' × '.join(f'{v:.2f}' for v in result['inner_dimensions_mm'])} mm，盒体外尺寸为 **{' × '.join(f'{v:.2f}' for v in result['outer_dimensions_mm'])} mm**。

中央血管和四根延长段直接使用既有选中模型，未重新拟合或选路。内壁位置保持，侧壁向外侧减薄、底板从底面减薄。由于外壁退回，墙外可见接口增加至 {'、'.join(f"{p['port_id']} {p['outside_stub_mm']:.2f} mm" for p in extension)}，端口总长度没有改变。

- `{union_path.name}`：包含血管、端口和盒体的一体 STL，保留铸造坐标。
- `{oriented_path.name}`：相同模型，按已选旋转姿态重新居中、底面落在 Z=0，供导入切片器。
- `{box_path.name}`：单独盒体。
- `geometry_summary.json`、CSV：实际 STL 厚度、端口穿墙连接、接口和材料保留检查。

模型为单连通、闭合且流形的实体，血管材料损失 {audit['core_material_loss_mm3']:.8f} mm³。各端口穿墙检查和接口检查均通过。此前 ROI 未达到 60 × 50 mm 容差的限制仍保留。

本次仅重新生成 STL，没有针对新厚度切片。旧目录中的 3MF/G-code 对应原 3 mm 版本，不是本次 2 mm STL 的切片结果。

复现：`./run_python.sh s1-7_support_aware_roi_60x50.py --stage box-stl --box-thickness-mm {thickness:g}`。
''',encoding='utf8')
    save_json(out/'implementation_manifest.json',{str(p.relative_to(ROOT)):sha256(p) for p in (ROOT/'s1-7_support_aware_roi_60x50.py',Path(__file__))})
    print(result['status'],str(oriented_path),flush=True)
    return result
