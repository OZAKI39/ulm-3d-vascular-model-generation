#!/usr/bin/env python3
"""Export candidates 0 and 15 from verified, unchanged real Bambu slices."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import trimesh
from vascular_processing import abs_casting_mold as mold
from vascular_processing import sacrificial_fixture as fixture
from vascular_processing import bambu_manufacturing as bambu
from vascular_processing import casting_mold_review as review
from vascular_processing.sacrificial_fixture_review import save_json
from vascular_processing.project_paths import verify_migrated_snapshot

DEFAULT=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold'


def export(source,output):
    source=source.resolve();output=output.resolve()
    if output.exists():raise FileExistsError('Refusing to overwrite an existing candidate package: '+str(output))
    summary_path=source/'abs_casting_mold_summary.json'
    summary=json.loads(summary_path.read_text());cfg=mold.load_config(source/'resolved_config.yaml')
    orientations={r['candidate_id']:r for r in summary['orientation']['top_five']}
    orientations[0]=summary['baseline']['orientation']
    plans=[];protected={str(summary_path):fixture.sha256(summary_path)}
    for candidate_id,slug,title in [(0,'upright','正放：底板贴打印板，顶部开口朝上'),(15,'side','侧放：+X 侧壁贴打印板，开口朝侧方')]:
        orientation=orientations[candidate_id]
        records=[r for r in summary['bambu']['results'] if r['orientation_id']==candidate_id and r['support_mode']=='ON']
        if len(records)!=1:raise ValueError('Missing or ambiguous support-ON slice for candidate '+str(candidate_id))
        record=records[0];stl=Path(orientation['stl']);archive=Path(record['archive']);gcode=Path(record['gcode_file'])
        if record['status']!='BAMBU_SLICED' or not record['feature_preservation']['passed']:
            raise ValueError('Actual slicing or feature preservation failed')
        if fixture.sha256(stl)!=record['geometry_qc']['input_stl_sha256']:
            raise ValueError('STL differs from actual slicer input')
        if fixture.sha256(archive)!=record['geometry_qc']['sliced_3mf_sha256']:
            raise ValueError('3MF archive differs from audited slice')
        with zipfile.ZipFile(archive) as z:
            if z.testzip() is not None:raise ValueError('3MF CRC failure')
            settings=json.loads(z.read('Metadata/project_settings.config'))
            if str(settings.get('enable_support'))!='1':raise ValueError('Support is not enabled')
            names=[n for n in z.namelist() if n.endswith('.gcode')]
            if len(names)!=1 or hashlib.sha256(z.read(names[0])).hexdigest()!=fixture.sha256(gcode):
                raise ValueError('External G-code differs from the real 3MF payload')
        for path in (stl,archive,gcode,Path(record['native_preview'])):protected[str(path)]=fixture.sha256(path)
        plans.append((candidate_id,slug,title,orientation,record,stl,archive,gcode))
    output.mkdir(parents=True);items=[]
    for candidate_id,slug,title,orientation,record,stl,archive,gcode in plans:
        folder=output/f'candidate_{candidate_id:02}_{slug}';folder.mkdir()
        basename=f'BG001_RMCA_candidate_{candidate_id:02}_{slug}'
        copies={basename+'.stl':stl,basename+'_support_ON.3mf':archive,
                basename+'_support_ON.gcode':gcode,'native_bambu_preview.png':Path(record['native_preview'])}
        for name,path in copies.items():
            shutil.copy2(path,folder/name)
            if fixture.sha256(folder/name)!=fixture.sha256(path):raise ValueError('Export copy checksum mismatch')
        mesh=trimesh.load_mesh(folder/(basename+'.stl'));geometry=fixture.mesh_qc(mesh,cfg)
        if not geometry['passed']:raise ValueError('Exported STL failed mesh QC')
        if not np.allclose(mesh.bounds,[orientation['bbox_min_mm'],orientation['bbox_max_mm']],atol=2e-5,rtol=0):
            raise ValueError('Exported STL orientation mismatch')
        archive_geometry=bambu.sliced_geometry_qc(folder/(basename+'_support_ON.3mf'),folder/(basename+'.stl'))
        if not archive_geometry['rotation_scale_position_preserved']:raise ValueError('3MF geometry differs from exported STL')
        parsed=bambu.parse_sliced_archive(folder/(basename+'_support_ON.3mf'))
        if not parsed['slice_success']:raise ValueError('Exported 3MF has no real extrusion paths')
        save_json(folder/'print_transform.json',orientation)
        save_json(folder/'casting_restore_transform.json',dict(transform_4x4=orientation['casting_restore_transform'],casting_top_direction=[0,0,1]))
        save_json(folder/'original_slice_audit.json',record)
        fixture.write_csv(folder/'slicer_feature_preservation.csv',record['feature_preservation']['rows'])
        image=review.scene(mesh,bounds=mesh.bounds,bed=True,view=review.printing_view(orientation))
        review.caption(image,f'Candidate {candidate_id}：'+title,
            ['模型与对应真实 Bambu 切片均保持原样；文件已开启树状支撑。'],cfg).save(folder/'print_orientation.png')
        risk=orientation['trapped_risk_count']
        manifest=dict(candidate_id=candidate_id,source_rank=orientation['rank'],orientation=title,
            support_mode='ON',support_type=parsed['actual_project_settings']['support_type'],
            actual_project_settings=parsed['actual_project_settings'],print_time_seconds=parsed['estimated_print_time_seconds'],
            total_filament_g=parsed['filament_used_g'],estimated_internal_support_g=record['support_quantity']['estimated_internal_support_mass_g'],
            support_removal_proxy_risk_count=risk,real_support_removal_certified=False,
            status='SUPPORT_REMOVAL_PROXY_REVIEW_REQUIRED' if risk else 'GEOMETRY_PROXY_SUPPORT_REMOVAL_PASS',
            mesh_qc=geometry,archive_geometry_qc=archive_geometry,sampled_gcode_features_preserved=True,
            gcode_matches_3mf_payload=True,source_files={name:str(path) for name,path in copies.items()},
            file_sha256={p.name:fixture.sha256(p) for p in folder.iterdir() if p.is_file()},
            geometry_modified=False,resliced=False,printer_job_sent=False)
        save_json(folder/'manifest.json',manifest)
        seconds=parsed['estimated_print_time_seconds'];minutes=seconds/60
        text=f'''# Candidate {candidate_id}\n\n{title}。\n\n请在 Bambu Studio 中打开 `{basename}_support_ON.3mf` 检查真实切片；`{basename}.stl` 是同一打印姿态的整体模型，`{basename}_support_ON.gcode` 与 3MF 中的 G-code 完全一致。\n\n配置：Bambu Lab P1S、0.4 mm 喷嘴、Bambu ABS、0.20 mm Standard；树状支撑已开启。预计 {minutes:.1f} 分钟，总耗材 {parsed['filament_used_g']:.2f} g，盒内支撑估计 {manifest['estimated_internal_support_g']:.2f} g。\n\n顶部清理通路代理判定有 {risk} 个受阻区域。这个数值不等同于实物支撑一定被困；0 个也不代表已验证支撑可以拆除。两种姿态均需人工检查，尚未进行实物拆支撑验证。\n\n原始切片内部对象名使用历史排序编号；当前文件名和 manifest 中的 candidate_id 才是本包的候选编号。Candidate 15 对应历史 rank 4，未修改原始 3MF 的内部对象名。\n\nprint_transform.json 是倒模坐标到打印板坐标的刚性变换；casting_restore_transform.json 为逆变换。侧放件在打印后需转回开口朝上再灌注。血管、半径、端口、盒体和原始文件均未修改，也没有发送打印任务。\n'''
        (folder/'README.md').write_text(text,encoding='utf-8')
        items.append(dict(candidate_id=candidate_id,directory=folder.name,orientation=title,
            print_time_minutes=minutes,total_filament_g=parsed['filament_used_g'],
            internal_support_g=manifest['estimated_internal_support_g'],support_removal_proxy_risk_count=risk,
            stl=basename+'.stl',three_mf=basename+'_support_ON.3mf',gcode=basename+'_support_ON.gcode'))
        print(f'Candidate {candidate_id}: STL / 3MF / G-code verified; {minutes:.1f} min; {parsed["filament_used_g"]:.2f} g',flush=True)
    protection=fixture.verify_snapshot(protected)
    upstream=verify_migrated_snapshot(json.loads((source/'protected_before.json').read_text()))
    if not protection['all_unchanged'] or not upstream['passed']:raise ValueError('Protected source hash changed beyond recorded migration edits')
    save_json(output/'source_protection.json',dict(export_sources=protection,upstream_sources=upstream))
    save_json(output/'export_summary.json',dict(source_summary=str(summary_path),candidates=items,
        source_files_unchanged=True,real_slicer_archives_reused=True,printer_job_sent=False))
    fixture.write_csv(output/'candidate_comparison.csv',items)
    rows='\n'.join(f"| {r['candidate_id']} | {r['orientation']} | {r['print_time_minutes']:.1f} | {r['total_filament_g']:.2f} | {r['internal_support_g']:.2f} | {r['support_removal_proxy_risk_count']} |" for r in items)
    (output/'README.md').write_text('# Candidate 0 与 Candidate 15 独立打印文件\n\n两套文件均使用已开启树状支撑的真实 Bambu 切片。每个目录包含对应姿态的 STL、3MF、G-code、姿态图、变换、校验结果和说明。\n\n| 候选 | 姿态 | 分钟 | 总耗材 g | 盒内支撑估计 g | 代理受阻区域 |\n|---|---|---:|---:|---:|---:|\n'+rows+'\n\n请分别打开各目录中的 `_support_ON.3mf`。这里只导出两种姿态供比较，没有把 candidate 0 的代理检查结果改成通过；两种姿态均未完成实际支撑拆除验证。已逐字节核对 STL、3MF、G-code，重新核对了 3MF 模型姿态与 STL 的一致性。\n',encoding='utf-8')
    for item in items:
        folder=output/item['directory'];archive=output/(item['directory']+'.zip')
        with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for path in sorted(folder.rglob('*')):
                if path.is_file():z.write(path,path.relative_to(output))
        with zipfile.ZipFile(archive) as z:
            if z.testzip() is not None:raise ValueError('Candidate ZIP CRC failure')
    sums={str(p.relative_to(output)):fixture.sha256(p) for p in output.rglob('*') if p.is_file()}
    (output/'SHA256SUMS').write_text(''.join(f'{h}  {name}\n' for name,h in sorted(sums.items())),encoding='utf-8')
    return items


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=DEFAULT)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();export(args.source,args.output or args.source/'candidate_0_and_15')
