#!/usr/bin/env python3
"""Create ordinary portable ZIP files from verified candidate 0/15 deliveries.

Geometry and G-code are copied byte for byte. This is packaging, not slicing.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def sha256(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def write_json(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def repack(source_root,output):
    source_root=Path(source_root).resolve();output=Path(output).absolute()
    if output.exists():raise FileExistsError(output)
    plans=[]
    for cid in (0,15):
        source=(source_root/f'candidate_{cid}').resolve(strict=True)
        source_manifest=json.loads((source/'manifest.json').read_text())
        if source_manifest['candidate_id']!=cid:raise ValueError('CANDIDATE_ID_MISMATCH')
        files={suffix:list(source.glob('*'+suffix)) for suffix in ('.stl','.3mf','.gcode')}
        if any(len(v)!=1 for v in files.values()):raise ValueError('AMBIGUOUS_CANDIDATE_FILES')
        stl,three_mf,gcode=(files[s][0] for s in ('.stl','.3mf','.gcode'))
        hashes={p.name:sha256(p) for p in source.iterdir() if p.is_file()}
        for name,expected in source_manifest['file_sha256'].items():
            if hashes.get(name)!=expected:raise ValueError('SOURCE_MANIFEST_HASH_MISMATCH: '+name)
        with zipfile.ZipFile(three_mf) as z:
            if z.testzip() is not None:raise ValueError('SOURCE_3MF_CRC_FAILED')
            project=json.loads(z.read('Metadata/project_settings.config'))
            gcodes=[n for n in z.namelist() if n.endswith('.gcode')]
            if len(gcodes)!=1 or hashlib.sha256(z.read(gcodes[0])).hexdigest()!=sha256(gcode):
                raise ValueError('SOURCE_GCODE_PAYLOAD_MISMATCH')
        if project['printer_model']!='Bambu Lab X1 Carbon' or project['nozzle_diameter']!=['0.2']:
            raise ValueError('EXPECTED_X1_CARBON_02_MM')
        if source_manifest['box']['wall_thickness_mm']!=2 or source_manifest['box']['bottom_thickness_mm']!=2:
            raise ValueError('EXPECTED_TWO_MM_BOX')
        plans.append((cid,source,source_manifest,stl,three_mf,gcode,hashes))
    output.mkdir(parents=True)
    summary=[]
    for cid,source,original,stl,three_mf,gcode,source_hashes in plans:
        folder=output/f'candidate_{cid}';folder.mkdir()
        mapping={f'candidate_{cid}.stl':stl,
                 f'candidate_{cid}_X1C_02mm_support_ON.3mf':three_mf,
                 f'candidate_{cid}_X1C_02mm_support_ON.gcode':gcode}
        for name in ['print_orientation.png','native_bambu_preview.png','print_transform.json',
                     'casting_restore_transform.json','slicer_feature_preservation.csv','original_slice_audit.json']:
            mapping[name]=source/name
        mapping['source_manifest.json']=source/'manifest.json'
        for name,path in mapping.items():
            shutil.copyfile(path,folder/name);(folder/name).chmod(0o644)
            if sha256(folder/name)!=sha256(path):raise ValueError('COPY_HASH_MISMATCH')
        instructions=(f'Candidate {cid} — X1 Carbon / 0.2 mm / ABS\n\n'
            f'1. 请先将本 ZIP 完整解压到普通目录，例如 C:\\vascular\\candidate_{cid}。\n'
            f'2. 在 Bambu Studio 中打开 candidate_{cid}_X1C_02mm_support_ON.3mf 查看已有切片。\n'
            f'3. candidate_{cid}.stl 用于单独查看几何；STL 按 mm 解释。\n'
            '4. 3MF 本身是标准压缩容器，但使用时请直接在 Bambu Studio 打开，不要再改成 ZIP 或解压。\n\n'
            '本轮仅重新打包：侧壁和底板均为 2 mm，血管、端口、姿态、切片和 G-code 逐字节保持不变。\n'
            '工艺为 X1 Carbon、0.2 mm 喷嘴、0.10mm Standard Vessel、Bambu ABS Vessel。\n'
            f'预计时间 {original["print_time_seconds"]/60:.2f} 分钟，耗材 {original["total_filament_g"]:.2f} g。\n'
            f'支撑清理通路代理受阻区域 {original["support_removal_proxy_risk_count"]} 个，尚未进行实物拆支撑验证。\n\n'
            'source_manifest.json 与 original_slice_audit.json 保留原始来源记录，其中旧文件名指原始切片；\n'
            '本包的新文件名、来源映射和哈希以 manifest.json、SHA256SUMS.txt 为准。\n')
        (folder/'README.txt').write_text(instructions,encoding='utf-8-sig')
        manifest=dict(candidate_id=cid,operation='REPACKAGE_ONLY',geometry_modified=False,resliced=False,
            actual_project_settings=original['actual_project_settings'],box=original['box'],
            source_directory=str(source),source_files={name:str(p) for name,p in mapping.items()},
            print_time_seconds=original['print_time_seconds'],total_filament_g=original['total_filament_g'],
            support_removal_proxy_risk_count=original['support_removal_proxy_risk_count'],real_support_removal_certified=False,
            file_sha256={p.name:sha256(p) for p in folder.iterdir() if p.is_file()})
        write_json(folder/'manifest.json',manifest)
        sums={p.name:sha256(p) for p in folder.iterdir() if p.is_file()}
        (folder/'SHA256SUMS.txt').write_text(''.join(f'{h}  {n}\n' for n,h in sorted(sums.items())),encoding='ascii')
        target=output/f'candidate_{cid}_X1C02_2mm_repacked.zip';temporary=target.with_suffix('.zip.tmp')
        # DOS metadata, ASCII short paths, ordinary deflate, no links/ZIP64/encryption.
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,allowZip64=False) as z:
            for p in sorted(folder.iterdir()):
                name=folder.name+'/'+p.name
                if not name.isascii() or len(name)>120 or p.is_symlink():raise ValueError('NONPORTABLE_ZIP_MEMBER')
                info=zipfile.ZipInfo(name);info.create_system=0;info.external_attr=0x20
                info.compress_type=zipfile.ZIP_DEFLATED
                z.writestr(info,p.read_bytes(),compresslevel=6)
        with zipfile.ZipFile(temporary) as z:
            if z.testzip() is not None:raise ValueError('NEW_ZIP_CRC_FAILED')
            for p in folder.iterdir():
                if hashlib.sha256(z.read(folder.name+'/'+p.name)).hexdigest()!=sha256(p):
                    raise ValueError('NEW_ZIP_PAYLOAD_MISMATCH')
        os.replace(temporary,target)
        if any(sha256(source/name)!=h for name,h in source_hashes.items()):raise ValueError('SOURCE_CHANGED')
        summary.append(dict(candidate_id=cid,archive=str(target),bytes=target.stat().st_size,
            sha256=sha256(target),file_count=len(list(folder.iterdir())),ordinary_file=not target.is_symlink(),
            crc_and_payload_checks_passed=True,source_files_unchanged=True))
    write_json(output/'repack_summary.json',dict(status='PASS',operation='REPACKAGE_ONLY',candidates=summary))
    (output/'SHA256SUMS.txt').write_text(''.join(f"{v['sha256']}  {Path(v['archive']).name}\n" for v in summary),encoding='ascii')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,default=ROOT/'deliverables')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(json.dumps(repack(a.source_root,a.output),ensure_ascii=False,indent=2))
