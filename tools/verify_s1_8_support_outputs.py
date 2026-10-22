#!/usr/bin/env python3
"""Verify delivered files; NEEDS_ADJUSTMENT may still have intact artifacts."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pyvista as pv
import trimesh
from vascular_processing.support_inputs import ROOT,save,verify,digest
from vascular_processing.support_artifacts import verify_3mf,write_manifests

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_support_planner')
    out=parser.parse_args().output.resolve()
    summary=json.loads((out/'vascular_support_summary.json').read_text())
    protection=verify(json.loads((out/'protected_before.json').read_text()))
    meshes=[];archives=[];slices=[];failures=[]
    files=list((out/'geometry').glob('*.stl'))+[out/'coupon/support_breakaway_coupon.stl']
    for folder in sorted((out/'custom_candidates').glob('custom_*')):
        for name in ['custom_support.stl','final_mold_with_breakaway_support.stl','support_only_with_approved_bottom.stl']:
            files.append(folder/name)
        for mode in ['bambu_assembly','bambu_support_only']:
            path=folder/mode/'slice_result.json';result=json.loads(path.read_text())
            metrics=json.loads((folder/mode/'custom_slice_metrics.json').read_text())
            current=folder/('final_mold_with_breakaway_support.stl' if mode=='bambu_assembly' else 'support_only_with_approved_bottom.stl')
            geometry=result.get('production_geometry',{})
            matched=geometry.get('expected_sha256')==digest(current)
            ok=result.get('slicer_success') and geometry.get('passed') and matched and metrics.get('auto_support_off') and metrics.get('native_support_toolpath_count')==0 and metrics['features']['passed']
            slices.append(dict(candidate=folder.name,mode=mode,passed=bool(ok),source_hash_matches=matched,
                features_preserved=metrics['features']['passed'],auto_support_off=metrics['auto_support_off'],
                native_support_toolpaths=metrics['native_support_toolpath_count'],critical_warnings=metrics.get('critical_warnings',[])))
            if not ok:failures.append(str(path))
            archives.append(verify_3mf(folder/mode/'sliced.3mf'))
    for path in files:
        mesh=trimesh.load_mesh(path)
        meshes.append(dict(path=str(path.relative_to(out)),is_volume=bool(mesh.is_volume),faces=len(mesh.faces),volume_mm3=float(mesh.volume),sha256=digest(path)))
        if not mesh.is_volume:failures.append(str(path))
    for path in [out/'custom_support_auto_support_off.3mf',out/'custom_support_only_auto_support_off.3mf',out/'coupon/bambu/sliced.3mf']:
        archives.append(verify_3mf(path))
    for row in archives:
        if not row['zip_crc_passed'] or not row['gcode_members']:failures.append(row['path'])
    coupon=json.loads((out/'coupon/slice_validation.json').read_text())
    if not (coupon['slicer_success'] and coupon['auto_support_off'] and coupon['native_support_toolpath_count']==0 and coupon['all_coupon_features_preserved']):failures.append('coupon slice features')
    if coupon['production_geometry']['expected_sha256']!=digest(out/'coupon/support_breakaway_coupon.stl'):failures.append('coupon stale geometry')
    vtp={p.name:dict(points=pv.read(p).n_points,cells=pv.read(p).n_cells) for p in (out/'geometry').glob('*.vtp')}
    images=json.loads((out/'visual_qc_manifest.json').read_text())['main_images']
    if len(images)!=18 or not all(Path(p).is_file() for p in images):failures.append('18 visual QC images')
    result=dict(artifact_integrity_passed=not failures,failures=failures,source_protection=protection,
        scale_factor=summary['scale_factor'],production_status=summary['status'],
        production_ready=summary['production_selected'] is not None,hard_gate_failures=summary['selected']['hard_gate_failures'],
        meshes=meshes,archives=archives,real_slices=slices,coupon_slice=dict(succeeded=coupon['slicer_success'],auto_support_off=coupon['auto_support_off'],all_features_preserved=coupon['all_coupon_features_preserved']),
        vtp=vtp,qc_image_count=len(images),printer_upload=False,printer_start=False)
    save(out/'verification.json',result);write_manifests(out,summary)
    print(json.dumps({k:result[k] for k in ['artifact_integrity_passed','failures','source_protection','production_status','production_ready','qc_image_count']},ensure_ascii=False))
    return 0 if not failures else 1

if __name__=='__main__':raise SystemExit(main())
