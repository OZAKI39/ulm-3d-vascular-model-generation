"""Audit the actual regenerated models and slicer packages against the 3 mm run."""
import hashlib
import json
from pathlib import Path
import sys
import zipfile
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from vascular_processing import bambu_manufacturing as bambu

AUDIT=Path(__file__).resolve().parent
BASE=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi'
OLD=BASE/'final_abs_casting_mold'
NEW=BASE/'final_abs_casting_mold_2mm_X1C_02_20260930_run02'

def read(path):return json.loads(path.read_text())
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

old=read(OLD/'abs_casting_mold_summary.json');new=read(NEW/'abs_casting_mold_summary.json')
checks={};details={}
def check(name,value):
    checks[name]=bool(value)
    if not value:raise AssertionError(name)

check('run_no_failures',not new['failures'])
check('current_source_files_unchanged',new['protected_sources_unchanged'])
old_hashes=read(AUDIT/'original_files_sha256.json')
check('all_230_original_deliverable_files_unchanged',all(sha(ROOT/p)==h for p,h in old_hashes.items()))
details['original_deliverable_file_count']=len(old_hashes)
for name in ['s1-2_swc_roi_generate_human.py','s1-3_swc_roi_generate_MeVO.py']:
    original=read(BASE/'final_abs_casting_mold_2mm_X1C_02_20260930/protected_before.json')
    check(name+'_unchanged',sha(ROOT/name)==original[str(ROOT/name)])
box=new['box']
check('four_side_walls_and_bottom_3_to_2_mm',box['wall_thickness_mm']==box['bottom_thickness_mm']==2 and
      old['box']['wall_thickness_mm']==old['box']['bottom_thickness_mm']==3)
check('outer_bounds_unchanged',np.array_equal(box['outer_bounds_mm'],old['box']['outer_bounds_mm']))
check('cavity_expanded_by_2_2_1_mm',np.allclose(np.subtract(box['inner_dimensions_mm'],old['box']['inner_dimensions_mm']),[2,2,1],atol=1e-9))
check('open_top_intact_sidewalls_and_bottom',all(box[k] for k in ['top_fully_open','four_side_walls_closed','bottom_closed']))
check('watertight_manifold_single_component',new['union']['passed'] and new['union']['connected_components']==1)
check('accepted_vascular_core_hash_unchanged',new['input']['accepted_core_hash']==old['input']['accepted_core_hash']==sha(Path(new['input']['accepted_core_path'])))
check('zero_vascular_material_loss',new['union']['core_material_loss_mm3']==0)
check('four_ports_full_wall_connection',len(new['ports'])==4 and all(r['status']=='PASS' and abs(r['intersection_axial_span_mm']-2)<1e-5 for r in new['ports']))
check('frozen_interfaces_pass',all(r['status']=='PASS' for r in new['interfaces']))
check('O3_gap_preserved',abs(new['pdms']['o3_gap_mm']-old['pdms']['o3_gap_mm'])<1e-5)
check('profile_sources_unchanged',all(sha(Path(r['path']))==r['sha256'] for r in new['bambu_discovery']['profile_sources']))
rows={r['candidate_id']:r for r in new['orientation']['top_five']};rows[0]=new['baseline']['orientation']
old_rows={r['candidate_id']:r for r in old['orientation']['top_five']};old_rows[0]=old['baseline']['orientation']
candidates=[]
for cid,slug in [(0,'upright'),(15,'side')]:
    check(f'candidate_{cid}_original_pose_preserved',np.allclose(rows[cid]['transform_4x4'],old_rows[cid]['transform_4x4'],atol=1e-8))
    folder=NEW/'candidate_0_and_15'/f'candidate_{cid:02}_{slug}'
    manifest=read(folder/'manifest.json')
    archive=next(folder.glob('*.3mf'));gcode=next(folder.glob('*.gcode'));stl=next(folder.glob('*.stl'))
    parsed=bambu.parse_sliced_archive(archive);settings=parsed['actual_project_settings']
    check(f'candidate_{cid}_X1_Carbon_02_ABS',settings['printer_model']=='Bambu Lab X1 Carbon' and
          settings['nozzle_diameter']==['0.2'] and settings['filament_type']==['ABS'])
    check(f'candidate_{cid}_active_custom_process',settings['print_settings_id']=='0.10mm Standard Vessel' and
          settings['filament_settings_id']==['Bambu ABS Vessel'] and float(settings['layer_height'])==.1 and
          settings['nozzle_temperature'][0]=='265' and settings['nozzle_temperature_initial_layer'][0]=='270')
    check(f'candidate_{cid}_real_sliced_support_on',parsed['slice_success'] and settings['enable_support']=='1')
    check(f'candidate_{cid}_sliced_geometry_preserved',bambu.sliced_geometry_qc(archive,stl)['rotation_scale_position_preserved'])
    record=read(folder/'original_slice_audit.json')
    check(f'candidate_{cid}_real_features_preserved',record['feature_preservation']['passed'])
    with zipfile.ZipFile(archive) as z:
        payload=z.read(parsed['gcode_files'][0]);text=payload.decode(errors='replace')
        check(f'candidate_{cid}_external_gcode_identical',hashlib.sha256(payload).hexdigest()==sha(gcode))
        check(f'candidate_{cid}_real_Gcode_temperature',any('S270' in line and line.startswith(('M104','M109')) for line in text.splitlines()) and
              any('S265' in line and line.startswith(('M104','M109')) for line in text.splitlines()))
    with zipfile.ZipFile(folder.with_suffix('.zip')) as z:check(f'candidate_{cid}_ZIP_CRC',z.testzip() is None)
    candidates.append(dict(candidate_id=cid,print_time_minutes=parsed['estimated_print_time_seconds']/60,
        total_filament_g=parsed['filament_used_g'],support_removal_proxy_risk_count=manifest['support_removal_proxy_risk_count'],
        feature_count=len(record['feature_preservation']['rows']),sha256={p.name:sha(p) for p in [stl,archive,gcode]}))
details.update(box=box,core_material_loss_mm3=new['union']['core_material_loss_mm3'],o3_gap_mm=new['pdms']['o3_gap_mm'],
    candidates=candidates,source_summary=str(NEW/'abs_casting_mold_summary.json'))
result=dict(status='PASS',checks=checks,passed=len(checks),details=details)
(AUDIT/'real_regeneration_validation.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(dict(status='PASS',passed=len(checks),candidates=candidates),indent=2,ensure_ascii=False))
