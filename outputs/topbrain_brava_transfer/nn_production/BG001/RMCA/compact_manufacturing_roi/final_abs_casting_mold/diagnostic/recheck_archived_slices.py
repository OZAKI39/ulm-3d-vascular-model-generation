"""Recheck unchanged real Bambu archives after a documented coordinate audit fix."""
import sys,json,shutil,datetime
from pathlib import Path
import numpy as np
import trimesh
ROOT=Path('/home/lzy/projects/ulm_3D_vascular');sys.path.insert(0,str(ROOT))
from vascular_processing import abs_casting_mold as mold,casting_mold_qc as qc,support_removal_qc as access,bambu_casting_mold as bambu,sacrificial_fixture as old
from vascular_processing import casting_mold_review as review
from vascular_processing.sacrificial_fixture_review import save_json
out=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold'
backup=out/'diagnostic/initial_coordinate_audit';backup.mkdir(exist_ok=True)
for name in ('abs_casting_mold_summary.json','abs_casting_mold_report.md','validation_results.json','failure_summary.md'):
    path=out/name
    if path.exists() and not (backup/name).exists():shutil.copy2(path,backup/name)
for name in ('failure_evidence.png','15_bambu_slice_preview.png'):
    path=out/'QC'/name
    if path.exists() and not (backup/name).exists():shutil.copy2(path,backup/name)
s=json.loads((out/'abs_casting_mold_summary.json').read_text());cfg=mold.load_config(out/'resolved_config.yaml');inputs=mold.load_inputs(cfg)
box=mold.build_intact_box(inputs,cfg);union=trimesh.load_mesh(s['final_stl']);partition=qc.partition_vascular(union,inputs,box,cfg)
pdms=qc.ligament_qc(union,inputs,partition,cfg);walls=qc.wall_clearance(union,partition,inputs,box,cfg);_,contacts=qc.port_overlaps(inputs,box,cfg)
top=s['orientation']['top_five'];rowmap={r['candidate_id']:r for r in top};rowmap[0]=s['baseline']['orientation']
records=s['bambu']['results'];evidence=[]
for record in records:
    archive=Path(record['archive']);before=old.sha256(archive)
    evidence.append(dict(archive=str(archive),sha256=before,before_feature_pass=record['feature_preservation']['passed']))
    cached=json.loads((archive.parent/'casting_slice_audit.json').read_text())
    if cached.get('gcode_coordinate_frame',{}).get('offset_added_mm')==[0,2,0]:
        record.update(cached)
    else:
        bambu.audit_slice(record,cfg,rowmap[record['orientation_id']],union,partition,box)
    assert old.sha256(archive)==before
    evidence[-1]['after_feature_pass']=record['feature_preservation']['passed']
    print(record['orientation_id'],record['support_mode'],record['feature_preservation']['status'],record['support_quantity']['estimated_internal_support_mass_g'],flush=True)
for result_path in (out/'bambu').glob('*/slice_results.json'):
    r=json.loads(result_path.read_text());r['results']=[record for record in records if Path(record['archive']).parent.parent==result_path.parent];save_json(result_path,r)
bambu.write_slice_comparison(records,out/'tables/bambu_slice_comparison.csv')
rowmap={r['candidate_id']:r for r in top};eligible=[]
for record in records:
    row=rowmap.get(record['orientation_id'])
    if row and row['passed'] and record['support_mode']=='ON' and record['status']=='BAMBU_SLICED' and record['feature_preservation']['passed']:
        record['final_selection_rank']=row['lexicographic_rank'][:4]+[record['support_quantity']['total_support_extrusion_filament_mm']]+row['lexicographic_rank'][5:];eligible.append(record)
assert eligible,'No real sliced candidate passed feature checks'
chosen=min(eligible,key=lambda r:r['final_selection_rank']);selected=rowmap[chosen['orientation_id']]
regions=access.make_regions(union,partition,cfg);threshold=json.loads((out/'bambu/overhang_threshold.json').read_text());threshold['reference']='https://github.com/bambulab/BambuStudio/blob/master/src/libslic3r/PrintConfig.cpp';save_json(out/'bambu/overhang_threshold.json',threshold)
support=access.assess_regions(union,box,regions,np.array(selected['transform_4x4'])[:3,:3],threshold,cfg)
selected['support']=support;s['support']={k:v for k,v in support.items() if k!='risk_face_ids'};s['support']['real_support_removal_certified']=False
assert not support['trapped_risk_count']
s['orientation']['selected']={k:v for k,v in selected.items() if k!='support'}
s['bambu'].update(selected=chosen,results=records,status='BAMBU_SLICED',print_time_seconds=chosen['estimated_print_time_seconds'],filament_usage_g=chosen['filament_used_g'])
s['failures']=[f for f in s['failures'] if f not in ('SLICER_DROPPED_VASCULAR_FEATURE',"TypeError: render.<locals>.save() got multiple values for argument 'lines'")]
old.write_csv(out/'tables/internal_support_accessibility.csv',support['rows']);save_json(out/'diagnostic/selected_support_access.json',support)
old.write_csv(out/'tables/slicer_feature_preservation.csv',chosen['feature_preservation']['rows'])
shutil.copy2(selected['stl'],s['final_oriented_stl']);save_json(out/'print_transform_final.json',s['orientation']['selected'])
save_json(out/'casting_restore_transform.json',dict(transform_4x4=selected['casting_restore_transform'],casting_top_direction=[0,0,1],description='Inverse rigid transform restores the top opening to +Z; no scale change.'))
for key,extension,source in [('final_3mf','3mf','archive'),('final_gcode','gcode','gcode_file')]:
    path=out/'final'/(cfg['output']['basename']+'_ABS_casting_mold_with_vascular_core.'+extension);shutil.copy2(chosen[source],path);s[key]=str(path)
s['coordinate_audit_correction']=dict(status='RESOLVED',archived_extruder_offset_mm=[0,2,0],archive_hashes_unchanged=True,
    initial_evidence_directory=str(backup),explanation='The initial checker compared raw toolhead G-code with model coordinates. Actual archived extruder_offset=0x2 requires adding 2 mm to Y. All feature audits were recomputed without changing STL, 3MF, G-code or profiles.',
    reference='https://github.com/bambulab/BambuStudio/blob/master/src/libslic3r/GCode.cpp',archives=evidence)
s['png_files']=review.render(out,inputs,box,union,contacts,partition,pdms,walls,top,selected,s,cfg)
protected=json.loads((out/'protected_before.json').read_text());protection=old.verify_snapshot(protected);assert protection['all_unchanged'];save_json(out/'protected_source_hashes.json',protection)
s.update(protected_file_count=len(protected),hash_mismatches=protection['changed'],protected_sources_unchanged=True,status='READY_FOR_HUMAN_REVIEW' if not s['failures'] else 'NEEDS_ADJUSTMENT')
# Drop the bulky ndarray support payload added to the shared selected top row.
selected.pop('support')
save_json(out/'abs_casting_mold_summary.json',s);save_json(out/'validation_results.json',dict(status=s['status'],failures=s['failures'],warnings=s['warnings'],protected_sources_unchanged=True,human_review_pending=True))
review.write_report(out,s)
if not s['failures']:
    for path in (out/'failure_summary.md',out/'QC/failure_evidence.png'):
        if path.exists():path.unlink() # Copies retained in initial_coordinate_audit.
with (out/'run.log').open('a') as f:f.write(f'{datetime.datetime.now().isoformat()} Corrected actual extruder-offset coordinate audit; no reslicing or geometry changes. Final status {s["status"]}\n')
shutil.copy2(__file__,out/'diagnostic/recheck_archived_slices.py')
print(s['status'],selected['candidate_id'],s['bambu']['print_time_seconds'],s['bambu']['filament_usage_g'],flush=True)
