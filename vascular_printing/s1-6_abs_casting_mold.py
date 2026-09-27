#!/usr/bin/env python3
"""Frozen vascular core + intact five-wall casting box, tests and local slice evidence."""
import argparse
import json
import logging
from pathlib import Path
import shutil
import warnings
import numpy as np
import trimesh
import yaml
import cadquery as cq
from vascular_processing import abs_casting_mold as mold
from vascular_processing import casting_mold_qc as qc
from vascular_processing import support_removal_qc as access
from vascular_processing import bambu_casting_mold as bambu
from vascular_processing import sacrificial_fixture as old
from vascular_processing.sacrificial_print_frame import box_mesh
from vascular_processing.sacrificial_fixture_review import save_json,polydata


def small_row(row):return {k:v for k,v in row.items() if k!='support'}
def small_access(value):return {k:v for k,v in value.items() if k!='risk_face_ids'}


def run(config_path,output_root=None,skip_slicing=False):
    cfg=mold.load_config(config_path);inputs=mold.load_inputs(cfg)
    out=Path(output_root).resolve() if output_root else inputs['base']/cfg['output']['directory']
    if out.exists():raise FileExistsError('Refusing to overwrite any existing review directory: '+str(out))
    if out.is_relative_to(inputs['folder']):raise ValueError('ACCEPTED_RESULTS_ARE_FROZEN')
    out.mkdir(parents=True)
    for folder in ('reference','final','diagnostic','tables','QC','bambu','orientations'):(out/folder).mkdir()
    logger=logging.getLogger('abs_casting_mold');logger.setLevel(logging.INFO)
    handler=logging.FileHandler(out/'run.log',encoding='utf-8');handler.setFormatter(logging.Formatter('%(asctime)s %(message)s'));logger.addHandler(handler)
    protected=mold.snapshot(inputs,out);save_json(out/'protected_before.json',protected)
    (out/'resolved_config.yaml').write_text(yaml.safe_dump(cfg,allow_unicode=True,sort_keys=False),encoding='utf-8')
    summary=dict(status='NEEDS_ADJUSTMENT',production_type='OPEN_TOP_FIVE_WALL_CASTING_BOX',failures=[],warnings=[],
        input=dict(accepted_core_path=str(inputs['core_path']),accepted_core_hash=old.sha256(inputs['core_path']),
            summary_authority=str(inputs['summary_path']),report_authority=str(inputs['report_path']),status='ALL_PORT_ATTACHMENTS_ACCEPTED'),
        old_frames_preserved_as_historical_exploration=True,upstream_recomputed=False,vascular_geometry_modified=False,printer_job_sent=False)
    box=None;union=None;partition=None;pdms=None;wall_rows=[];selected=None;top=[];all_rows=[];contacts={};records=[];discovery={}
    try:
        logger.info('Read immutable accepted vascular core; generate intact box from audited dimensions')
        shutil.copy2(inputs['core_path'],out/'reference'/inputs['core_path'].name)
        box=mold.build_intact_box(inputs,cfg);box_path=out/'reference/BG001_RMCA_intact_open_top_box.stl'
        old.export_checked_stl(box['mesh'],box_path,cfg)
        cq.exporters.export(box['solid'],str(out/'reference/BG001_RMCA_intact_open_top_box.step'))
        summary['box']=qc.intact_box_qc(box,cfg);summary['box'].update(inner_bounds_mm=box['inner'],outer_bounds_mm=box['outer'])
        if not all(summary['box'][k] for k in ('four_side_walls_closed','bottom_closed','top_fully_open','no_assembly_clearance_holes')):
            raise ValueError('FIVE_WALL_BOX_STRUCTURE_FAILED')
        port_rows,contacts=qc.port_overlaps(inputs,box,cfg);summary['ports']=port_rows;old.write_csv(out/'tables/port_wall_union.csv',port_rows)
        failed=[r['status'] for r in port_rows if r['status']!='PASS']
        if failed:raise ValueError('; '.join(failed))
        union,boolean_info=mold.combine(inputs['core'],box['mesh'],cfg)
        basename=cfg['output']['basename'];final_path=out/'final'/(basename+'_ABS_casting_mold_with_vascular_core.stl')
        old.export_checked_stl(union,final_path,cfg);union=trimesh.load_mesh(final_path)
        summary['final_stl']=str(final_path);summary['union']=qc.union_audit(inputs['core'],box['mesh'],union,cfg)
        summary['union']['numerical_cleanup_used']=boolean_info['boolean_cleanup_used'];save_json(out/'boolean_union_audit.json',summary['union'])
        summary['failures'].extend(summary['union']['failures'])
        interfaces=qc.interface_geometry_checks(inputs['core'],union,inputs);summary['interfaces']=interfaces;old.write_csv(out/'tables/frozen_interface_qc.csv',interfaces)
        if any(r['status']!='PASS' for r in interfaces):summary['failures'].append('ACCEPTED_INTERFACE_GEOMETRY_CHANGED')
        logger.info('Union complete: components=%s core loss=%g box loss=%g',summary['union']['connected_components'],summary['union']['core_material_loss_mm3'],summary['union']['box_material_loss_mm3'])
        partition=qc.partition_vascular(union,inputs,box,cfg)
        np.savez_compressed(out/'diagnostic/vascular_surface_partition.npz',face_ids=partition['face_ids'],labels=partition['labels'],interior=partition['interior'])
        pdms=qc.ligament_qc(union,inputs,partition,cfg);wall_rows=qc.wall_clearance(union,partition,inputs,box,cfg)
        summary['pdms']={k:v for k,v in pdms.items() if k!='all_pairs'};summary['pdms'].update(minimum_core_wall_mm=wall_rows[0]['distance_mm'],
            minimum_core_wall_outside_legal_band_mm=wall_rows[0]['distance_mm'],wall_crossing_exemption_mm=cfg['pdms']['wall_crossing_exemption_mm'],
            minimum_native_core_wall_mm=min(r['distance_mm'] for r in wall_rows if not r['vascular_branch'].startswith(('I','O'))))
        old.write_csv(out/'tables/pdms_ligament_top10.csv',pdms['top10']);old.write_csv(out/'tables/pdms_ligament_all_pairs.csv',pdms['all_pairs'])
        old.write_csv(out/'tables/core_to_inner_wall_clearance.csv',wall_rows[:10]);old.write_csv(out/'tables/core_to_inner_wall_clearance_all.csv',wall_rows)
        if pdms['o3_status']!='PASS':summary['failures'].append('O3_PDMS_GAP_REGRESSION')
        logger.info('PDMS minimum ligament=%.9f O3 deviation=%.9g',pdms['minimum_ligament_mm'],pdms['o3_deviation_mm'])
        future=trimesh.boolean.difference([box_mesh(box['inner']),inputs['core']],engine='manifold');polydata(future).save(out/'diagnostic/future_PDMS_volume.vtp')
        summary['pdms']['future_volume_mm3']=abs(float(future.volume));summary['pdms']['diagnostic_only']=True
        discovery=bambu.discover(cfg,out/'bambu');summary['bambu_discovery']={k:v for k,v in discovery.items() if k!='flattened'}
        threshold=access.overhang_threshold(discovery,cfg);save_json(out/'bambu/overhang_threshold.json',threshold)
        if threshold['source']=='OVERHANG_THRESHOLD_HEURISTIC':summary['warnings'].append('OVERHANG_THRESHOLD_HEURISTIC')
        build=discovery.get('build_volume_mm') or cfg['printer']['provisional_build_volume_mm']
        regions=access.make_regions(union,partition,cfg);save_json(out/'diagnostic/support_regions.json',regions)
        cache={};manager=trimesh.collision.CollisionManager();manager.add_object('mold',union)
        baseline_access=access.assess_regions(union,box,regions,np.eye(3),threshold,cfg,cache,manager)
        baseline=access.orientation_row(0,'TOP_OPENING_UP_yaw0',np.eye(3),union,box,build,baseline_access,cfg);baseline['rank']=0
        basefile=out/'orientations/baseline_upright.stl';bm=union.copy();bm.apply_transform(baseline['transform_4x4']);old.export_checked_stl(bm,basefile,cfg);baseline['stl']=str(basefile)
        save_json(out/'diagnostic/upright_support_access.json',baseline_access)
        summary['baseline']=dict(orientation=small_row(baseline),support=small_access(baseline_access))
        logger.info('Mandatory upright baseline: risk regions=%s trapped proxy=%s. Slice OFF before ON.',baseline_access['candidate_regions'],baseline_access['trapped_risk_count'])
        if not skip_slicing:
            result=bambu.slice_one(discovery,cfg,baseline,out/'bambu/baseline_upright',logger,union,partition,box);records.extend(result['results'])
        # Search is conditional on the actual baseline risk and slicer evidence.
        needs_search=baseline_access['trapped_risk_count']>0 or baseline_access['internal_support_risk_area_mm2']>cfg['orientation']['baseline_search_trigger_area_mm2']
        rotations=access.rotations(cfg) if needs_search else [('TOP_OPENING_UP_yaw0',np.eye(3))]
        support_cache={tuple(np.round(np.eye(3)[2],8)):baseline_access}
        for index,(name,rotation) in enumerate(rotations):
            key=tuple(np.round(rotation[2],8))
            if key not in support_cache:support_cache[key]=access.assess_regions(union,box,regions,rotation,threshold,cfg,cache,manager)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',message='Gimbal lock detected.*')
                row=access.orientation_row(index,name,rotation,union,box,build,support_cache[key],cfg)
            all_rows.append(row)
        placed=sorted([r for r in all_rows if r['geometry_placement_pass']],key=lambda r:r['lexicographic_rank'])
        top=placed[:cfg['orientation']['top_count']]
        if not top:raise ValueError('NO_BUILD_VOLUME_AND_BED_STABILITY_VALID_ORIENTATION')
        for rank,row in enumerate(top,1):
            row['rank']=rank;path=out/'orientations'/f'top_{rank:02}.stl';mesh=union.copy();mesh.apply_transform(row['transform_4x4']);old.export_checked_stl(mesh,path,cfg);row['stl']=str(path)
        save_json(out/'orientations/all_orientations.json',[small_row(r) for r in all_rows])
        fields=['candidate_id','method','rank','fits_build_volume','stable','trapped_risk_count','internal_support_risk_area_mm2','top_accessible_count','probe_accessible_count','external_support_area_mm2','downward_overhang_area_mm2','z_height_mm','bed_footprint_mm2','passed']
        old.write_csv(out/'tables/orientation_candidates.csv',[{k:r.get(k) for k in fields} for r in all_rows])
        logger.info('Orientation search: %s candidates; Top5 ids=%s',len(all_rows),[r['candidate_id'] for r in top])
        for row in top:
            if row['candidate_id']==0:continue # Baseline was actually sliced first.
            if not skip_slicing:
                result=bambu.slice_one(discovery,cfg,row,out/'bambu'/f"rank_{row['rank']:02}",logger,union,partition,box);records.extend(result['results'])
        bambu.write_slice_comparison(records,out/'tables/bambu_slice_comparison.csv')
        rowmap={r['candidate_id']:r for r in top};eligible=[]
        for record in records:
            row=rowmap.get(record['orientation_id'])
            if row and row['passed'] and record['support_mode']=='ON' and record['status']=='BAMBU_SLICED' and record['feature_preservation']['passed']:
                # Refine the support-amount tie using actual filament extrusion;
                # geometry hard gates retain their lexicographic priority.
                record['final_selection_rank']=row['lexicographic_rank'][:4]+[record['support_quantity']['total_support_extrusion_filament_mm']]+row['lexicographic_rank'][5:]
                eligible.append(record)
        chosen_slice=min(eligible,key=lambda r:r['final_selection_rank']) if eligible else None
        selected=rowmap[chosen_slice['orientation_id']] if chosen_slice else next((r for r in top if r['passed']),top[0])
        summary['support']=small_access(selected['support']);summary['support']['real_support_removal_certified']=False
        save_json(out/'diagnostic/selected_support_access.json',selected['support'])
        old.write_csv(out/'tables/internal_support_accessibility.csv',selected['support']['rows'])
        if selected['support']['trapped_risk_count']:summary['failures'].append('TRAPPED_SUPPORT_RISK')
        summary['orientation']=dict(candidate_count=len(all_rows),search_triggered_by_baseline_risk=needs_search,
            hard_gate_pass_count=sum(r['passed'] for r in all_rows),selected=small_row(selected),top_five=[small_row(r) for r in top],
            build_volume_mm=build,build_volume_confirmed=discovery['status']=='BAMBU_PROFILE_CONFIRMED',
            ranking_order=['trapped risk','internal support risk area','top probe access','external support area','total support amount','Z height','footprint'],
            mandatory_baseline_sliced_first=True)
        oriented=out/'final'/(basename+'_ABS_casting_mold_with_vascular_core_oriented.stl');shutil.copy2(selected['stl'],oriented);summary['final_oriented_stl']=str(oriented)
        save_json(out/'print_transform_final.json',small_row(selected));save_json(out/'casting_restore_transform.json',dict(transform_4x4=selected['casting_restore_transform'],casting_top_direction=[0,0,1],description='Inverse rigid transform restores the top opening to +Z; no scale change.'))
        if chosen_slice:
            path=out/'final'/(basename+'_ABS_casting_mold_with_vascular_core.3mf');shutil.copy2(chosen_slice['archive'],path);summary['final_3mf']=str(path)
            path=out/'final'/(basename+'_ABS_casting_mold_with_vascular_core.gcode');shutil.copy2(chosen_slice['gcode_file'],path);summary['final_gcode']=str(path)
            old.write_csv(out/'tables/slicer_feature_preservation.csv',chosen_slice['feature_preservation']['rows'])
        elif any(r['status']=='BAMBU_SLICED' and not r.get('feature_preservation',{}).get('passed',False) for r in records):
            summary['failures'].append('SLICER_DROPPED_VASCULAR_FEATURE')
        summary['bambu']=dict(binary=discovery.get('binary'),version=discovery.get('version'),printer=discovery.get('model'),
            nozzle_mm=discovery.get('nozzle_diameter_mm'),filament=discovery.get('active_gui_presets',{}).get('filaments'),
            process=discovery.get('active_gui_presets',{}).get('process'),results=records,selected=chosen_slice,
            support_geometry_available=False,support_geometry_status='BAMBU_SUPPORT_GEOMETRY_UNAVAILABLE',
            support_off_result=[dict(candidate_id=r['orientation_id'],status=r['status']) for r in records if r['support_mode']=='OFF'],
            support_on_result=[dict(candidate_id=r['orientation_id'],status=r['status']) for r in records if r['support_mode']=='ON'],
            print_time_seconds=chosen_slice.get('estimated_print_time_seconds') if chosen_slice else None,
            filament_usage_g=chosen_slice.get('filament_used_g') if chosen_slice else None,
            status='BAMBU_SLICED' if chosen_slice else discovery['status'] if discovery['status']!='BAMBU_PROFILE_CONFIRMED' else 'BAMBU_VALIDATION_PENDING')
        summary['warnings'].append('BAMBU_SUPPORT_GEOMETRY_UNAVAILABLE: actual toolpaths inspected, no independent support solid; top removal remains a geometry proxy.')
        if any(r.get('warnings') for r in records):summary['warnings'].append('BAMBU_METADATA_WARNINGS: anomalous first-layer time metadata retained but not used.')
        from vascular_processing.casting_mold_review import render
        summary['png_files']=render(out,inputs,box,union,contacts,partition,pdms,wall_rows,top,selected,summary,cfg)
    except Exception as exc:
        logger.exception('Casting mold failed');summary['failures'].append(type(exc).__name__+': '+str(exc))
        from vascular_processing.casting_mold_review import failure_image
        failure_image(out,inputs,box,union,str(exc),cfg)
    protection=old.verify_snapshot(protected);save_json(out/'protected_source_hashes.json',protection)
    summary.update(protected_file_count=len(protected),hash_mismatches=protection['changed'],protected_sources_unchanged=protection['all_unchanged'])
    if not protection['all_unchanged']:summary['failures'].append('PROTECTED_SOURCE_CHANGED')
    if not summary['failures']:summary['status']='READY_FOR_HUMAN_REVIEW' if summary.get('bambu',{}).get('selected') else 'READY_FOR_HUMAN_REVIEW_WITH_BAMBU_PENDING'
    save_json(out/'abs_casting_mold_summary.json',summary);save_json(out/'validation_results.json',dict(status=summary['status'],failures=summary['failures'],warnings=summary['warnings'],protected_sources_unchanged=protection['all_unchanged'],human_review_pending=True))
    from vascular_processing.casting_mold_review import write_report
    write_report(out,summary)
    if summary['failures']:(out/'failure_summary.md').write_text('# 当前模型需要调整\n\n'+'\n'.join('- '+r for r in summary['failures'])+'\n\n本次几何、切片、CSV 和图像均保留；没有修改冻结血管或既有结果。\n',encoding='utf-8')
    logger.info('Final status: %s',summary['status']);handler.close();logger.removeHandler(handler)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',type=Path,default=mold.ROOT/'config/abs_casting_mold_BG001.yaml');parser.add_argument('--output-root',type=Path)
    parser.add_argument('--skip-slicing',action='store_true',help='Explicit geometry-only diagnostic; default performs real local slicing')
    args=parser.parse_args();result=run(args.config,args.output_root,args.skip_slicing);print(result['status']);return 0 if result['status'].startswith('READY_FOR_HUMAN_REVIEW') else 2

if __name__=='__main__':raise SystemExit(main())
