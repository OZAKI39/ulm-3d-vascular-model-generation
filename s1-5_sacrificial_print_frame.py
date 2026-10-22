#!/usr/bin/env python3
"""Build and review an open temporary ABS print frame around accepted BG001 geometry."""
import argparse
import json
import logging
from pathlib import Path
import shutil
import sys
import warnings

import numpy as np
import trimesh
import yaml

from vascular_processing import sacrificial_print_frame as frame
from vascular_processing import print_frame_qc as qc
from vascular_processing import support_access_qc as access
from vascular_processing import bambu_slice_adapter as bambu
from vascular_processing import sacrificial_fixture as legacy
from vascular_processing.sacrificial_fixture_review import save_json


def public_design(design):
    return {k:v for k,v in design.items() if k not in ('mesh','ports_for_ranking','core','intrusion_mesh','contact_mesh')}


def run(config_path,output_root=None,skip_slicing=False):
    cfg=frame.load_config(config_path);inputs=frame.load_accepted(cfg)
    out=Path(output_root).resolve() if output_root else inputs['base']/cfg['output']['directory']
    if out.exists():raise FileExistsError('Refusing to overwrite print-frame review: '+str(out))
    if out.is_relative_to(inputs['folder']):raise ValueError('ACCEPTED_RESULTS_ARE_FROZEN')
    out.mkdir(parents=True)
    for name in ('candidates','final','tables','QC','bambu','orientations'):(out/name).mkdir()
    logger=logging.getLogger('sacrificial_print_frame');logger.setLevel(logging.INFO)
    handler=logging.FileHandler(out/'run.log',encoding='utf-8');logger.addHandler(handler)
    snapshot=frame.protected_snapshot(inputs,out);save_json(out/'protected_before.json',snapshot)
    (out/'resolved_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False,allow_unicode=True),encoding='utf-8')
    summary=dict(status='NEEDS_ADJUSTMENT',input_status='ALL_PORT_ATTACHMENTS_ACCEPTED',source_hash=legacy.sha256(inputs['core_path']),
        source_path=str(inputs['core_path']),selected_candidate=None,output_directory=str(out),failures=[],warnings=[],
        vascular_geometry_changed=False,ports_changed=False,upstream_recomputed=False,
        whole_abs_object_enters_separate_pdms_container=True,frame_is_pdms_container=False,
        printer_job_sent=False,geometry_support_removal_certified=False)
    basename=cfg['output']['basename'];chosen={};all_layouts=[];top=[];selected=None
    try:
        vascular_copy=out/'final'/(basename+'_vascular_core_only.stl');shutil.copy2(inputs['core_path'],vascular_copy)
        if legacy.sha256(vascular_copy)!=summary['source_hash']:raise ValueError('VASCULAR_COPY_HASH_MISMATCH')
        summary['vascular_only_file']=str(vascular_copy)
        keep=frame.keep_zone_bounds(inputs['central'].bounds,cfg)
        summary['pdms_keep_zone']=dict(basis=cfg['pdms_keep_zone']['basis'],authorization=cfg['pdms_keep_zone']['authorization'],
            bounds_mm=keep,dimensions_mm=keep[1]-keep[0],central_reference=str(inputs['central_path']))
        discovery=bambu.discover_current(cfg,out/'bambu');summary['bambu_discovery']={k:v for k,v in discovery.items() if k!='flattened'}
        summary['warnings'].extend(discovery.get('warnings',[]))
        logger.info('Bambu discovery: %s, %s, nozzle %s',discovery['status'],discovery.get('model'),discovery.get('nozzle_diameter_mm'))
        regions=access.vascular_regions(inputs,cfg);save_json(out/'vascular_regions.json',regions)
        for kind in ('OPEN_PANEL','SPARSE_FRAME'):
            layouts=[]
            for axis in (0,1):
                directory=out/'candidates'/kind/('OPEN_X' if axis==0 else 'OPEN_Y');directory.mkdir(parents=True)
                logger.info('Build %s %s',kind,directory.name)
                design=frame.build_frame(kind,axis,inputs,cfg);design['ports_for_ranking']=inputs['ports']
                checks,intrusion,contact=qc.frame_checks(design,inputs,keep,cfg)
                design.update(qc=checks,intrusion_mesh=intrusion,contact_mesh=contact)
                frame_path=directory/'frame_only.stl';design['mesh'].export(frame_path)
                checks['exported_frame_mesh']=legacy.mesh_qc(trimesh.load_mesh(frame_path),cfg)
                if not checks['exported_frame_mesh']['passed']:checks['failures'].append('FRAME_STL_ROUNDTRIP_FAILED')
                try:
                    core,union_qc=frame.combine(inputs['core'],design['mesh'],cfg)
                    core_path=directory/'vascular_with_frame.stl'
                    export_qc=legacy.export_checked_stl(core,core_path,cfg)
                    design['core']=trimesh.load_mesh(core_path);checks['union']=union_qc;checks['exported_union']=export_qc
                    design['union_file']=str(core_path)
                except Exception as exc:
                    checks['failures'].append('PORT_ONLY_FRAME_NOT_STRUCTURALLY_CONNECTED: '+str(exc));logger.exception('Candidate union failure')
                proxy,rows,paths=access.evaluate_access(design,inputs['core'],regions,cfg)
                design.update(access=proxy,access_rows=rows,access_paths=paths,frame_file=str(frame_path))
                if proxy['accessible_fraction']<cfg['support_access']['minimum_accessible_fraction']:checks['failures'].append('GEOMETRY_SUPPORT_ACCESS_RISK')
                if proxy['trapped_proxy_region_count']:checks['failures'].append('TRAPPED_SUPPORT_GEOMETRY_PROXY')
                checks['passed']=not checks['failures']
                design['lexicographic_rank']=qc.layout_rank(design,checks,proxy)
                legacy.write_csv(directory/'support_accessibility.csv',rows)
                save_json(directory/'candidate_qc.json',public_design(design))
                logger.info('%s %s volume=%.6f intrusion=%.9g access=%s/%s probe=%s/%s failures=%s',kind,design['layout'],
                    checks['frame_volume_mm3'],checks['frame_intrusion_keep_zone_mm3'],proxy['accessible_region_count'],proxy['region_count'],
                    proxy['probe_accessible_region_count'],proxy['region_count'],checks['failures'])
                layouts.append(design);all_layouts.append(design)
            feasible=[d for d in layouts if d['qc']['passed']]
            chosen[kind]=min(feasible or layouts,key=lambda d:d['lexicographic_rank'])
        comparison=[]
        for d in all_layouts:
            q=d['qc'];a=d['access'];faces=q['open_faces']
            comparison.append(dict(candidate=d['candidate'],open_faces='/'.join(d['open_faces']),selected_layout=d is chosen[d['candidate']],
                frame_volume_mm3=q['frame_volume_mm3'],frame_intrusion_keep_zone_mm3=q['frame_intrusion_keep_zone_mm3'],
                open_face_A_free_fraction=faces[0]['free_area_fraction'],open_face_B_free_fraction=faces[1]['free_area_fraction'],
                support_access_score=a['accessible_fraction'],probe_access_fraction=a['probe_accessible_fraction'],
                port_connection_count=sum(p['status']=='PASS' for p in q['port_connections']),
                connected_components=q.get('union',{}).get('connected_components'),estimated_support_risk=a['status'],
                status='PASS' if q['passed'] else 'FAIL'))
        legacy.write_csv(out/'tables/frame_candidate_comparison.csv',comparison)
        selected=chosen['SPARSE_FRAME'];q=selected['qc'];a=selected['access']
        summary['candidates']={k:public_design(v) for k,v in chosen.items()}
        summary['layout_selection']=dict(method='LEXICOGRAPHIC',order=['accessible region fraction','4 mm probe accessible fraction',
            'ports on selected open faces','keep-zone intrusion','frame volume','downward frame surface area'],
            reason='Compare OPEN_X and OPEN_Y within each family. Sparse geometry is never replaced by the panel candidate.',
            evaluated=[dict(candidate=d['candidate'],layout=d['layout'],rank=d['lexicographic_rank']) for d in all_layouts])
        summary['frame']=dict(type='SPARSE_FRAME',open_faces=selected['open_faces'],dimensions_mm=np.diff(selected['outer_bounds'],axis=0)[0],
            outer_bounds_mm=selected['outer_bounds'],volume_mm3=q['frame_volume_mm3'],connected_components=q['frame_mesh']['connected_components'],
            port_connections=q['port_connections'],material_reduction_vs_open_panel_fraction=1-q['frame_volume_mm3']/chosen['OPEN_PANEL']['qc']['frame_volume_mm3'],
            additional_ghost_void_volume_mm3=q['additional_ghost_void_volume_mm3'],unintended_mid_vessel_contact_mm3=q['unintended_mid_vessel_contact_mm3'])
        summary['pdms_keep_zone'].update(frame_intrusion_mm3=q['frame_intrusion_keep_zone_mm3'],ghost_void_volume_mm3=q['ghost_void_volume_inside_keep_zone_mm3'])
        summary['access']=dict(**a,open_face_metrics=q['open_faces'],inaccessible_regions=[r['region_id'] for r in selected['access_rows'] if r['status']!='GEOMETRY_SUPPORT_ACCESS_PASS'])
        legacy.write_csv(out/'tables/port_frame_connection.csv',q['port_connections'])
        legacy.write_csv(out/'tables/support_accessibility.csv',selected['access_rows'])
        legacy.write_csv(out/'tables/open_face_access.csv',q['open_faces'])
        save_json(out/'support_probe_paths.json',selected['access_paths'])
        if not q['passed']:
            summary['failures'].extend(q['failures']);summary['failures'].append('SPARSE_FRAME_NEEDS_ADJUSTMENT')
        else:
            summary['selected_candidate']='SPARSE_FRAME'
            final=out/'final'/(basename+'_vascular_with_sparse_frame.stl');shutil.copy2(selected['union_file'],final)
            summary['final_stl']=str(final);summary['final_mesh_qc']=legacy.mesh_qc(selected['core'],cfg)
            summary['source_volume_loss_mm3']=q['union']['source_subtracted_volume_mm3']
            interfaces=qc.interface_geometry_checks(inputs['core'],selected['core'],inputs)
            legacy.write_csv(out/'tables/frozen_interface_qc.csv',interfaces);summary['interface_qc']=interfaces
            if any(r['status']!='PASS' for r in interfaces):summary['failures'].append('ACCEPTED_INTERFACE_GEOMETRY_CHANGED')
            build=discovery['build_volume_mm'] if discovery['status']=='BAMBU_PROFILE_CONFIRMED' else cfg['printer']['provisional_build_volume_mm']
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore',message='Gimbal lock detected.*')
                orientations,top=qc.rank_orientations(selected['core'],selected['mesh'],build,a,cfg)
            fields=['candidate_id','method','fits_build_volume','stable','bed_contact_area_mm2','support_polygon_margin_mm',
                'downward_overhang_area_mm2','z_height_mm','bed_footprint_mm2','passed','status']
            legacy.write_csv(out/'tables/orientation_candidates.csv',[{k:r[k] for k in fields} for r in orientations])
            save_json(out/'orientations/all_orientations.json',orientations)
            for row in top:
                path=out/'orientations'/f"top_{row['rank']:02}.stl";mesh=selected['core'].copy();mesh.apply_transform(row['transform_4x4'])
                legacy.export_checked_stl(mesh,path,cfg);row['stl']=str(path)
            summary['orientation']=dict(candidate_count=len(orientations),feasible_count=sum(r['passed'] for r in orientations),
                top_five=top,build_volume_mm=build,build_volume_confirmed=discovery['status']=='BAMBU_PROFILE_CONFIRMED')
            if not top:summary['failures'].append('NO_STABLE_PRINT_ORIENTATION')
            else:
                slicing=bambu.slice_top_five(discovery,cfg,top,out/'bambu',logger) if not skip_slicing else dict(status='BAMBU_SKIPPED_BY_CLI_OPTION',results=[])
                summary['bambu']=slicing
                selected_slice=slicing.get('selected')
                oriented=next((r for r in top if selected_slice and r['candidate_id']==selected_slice['orientation_id']),top[0])
                summary['orientation']['selected']=oriented
                transform_path=out/'print_transform_frame.json';save_json(transform_path,oriented)
                final_oriented=out/'final'/(basename+'_vascular_with_sparse_frame_oriented.stl');shutil.copy2(oriented['stl'],final_oriented)
                summary['final_oriented_stl']=str(final_oriented)
                if selected_slice:
                    final_3mf=out/'final'/(basename+'_vascular_with_sparse_frame.3mf');shutil.copy2(selected_slice['archive'],final_3mf);summary['final_3mf']=str(final_3mf)
                    summary['warnings'].append('BAMBU_SUPPORT_GEOMETRY_NOT_AVAILABLE: actual support toolpaths do not certify removable solid supports.')
                else:summary['warnings'].append('Real slicer validation incomplete: '+slicing['status'])
        from vascular_processing.print_frame_review import render_review
        summary['png_files']=render_review(out,inputs,chosen,keep,regions,top,cfg,summary)
    except Exception as exc:
        logger.exception('Print-frame stage failed');summary['failures'].append(type(exc).__name__+': '+str(exc))
        from vascular_processing.print_frame_review import render_failure
        render_failure(out,inputs,chosen,cfg,str(exc))
    protection=legacy.verify_snapshot(snapshot);save_json(out/'protected_source_hashes.json',protection)
    summary['protected_file_count']=len(snapshot);summary['hash_mismatches']=protection['changed'];summary['protected_sources_unchanged']=protection['all_unchanged']
    if not protection['all_unchanged']:summary['failures'].append('PROTECTED_INPUT_CHANGED')
    if not summary['failures']:
        summary['status']='READY_FOR_HUMAN_REVIEW' if summary.get('bambu',{}).get('selected') else 'READY_FOR_HUMAN_REVIEW_WITH_BAMBU_PENDING'
    save_json(out/'sacrificial_print_frame_summary.json',summary)
    save_json(out/'validation_results.json',dict(status=summary['status'],failures=summary['failures'],warnings=summary['warnings'],
        protected_sources_unchanged=summary['protected_sources_unchanged'],geometry=summary.get('final_mesh_qc'),human_review_pending=True))
    from vascular_processing.print_frame_review import write_report
    write_report(out,summary)
    if summary['failures']:(out/'failure_summary.md').write_text('# 框架需要调整\n\n'+'\n'.join('- '+v for v in summary['failures'])+'\n\n候选和图像保留在 candidates/ 与 QC/。\n',encoding='utf-8')
    logger.info('Final status: %s',summary['status']);handler.close();logger.removeHandler(handler)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=frame.ROOT/'config/sacrificial_print_frame_BG001.yaml')
    parser.add_argument('--output-root',type=Path)
    parser.add_argument('--skip-slicing',action='store_true',help='Explicit geometry-only diagnostic; default attempts real local Bambu slicing')
    args=parser.parse_args();result=run(args.config,args.output_root,args.skip_slicing)
    print(result['status']);return 0 if result['status'].startswith('READY_FOR_HUMAN_REVIEW') else 2


if __name__=='__main__':raise SystemExit(main())
