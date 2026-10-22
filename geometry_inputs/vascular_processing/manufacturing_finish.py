"""Finish the manufacturing derivative: closed mesh, orientations, real slicing, QC."""
from pathlib import Path
import json
import shutil

import numpy as np
import pyvista as pv

from .bambu_manufacturing import discover, slice_candidates
from .boundary_review_export import protect_inputs
from .manufacturing_roi import load_semantic
from .manufacturing_qc import write_figures, write_report
from .mevo_refinement import frozen_upstream
from .print_orientation import cap_native_ports, optimize, calibration_coupon
from .swc_export import read_source
from .topbrain_qc import sha256


def choose_candidate(orientation,slicer,config):
    candidates=orientation['top_candidates']
    by_rank={r['rank']:r for r in candidates}
    successes=[r for r in slicer['results'] if r['status']=='BAMBU_SLICED']
    usable=[r for r in successes if r.get('estimated_print_time_seconds') is not None]
    if not usable:
        return candidates[0],dict(method='SURFACE_HEURISTIC',reason='No successful slice with exposed numeric time')
    times=np.array([r['estimated_print_time_seconds'] for r in usable])
    # Total material is an explicitly labelled proxy, never reported as measured
    # support mass. Model volume is unchanged across rigid orientations.
    have_material=all(r.get('filament_used_g') is not None for r in usable)
    burdens=np.array([r['filament_used_g'] if have_material else by_rank[r['rank']]['downward_support_area_mm2'] for r in usable])
    normalize=lambda a:(a-a.min())/max(float(np.ptp(a)),1e-12)
    weight=config['slicer']
    score=weight['support_selection_weight']*normalize(burdens)+weight['time_selection_weight']*normalize(times)
    index=min(range(len(usable)),key=lambda i:(score[i],by_rank[usable[i]['rank']]['score'],usable[i]['rank']))
    selected=usable[index]
    return by_rank[selected['rank']],dict(method='REAL_SLICE_TIME_AND_MATERIAL_PROXY' if have_material else 'REAL_SLICE_TIME_AND_HEURISTIC_AREA',
        support_burden_field='total filament g (proxy, not isolated support mass)' if have_material else 'surface heuristic mm2',
        weighted_scores=[dict(rank=r['rank'],score=float(score[i])) for i,r in enumerate(usable)],
        selected_slice=selected)


def finish(output):
    from tools.manufacture_bg001_rmca import save,check_protection
    output=Path(output);m=json.loads((output/'manufacturing_manifest.json').read_text())
    snapshot=json.loads((output/'protected_sources.json').read_text());check_protection(output)
    config=m['profile']
    for export in m['exports'].values():
        if sha256(export['path'])!=export['sha256']:
            raise ValueError('Manufacturing input was changed: '+export['path'])
    with protect_inputs(snapshot),frozen_upstream():
        discovery_path=output/'printer_discovery.json'
        if discovery_path.is_file():discovery=json.loads(discovery_path.read_text())
        else:
            discovery=discover(config);save(discovery_path,discovery)
        m['printer']={k:v for k,v in discovery.items() if k!='flattened'}
        config['printer']['model']=discovery['model']
        config['printer']['build_volume_mm']=discovery['build_volume_mm']
        config['effective']['usable_volume_mm']=(np.array(discovery['build_volume_mm'])-np.array([
            2*config['printer']['bed_margin_xy_mm'],2*config['printer']['bed_margin_xy_mm'],config['printer']['top_margin_mm']])).tolist()
        m['profile']=config;save(output/'effective_profile.json',config)
        comp=m['models']['compensated']
        if comp['status']!='success' or not comp.get('topology_preserved'):
            m['status']='VASCULARMD_MANUFACTURING_SURFACE_UNAVAILABLE'
            save(output/'manufacturing_manifest.json',m)
            check_protection(output)
            return m
        vtk=Path(comp['outputs']['surface_vtk']);surface_path=output/'BG001_RMCA_print_surface_closed.vtp'
        if m.get('native_surface_sha256') and sha256(vtk)!=m['native_surface_sha256']:
            raise ValueError('Cached native surface changed; use a new output run')
        if surface_path.is_file():
            surface=pv.read(surface_path);surface_qc=json.loads((output/'surface_qc.json').read_text())
            if m.get('orientation_input') and sha256(surface_path)!=m['orientation_input']['sha256']:
                raise ValueError('Cached print surface changed; use a new output run')
        else:
            surface,surface_qc=cap_native_ports(pv.read(vtk))
            surface.save(surface_path);surface.save(output/'BG001_RMCA_print_surface_closed.stl')
            save(output/'surface_qc.json',surface_qc)
        m['surface_qc']=surface_qc;m['native_surface_sha256']=sha256(vtk)
        m['orientation_input']=dict(path=str(surface_path),sha256=sha256(surface_path),
            source='Compensated native VascularMD sidewall, original ports capped only on print derivative')
        graph=read_source(Path(m['exports']['compensated']['path'])).graph
        orientation_path=output/'orientation/orientation_results.json'
        if orientation_path.is_file():orientation=json.loads(orientation_path.read_text())
        else:
            print('ORIENTATION SEARCH START',flush=True)
            orientation=optimize(surface,graph,config,orientation_path.parent)
            save(orientation_path,orientation)
        m['orientation']=orientation
        if not orientation['top_candidates']:
            m['status']='ROI_EXCEEDS_BUILD_VOLUME';save(output/'manufacturing_manifest.json',m);check_protection(output)
            return m
        if 'PRINT_ORIENTATION_OPTIMIZED' not in m['completed_states']:m['completed_states'].append('PRINT_ORIENTATION_OPTIMIZED')
        m['status']='PRINT_ORIENTATION_OPTIMIZED';save(output/'manufacturing_manifest.json',m)
        print('ORIENTATION SEARCH END',orientation['tested'],'candidates',flush=True)
        slicer_path=output/'BambuStudio/slice_results.json'
        if slicer_path.is_file():slicer=json.loads(slicer_path.read_text())
        else:
            slicer=slice_candidates(discovery,config,orientation['top_candidates'],slicer_path.parent)
            save(slicer_path,slicer)
        m['slicer']=slicer
        selected,selection=choose_candidate(orientation,slicer,config)
        m['selected_orientation']=selected;m['orientation_selection']=selection
        baseline=orientation['baseline']
        m['orientation_comparison']=dict(
            height_reduction_mm=baseline['height_mm']-selected['height_mm'],
            height_reduction_percent=100*(1-selected['height_mm']/baseline['height_mm']),
            unsupported_area_reduction_mm2=baseline['downward_support_area_mm2']-selected['downward_support_area_mm2'],
            unsupported_area_reduction_percent=100*(1-selected['downward_support_area_mm2']/baseline['downward_support_area_mm2']),
            support_area_is_heuristic=True)
        destination=output/'BG001_RMCA_print_candidate.stl';shutil.copy2(selected['stl'],destination)
        m['selected_stl']=str(destination)
        save(output/'print_transform.json',dict(**selected,
            convention='column vector p_print = R @ p_native + translation; units mm',
            source_surface_sha256=sha256(vtk),global_scale=1.0))
        if slicer['status']=='BAMBU_SLICED':
            if 'BAMBU_SLICED' not in m['completed_states']:m['completed_states'].append('BAMBU_SLICED')
            matching=[r for r in slicer['results'] if r['rank']==selected['rank'] and r['status']=='BAMBU_SLICED']
            if matching:
                if not matching[0].get('geometry_qc',{}).get('rotation_scale_position_preserved'):
                    raise ValueError('Selected sliced archive has no successful geometry verification')
                project=output/'BG001_RMCA_print_candidate.3mf';shutil.copy2(matching[0]['archive'],project)
                m['selected_3mf']=str(project)
                if m['roi_size_gate']['passed'] and m['post_model_radius_qc']['hard_minimum_pass']:
                    m['status']='PRINT_READY_CANDIDATE'
                    if m['status'] not in m['completed_states']:m['completed_states'].append(m['status'])
        if not m['post_model_radius_qc']['hard_minimum_pass']:m['status']='VMD_RADIUS_FLOOR_VIOLATION'
        m['manufacturing_validated']=False;m['print_jobs_submitted']=0
        m['warnings']=[]
        if not m['post_model_radius_qc']['preferred_floor_pass']:
            m['warnings'].append('VMD_PREFERRED_FLOOR_UNDERSHOOT: native smoothing lowered the preferred floor; hard minimum checked independently')
        if not m['roi_size_gate']['preferred_size']:
            m['warnings'].append('AABB longest dimension below preferred 100 mm; minimum size gate passed; no scaling')
        coupon=output/'bambu_abs_vessel_calibration_coupon.stl'
        if coupon.exists():m['calibration_coupon']=json.loads((output/'calibration_coupon.json').read_text())
        else:
            m['calibration_coupon']=calibration_coupon(config,coupon);save(output/'calibration_coupon.json',m['calibration_coupon'])
        cache,raw,semantic,_,_=load_semantic(Path(m['semantic_source']))
        core=raw.graph.subgraph(m['core_original_node_ids']).copy()
        from .boundary_review_export import read_csv
        raw_path=Path(m['exports']['original']['path'])
        mapping={int(r['swc_id']):int(r['original_swc_id']) for r in read_csv(raw_path.with_name(raw_path.stem+'_mapping.csv'))}
        original=raw.graph.subgraph(mapping.values()).copy()
        compgraph=original.copy()
        for n,old in mapping.items():compgraph.nodes[old]['coords']=graph.nodes[n]['coords'].copy()
        write_figures(output/'QC',m,core,original,compgraph,surface)
        m['protection']=check_protection(output)
        save(output/'manufacturing_manifest.json',m)
        write_report(output,m)
        (output/'README.md').write_text(
            '# BG001 RMCA manufacturing candidate\n\n'
            'NOT manufacturer guaranteed limits. Must be calibrated experimentally.\n\n'
            f"Status: **{m['status']}**. Physical printing/PDMS/ABS dissolution have NOT been validated. No print jobs submitted.\n\n"
            '- Layer A: immutable `../refined_roi/semantic_refined/`. Historical `diameter_075/` is retained and not used here.\n'
            '- Layer B: `BG001_RMCA_print_roi_original_radius.swc`: native manufacturing topology and raw radii.\n'
            '- Layer C: `BG001_RMCA_print_roi_compensated.swc`: MANUFACTURING_COMPENSATED_GEOMETRY.\n'
            '- `VascularMD/`: both native fits, QC and original native surfaces.\n'
            '- `BG001_RMCA_print_candidate.stl`: selected compensated print derivative, ports capped, rotated/translated, scale=1.\n'
            '- `BG001_RMCA_print_candidate.3mf`: present only when real Bambu slicing succeeded.\n'
            '- `BambuStudio/`: exact commands, profile inheritance, actual slice metadata and all candidate 3MFs.\n'
            '- `orientation/`: preserved identity baseline, Top 5 STL, all transforms and scores.\n'
            '- `bambu_abs_vessel_calibration_coupon.stl`: 25 separate calibration cylinders; CSV maps diameter/angle to row/column.\n\n'
            'Run from project root (choose a new output directory):\n\n```bash\n'
            '/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python tools/manufacture_bg001_rmca.py --output outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/manufacturing_run02\n```\n\n'
            'Alternative detail profile: `--config config/bambu_abs_0p2_detail.yaml`. This does not change your physical nozzle. '
            'Stages `--stage prepare`, `--stage model`, `--stage finish` allow checkpointed operation; `all` is the default. '
            'A new full run refuses an existing output directory. Existing semantic/refinement results are always protected.\n',encoding='utf-8')
    check_protection(output)
    root=Path(__file__).resolve().parents[1]
    code_files=[*root.glob('vascular_processing/manufacturing*.py'),root/'vascular_processing/print_orientation.py',
                root/'vascular_processing/bambu_manufacturing.py',root/'tools/manufacture_bg001_rmca.py',
                root/'config/bambu_manufacturing_profile.yaml',root/'config/bambu_abs_0p2_detail.yaml']
    save(output/'artifact_manifest.json',dict(files={str(p.relative_to(output)):sha256(p) for p in sorted(output.rglob('*'))
        if p.is_file() and p.name!='artifact_manifest.json'},code={str(p.relative_to(root)):sha256(p) for p in code_files}))
    print('MANUFACTURING COMPLETE',m['status'],output,flush=True)
    return m
