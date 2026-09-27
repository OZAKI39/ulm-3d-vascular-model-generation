#!/usr/bin/env python3
"""Default compact prototype: one frozen semantic component, three bounded derivatives."""
from __future__ import annotations

import argparse
from contextlib import ExitStack, contextmanager
import json
from pathlib import Path
import shutil
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from vascular_processing.boundary_review_export import protect_inputs
from vascular_processing.compact_roi import (MODES, load_sources, choose_source, export_candidate, metrics)
from vascular_processing.manufacturing_roi import load_semantic, profile
from vascular_processing.mevo_refinement import branch_edges, frozen_upstream
from vascular_processing.swc_export import read_source
from vascular_processing.topbrain_qc import sha256
from tools.manufacture_bg001_rmca import save, check_protection


@contextmanager
def freeze():
    with frozen_upstream(), ExitStack() as stack:
        for name in ['refine_component', 'prune_diameter']:
            stack.enter_context(patch('vascular_processing.mevo_refinement.' + name,
                                     side_effect=AssertionError('FROZEN_REFINEMENT_MUST_NOT_RUN')))
        yield


def prepare(source, output, config_path, reference):
    if output.exists():
        raise FileExistsError('Choose a new output directory, or resume with --stage model/finish')
    cache, raw, semantic, scores, snapshot = load_semantic(source)
    config = profile(config_path)
    if config['compact_phantom']['default_mode'] != 'BALANCED':
        raise ValueError('Default compact mode must be BALANCED')
    # The complete earlier manufacturing run is a frozen reference, including its manifests.
    for path in reference.rglob('*'):
        if path.is_file():
            snapshot[str(path.resolve())] = sha256(path)
    output.mkdir(parents=True); save(output / 'protected_sources.json', snapshot)
    with protect_inputs(snapshot), freeze():
        source_item, candidates, ranking = choose_source(load_sources(source, raw), semantic, raw, scores, cache, config)
        component = source_item['component']; edges = branch_edges(cache)
        exports = {mode: export_candidate(output / 'candidates' / mode, candidate, source_item, raw,
            raw.graph.subgraph(semantic[component]), edges) for mode, candidate in candidates.items()}
        for mode in ('MINI', 'RICH'):
            if exports[mode]['exports']['compensated']['sha256'] == exports['BALANCED']['exports']['compensated']['sha256']:
                exports[mode]['same_geometry_as'] = 'BALANCED'
        old = json.loads((reference / 'manufacturing_manifest.json').read_text())
        old_graph = read_source(Path(old['exports']['original']['path'])).graph
        old_stats = metrics(old_graph); balanced = exports['BALANCED']['stats']
        manifest = dict(status='COMPACT_CANDIDATES_PREPARED', default_print_candidate='BALANCED',
            candidate_label='COMPACT_PRINT_CANDIDATE', manufacturing_validated=False, print_jobs_submitted=0,
            semantic_source=str(source), source_component=component, source_derived_component=source_item['derived_component'],
            source_level=source_item['level'], source_swc=source_item['path'], source_sha256=source_item['source_sha256'],
            used_part04=component == 4, fallback_to_part03=component == 3,
            selection_reason='First structurally eligible configured component; frozen support minus radius-floor length penalty',
            selection_trials=ranking, profile=config, profile_path=str(config_path), profile_sha256=sha256(config_path),
            candidates=exports, old_reference=dict(role='FULL_CONTEXT_REFERENCE', path=str(reference),
                manifest_sha256=sha256(reference / 'manufacturing_manifest.json'), stats=old_stats),
            comparison=dict(centerline_reduction_percent=100 * (1 - balanced['centerline_length_mm'] / old_stats['centerline_length_mm']),
                outlet_reduction=old_stats['outlet_count'] - balanced['outlet_count'],
                bifurcation_reduction=old_stats['bifurcation_count'] - balanced['bifurcation_count'],
                branch_reduction=old_stats['branch_count'] - balanced['branch_count']),
            source_hashes=snapshot, scale=1., registration_rerun=False, nn_rerun=False,
            semantic_refinement_rerun=False, branch_aggregation_rerun=False,
            constraint_precedence='User-authorized adaptive paths replace fixed 60/80 mm clipping. Single component, branch/outlet/depth, 200 mm total length and 110/140 mm spatial caps remain hard.',
            warnings=['COMPACT_SIZE_BELOW_TARGET'] if exports['BALANCED']['gate']['undersized'] else [])
        save(output / 'compact_manifest.json', manifest)
    check_protection(output)
    print('PREPARED', component, source_item['level'], {k: dict(branches=v['stats']['branch_count'],
        outlets=v['stats']['outlet_count'], length=v['stats']['centerline_length_mm'],
        bbox=v['stats']['bbox']['extents_mm'], gate=v['gate']) for k,v in exports.items()}, flush=True)


def model(output):
    import pyvista as pv
    from tools.refine_bg001_rmca import model_one
    from vascular_processing.print_orientation import cap_native_ports
    m = json.loads((output / 'compact_manifest.json').read_text()); check_protection(output)
    with protect_inputs(m['source_hashes']), freeze():
        for mode in MODES:
            item = m['candidates'][mode]; directory = output / 'candidates' / mode
            export = item['exports']['compensated']
            if sha256(export['path']) != export['sha256']:
                raise ValueError('Compensated SWC changed after preparation')
            if 'model' not in item:
                item['model'] = model_one(Path(export['path']), directory / 'VascularMD')
                save(output / 'compact_manifest.json', m)
            result = item['model']
            if result['status'] != 'success' or not result.get('topology_preserved'):
                item['status'] = 'VASCULARMD_MODEL_FAILED'; continue
            native = Path(result['outputs']['surface_vtk'])
            target = directory / 'surface.vtk'
            if not target.exists():
                shutil.copy2(native, target)
            if sha256(target) != sha256(native):
                raise ValueError('Copied native surface differs')
            closed, qc = cap_native_ports(pv.read(native))
            closed.save(directory / 'surface.stl'); closed.save(directory / 'print_surface.vtp')
            radius = read_source(Path(result['outputs']['swc'])).rows[:, 5]
            minimum = float(2 * radius.min())
            item['post_model_radius_qc'] = dict(min_diameter_mm=minimum, median_diameter_mm=float(2*np.median(radius)),
                max_diameter_mm=float(2*radius.max()), hard_minimum_pass=minimum >= m['profile']['effective']['hard_min_diameter_mm'],
                preferred_minimum_pass=minimum >= m['profile']['effective']['preferred_min_diameter_mm'],
                note='Native fitted smooth SWC; no post-fit clamp')
            item['surface_qc'] = qc; item['surface_vtk'] = str(target); item['surface_stl'] = str(directory / 'surface.stl')
            item['status'] = 'COMPACT_SURFACE_READY'
            save(output / 'compact_manifest.json', m)
    check_protection(output)


def finish(output):
    import pyvista as pv
    from vascular_processing.bambu_manufacturing import discover, slice_candidates, json_object
    from vascular_processing.compact_orientation import optimize_compact
    from vascular_processing.compact_qc import figures, write_report
    m = json.loads((output / 'compact_manifest.json').read_text()); check_protection(output)
    config = m['profile']; balanced = m['candidates']['BALANCED']
    with protect_inputs(m['source_hashes']), freeze():
        discovery_path = output / 'printer_discovery.json'
        if discovery_path.exists():
            discovery = json.loads(discovery_path.read_text())
        else:
            discovery = discover(config); save(discovery_path, discovery)
        m['printer'] = {k:v for k,v in discovery.items() if k != 'flattened'}
        config['printer'].update(model=discovery['model'], build_volume_mm=discovery['build_volume_mm'])
        config['effective']['usable_volume_mm'] = (np.array(discovery['build_volume_mm']) - np.array([
            2*config['printer']['bed_margin_xy_mm'], 2*config['printer']['bed_margin_xy_mm'], config['printer']['top_margin_mm']])).tolist()
        if balanced['status'] != 'COMPACT_SURFACE_READY':
            raise ValueError('BALANCED_NATIVE_SURFACE_UNAVAILABLE: ' + balanced['status'])
        surface = pv.read(output / 'candidates/BALANCED/print_surface.vtp')
        graph = read_source(Path(balanced['exports']['compensated']['path'])).graph
        opath = output / 'orientation/orientation_results.json'
        if opath.exists():
            orientation = json.loads(opath.read_text())
        else:
            orientation = optimize_compact(surface, graph, config, opath.parent); save(opath, orientation)
        m['orientation'] = orientation
        from vascular_processing.swc_export import validate_tree
        orientation['inlet_native_mm'] = graph.nodes[validate_tree(graph)]['coords'][:3].tolist()
        if not orientation['top_candidates']:
            raise ValueError('COMPACT_ROI_EXCEEDS_BUILD_VOLUME')
        selected = orientation['top_candidates'][0]; m['selected_orientation'] = selected
        path = output / 'BG001_RMCA_BALANCED_print_candidate.stl'; shutil.copy2(selected['stl'], path)
        m['selected_stl'] = str(path); save(output / 'print_transform.json', selected)
        spath = output / 'BambuStudio/slice_results.json'
        if spath.exists():
            slicer = json.loads(spath.read_text())
        else:
            # Exactly one orientation of one candidate. MINI/RICH are never sent to CLI.
            slicer = slice_candidates(discovery, config, [selected], spath.parent); save(spath, slicer)
        m['slicer'] = slicer; m['actual_slice_attempts'] = len(slicer['results'])
        if m['actual_slice_attempts'] > 1:
            raise ValueError('ONLY_BALANCED_TOP1_MAY_BE_SLICED')
        m['status'] = 'COMPACT_PRINT_CANDIDATE'
        if slicer['status'] == 'BAMBU_SLICED':
            sliced = slicer['results'][0]
            project = output / 'BG001_RMCA_BALANCED_print_candidate.3mf'
            shutil.copy2(sliced['archive'], project); m['selected_3mf'] = str(project)
            m['status'] = ('PRINT_READY_CANDIDATE' if balanced['gate']['passed'] and balanced['post_model_radius_qc']['hard_minimum_pass']
                else 'COMPACT_SLICED_REQUIRES_ENGINEERING_REVIEW')
        if not balanced['post_model_radius_qc']['hard_minimum_pass']:
            m['warnings'].append('VMD_RADIUS_FLOOR_VIOLATION')
        for mode, candidate in m['candidates'].items():
            if candidate['radius_compensation']['warning']:
                m['warnings'].append(mode + ': LARGE_MANUFACTURING_COMPENSATION')
            if not candidate['post_model_radius_qc']['preferred_minimum_pass']:
                m['warnings'].append(mode + ': VMD_PREFERRED_FLOOR_UNDERSHOOT; hard minimum checked separately')
            if candidate['stats']['centerline_length_mm'] > config['compact_roi']['preferred_total_centerline_mm']['max']:
                m['warnings'].append(mode + ': ABOVE_PREFERRED_TOTAL_LENGTH; below hard 200 mm cap')
            if mode != 'BALANCED' and candidate['exports']['compensated']['sha256'] == balanced['exports']['compensated']['sha256']:
                candidate['same_geometry_as'] = 'BALANCED'
        stem_limits = config['orientation']['preferred_stem_angle_deg']
        if not stem_limits['min'] <= selected['proximal_stem_angle_from_plate_deg'] <= stem_limits['max']:
            m['warnings'].append('STEM_ANGLE_OUTSIDE_PREFERRED_RANGE: soft preference traded against inlet height and support')
        m['warnings'] = list(dict.fromkeys(m['warnings']))
        conf = discovery.get('active_preset_source')
        if conf:
            now = json_object(conf).get('presets', {})
            m['gui_presets_unchanged'] = all(now.get(k) == v for k,v in discovery['active_gui_presets'].items())
            if not m['gui_presets_unchanged']:
                raise ValueError('GUI presets changed')
        figures(output, m)
        m['protection'] = check_protection(output); save(output / 'compact_manifest.json', m)
        write_report(output, m)
    finalize_artifacts(output)
    print('FINISHED', m['status'], m['actual_slice_attempts'], flush=True)


def finalize_artifacts(output):
    m = json.loads((output / 'compact_manifest.json').read_text())
    from vascular_processing.compact_qc import write_report
    m['protection'] = check_protection(output)
    if (output / 'QC/test_results.json').exists():
        m['tests'] = json.loads((output / 'QC/test_results.json').read_text())
    save(output / 'compact_manifest.json', m); write_report(output, m)
    files = {str(p.relative_to(output)): sha256(p) for p in output.rglob('*') if p.is_file() and p.name != 'artifact_manifest.json'}
    code = [ROOT/'tools/compact_bg001_rmca.py', ROOT/'config/compact_manufacturing_profile.yaml',
            *ROOT.glob('vascular_processing/compact_*.py'), ROOT/'tests/test_compact_manufacturing.py']
    save(output / 'artifact_manifest.json', dict(files=files, code={str(p.relative_to(ROOT)):sha256(p) for p in code if p.is_file()}))


def main():
    base = ROOT / 'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA'
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=base/'refined_roi')
    p.add_argument('--output', type=Path, default=base/'compact_manufacturing_roi')
    p.add_argument('--reference', type=Path, default=base/'manufacturing_roi')
    p.add_argument('--config', type=Path, default=ROOT/'config/compact_manufacturing_profile.yaml')
    p.add_argument('--stage', choices=['prepare','model','finish','all','verify'], default='all')
    args = p.parse_args(); output = args.output.resolve()
    if args.stage in ['prepare','all']:
        prepare(args.input.resolve(), output, args.config.resolve(), args.reference.resolve())
    if args.stage in ['model','all']: model(output)
    if args.stage in ['finish','all']: finish(output)
    if args.stage == 'verify':
        check_protection(output)
        artifact = json.loads((output/'artifact_manifest.json').read_text())
        assert all(sha256(output/rel) == h for rel,h in artifact['files'].items())
        assert all(sha256(ROOT/rel) == h for rel,h in artifact['code'].items())
        print('All compact artifacts and protected sources verified')


if __name__ == '__main__':
    main()
