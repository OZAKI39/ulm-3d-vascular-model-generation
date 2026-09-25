#!/usr/bin/env python3
"""Frozen semantic core -> native connected topology -> compensated print geometry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from vascular_processing.boundary_review_export import protect_inputs
from vascular_processing.mevo_refinement import branch_edges, frozen_upstream
from vascular_processing.manufacturing_roi import (
    compensate_radius, coordinates, export_manufacturing, extend_proximal, geometry_stats,
    load_semantic, minimal_subtree, profile, prune_twigs, roi_size_gate, select_core)
from vascular_processing.swc_export import read_source
from vascular_processing.topbrain_qc import sha256


def save(path, value):
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def check_protection(output):
    snapshot = json.loads((output/'protected_sources.json').read_text())
    changed = [p for p,h in snapshot.items() if not Path(p).is_file() or sha256(p)!=h]
    result = dict(protected_files=len(snapshot),all_unchanged=not changed,changed=changed,
        registration_rerun=False,nn_rerun=False,branch_aggregation_rerun=False)
    save(output/'protection_verification.json',result)
    if changed:
        raise ValueError('Protected manufacturing source changed: '+str(changed))
    return result


def prepare(source_dir, output, config_path):
    if output.exists():
        raise FileExistsError('Output exists: choose a new directory or resume --stage model/finish')
    cache,raw,semantic,scores,snapshot = load_semantic(source_dir)
    config = profile(config_path)
    output.mkdir(parents=True)
    save(output/'protected_sources.json',snapshot)
    with protect_inputs(snapshot),frozen_upstream():
        selected,ranking = select_core(semantic,scores,raw,config)
        core_nodes = set().union(*(semantic[i] for i in selected))
        optional = sorted(set(config['roi'].get('optional_context_components',[]))-set(selected))
        if any(i not in semantic for i in optional):
            raise ValueError('Unknown optional semantic component ID')
        selected_nodes = core_nodes | set().union(*(semantic[i] for i in optional))
        core_graph = raw.graph.subgraph(core_nodes).copy()
        rmca = raw.graph.subgraph(cache.node_graph).copy()
        connected = minimal_subtree(rmca,selected_nodes)
        extended,extension = extend_proximal(connected,rmca,config['roi'])
        topology,removed = prune_twigs(extended,config)
        compensated,blend = compensate_radius(topology,config['effective']['radius_floor_mm'],
                                               config['manufacturing']['radius_blend_length_mm'])
        edge_branch = branch_edges(cache)
        context = set(topology)-core_nodes
        context_branches = sorted({edge_branch[e] for e in topology.edges if not all(n in core_nodes for n in e)})
        paths = dict(original=output/'BG001_RMCA_print_roi_original_radius.swc',
                     compensated=output/'BG001_RMCA_print_roi_compensated.swc')
        exports = {name:export_manufacturing(path,topology if name=='original' else compensated,
            raw,core_nodes,edge_branch,name=='compensated') for name,path in paths.items()}
        delta = np.array([compensated.nodes[n]['coords'][3]-topology.nodes[n]['coords'][3] for n in topology])
        manifest = dict(status='MANUFACTURING_RADIUS_COMPENSATED',
            completed_states=['ANATOMICAL_TRANSFER_CANDIDATE','MANUFACTURING_ROI_CONNECTED','MANUFACTURING_RADIUS_COMPENSATED'],
            geometry_label='MANUFACTURING_COMPENSATED_GEOMETRY',manufacturing_validated=False,
            semantic_source=str(source_dir),source_hashes=cache.provenance,
            semantic_source_hashes={p:h for p,h in snapshot.items() if str(source_dir) in p},
            profile_path=str(config_path),profile=config,selected_components=selected,component_selection=ranking,
            included_optional_components=optional,optional_component_label='OPTIONAL_SEMANTIC_COMPONENT',
            core_original_node_ids=sorted(core_nodes),context_original_node_ids=sorted(context),
            context_branches=context_branches,proximal_context=extension,
            semantic_core=geometry_stats(core_graph),before_twig_removal=geometry_stats(extended),
            manufacturing_topology=geometry_stats(topology),manufacturing_geometry=geometry_stats(compensated),
            roi_size_gate=roi_size_gate(geometry_stats(topology),config),
            radius_floor_mm=config['effective']['radius_floor_mm'],removed_terminal_twigs=removed,
            radius_inflation=dict(modified_samples=int(np.count_nonzero(delta>1e-12)),
                sum_added_radius_mm=float(delta.sum()),mean_added_radius_mm=float(delta.mean()),
                maximum_added_radius_mm=float(delta.max()),
                note='Sum is a sample-based statistic, not added volume or centerline length'),
            blend_windows=blend,exports=exports,scale=1.0,
            retired_default='diameter_075 retained unchanged as historical experiment',
            orientation_baseline=dict(method='IDENTITY_EXPORTED_GEOMETRY',rotation=np.eye(3).tolist(),
                audit='Existing BraVa viewer uses camera.Azimuth only; doc_visualize_v3 uses PCA camera framing, not exported mesh rotation. Native VMD mesh is saved in source coordinates.'),
            models={},files={})
        save(output/'manufacturing_manifest.json',manifest)
    check_protection(output)
    print('PREPARED',json.dumps(dict(selected=selected,context=context_branches,
         size=manifest['manufacturing_topology'],inflation=manifest['radius_inflation'],twigs=len(removed))),flush=True)


def model(output):
    from tools.refine_bg001_rmca import model_one
    manifest=json.loads((output/'manufacturing_manifest.json').read_text())
    snapshot=json.loads((output/'protected_sources.json').read_text())
    check_protection(output)
    with protect_inputs(snapshot),frozen_upstream():
        for name,export in manifest['exports'].items():
            if name in manifest['models']:
                continue
            assert sha256(export['path'])==export['sha256']
            manifest['models'][name]=model_one(Path(export['path']),output/'VascularMD'/name)
            save(output/'manufacturing_manifest.json',manifest)
        comp=manifest['models']['compensated']
        if comp['outputs'].get('swc'):
            rows=read_source(Path(comp['outputs']['swc'])).rows
            minimum=float(rows[:,5].min());cfg=manifest['profile']
            hard=cfg['effective']['hard_min_diameter_mm'];floor=cfg['effective']['radius_floor_mm']
            tolerance=cfg['manufacturing']['post_model_radius_tolerance_mm']
            manifest['post_model_radius_qc']=dict(min_radius_mm=minimum,min_diameter_mm=2*minimum,
                median_diameter_mm=float(2*np.median(rows[:,5])),max_diameter_mm=float(2*rows[:,5].max()),
                floor_undershoot_mm=max(0.,floor-minimum),preferred_floor_pass=minimum>=floor-tolerance,
                hard_minimum_pass=2*minimum>=hard,
                status='VMD_RADIUS_FLOOR_VIOLATION' if 2*minimum<hard else 'PASS_HARD_MINIMUM',
                source='Native smooth SWC; surface ports/closedness checked separately, no post-fit radius clamp')
        else:
            manifest['post_model_radius_qc']=dict(status='NATIVE_MODEL_UNAVAILABLE',hard_minimum_pass=False)
        save(output/'manufacturing_manifest.json',manifest)
    check_protection(output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/refined_roi')
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/manufacturing_roi')
    parser.add_argument('--config',type=Path,default=ROOT/'config/bambu_manufacturing_profile.yaml')
    parser.add_argument('--stage',choices=['prepare','model','finish','all'],default='all')
    args=parser.parse_args();output=args.output.resolve()
    if args.stage in ['prepare','all']:
        prepare(args.input.resolve(),output,args.config.resolve())
    if args.stage in ['model','all']:
        model(output)
    if args.stage in ['finish','all']:
        from vascular_processing.manufacturing_finish import finish
        finish(output)


if __name__=='__main__':
    main()
