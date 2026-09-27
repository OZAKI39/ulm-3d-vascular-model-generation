"""Compact BraVa ROI viewer; full global vessels and the human renderer remain intact."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import yaml

from vascular_processing.topbrain_dataset import discover, load_case
from vascular_processing.topbrain_pipeline import TopBrainOptions, process_case
from vascular_processing.topbrain_qc import TopBrainError, write_json

PROJECT_ROOT = Path(__file__).resolve().parent


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',choices=['compact-brava','manufacturing-stl','topbrain','topbrain-brava'],default=None,
                        help='Default: current compact BraVa candidates (BALANCED selected)')
    parser.add_argument('--stl-file',type=Path,default=PROJECT_ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/manufacturing_roi/BG001_RMCA_print_candidate.stl',
                        help='Saved print candidate for --source manufacturing-stl')
    parser.add_argument('--compact-results',type=Path,default=PROJECT_ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi')
    parser.add_argument('--brava-results',type=Path,default=PROJECT_ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001')
    parser.add_argument('--brava-side',choices=['LMCA','RMCA'])
    parser.add_argument('--brava-view',choices=['strict','surface'],default='strict')
    parser.add_argument('--config', type=Path, default=PROJECT_ROOT/'configs/topbrain_mevo.yaml')
    parser.add_argument('--topbrain-root', type=Path)
    parser.add_argument('--case-id')
    parser.add_argument('--modality', choices=['mr','ct'])
    parser.add_argument('--dataset-version', choices=['auto','2025','ta36'])
    parser.add_argument('--labelmap', type=Path)
    parser.add_argument('--annotation-source', choices=['gold_annotation','model_prediction'],default='gold_annotation')
    parser.add_argument('--roi')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--recompute', action='store_true')
    parser.add_argument('--discover', action='store_true', help='List actual paired cases without geometry processing')
    parser.add_argument('--show-native-labels', action='store_true')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--no-gui',action='store_true',help='Process and export only')
    group.add_argument('--off-screen',action='store_true',help='Render an off-screen dual-viewport preview')
    parser.add_argument('--smoke-gui-seconds',type=float,default=0.,help='Actual GUI smoke test with ROI switch and timed close')
    args=parser.parse_args(argv)
    # Preserve explicit legacy TopBrain requests, while launching current compact
    # candidates when the script is run without dataset/source arguments.
    if args.source is None:
        args.source='topbrain' if any((args.topbrain_root,args.case_id,args.modality,args.dataset_version,args.labelmap,args.discover,args.show_native_labels)) else 'compact-brava'
    return args


def main(argv=None):
    args = arguments(argv)
    if args.source=='compact-brava':
        return show_compact(args)
    if args.source=='manufacturing-stl':
        return show_manufacturing_stl(args)
    if args.source=='topbrain-brava':
        return show_brava(args)
    config = yaml.safe_load(args.config.read_text()) if args.config.is_file() else {}
    root = args.topbrain_root or Path(os.environ.get('TOPBRAIN_ROOT',config.get('topbrain_root','data/TopBrain')))
    if not root.is_absolute():
        root = PROJECT_ROOT/root
    output = args.output_dir or PROJECT_ROOT/config.get('output_dir','outputs/topbrain_mevo')
    modality = args.modality or config.get('modality','mr')
    version = args.dataset_version or config.get('dataset_version','2025')
    try:
        pairs = discover(root,version)
        if args.discover:
            print(json.dumps([{'case_id':p.case_id,'modality':p.modality,'version':p.version,
                               'image':str(p.image),'label':str(p.label)} for p in pairs],indent=2))
            return 0
        selected = [p for p in pairs if p.modality==modality]
        if not selected:
            raise TopBrainError('CASE_NOT_FOUND',f'No {modality} cases found')
        case_id = args.case_id or config.get('case_id') or selected[0].case_id
        current = next((i for i,p in enumerate(selected) if p.case_id==str(case_id).zfill(3)),None)
        if current is None:
            raise TopBrainError('CASE_NOT_FOUND',f'No case {case_id} for {modality}')
        options = TopBrainOptions()
        options.validate()
        while True:
            case = load_case(root,selected[current].case_id,modality,version=version,labelmap=args.labelmap,
                             annotation_source=args.annotation_source)
            print(f'[TopBrain]\nDataset root: {root}\nDataset version: {case.paths.version}\nCase: {case.paths.case_id}\nModality: {modality}\nImage: {case.paths.image}\nLabel: {case.paths.label}\nShape: {case.labels.shape}\nSpacing (mm): {case.spacing_mm}\nOrientation: {case.provenance["orientation"]}',flush=True)
            run,rois,manifest = process_case(case,output,options=options,recompute=args.recompute,
                                           cache_root=PROJECT_ROOT/'outputs/topbrain_mevo_cache')
            print(f'Output: {run}\nCache hit: {manifest["cache_hit"]}\nStatus: {manifest["status"]}',flush=True)
            if args.no_gui:
                break
            from vascular_processing.topbrain_display_adapter import KEY_BINDINGS, TopBrainViewer
            print('Key bindings: '+json.dumps(KEY_BINDINGS,ensure_ascii=False),flush=True)
            display = bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY') or sys.platform in {'win32','darwin'})
            show = display and not args.off_screen
            if not display and not args.off_screen:
                print('DISPLAY_NOT_AVAILABLE / SKIPPED_NO_DISPLAY: attempting off-screen rendering',flush=True)
            viewer = TopBrainViewer(case,rois,run,manifest,initial_roi=args.roi,show=show,show_native_labels=args.show_native_labels)
            report = viewer.run_window(show=show,smoke_seconds=args.smoke_gui_seconds)
            if not display:
                report['gui_test_status'] = 'SKIPPED_NO_DISPLAY'
                write_json(run/'topbrain_ui_compatibility.json',report)
            if not report['case_delta']:
                break
            current = (current+report['case_delta'])%len(selected)
        return 0 if manifest['status']=='PASS' else 2
    except TopBrainError as exc:
        print(str(exc),file=sys.stderr,flush=True)
        return 2
    except Exception as exc:
        print(f'TOPBRAIN_FAILED: {type(exc).__name__}: {exc}',file=sys.stderr,flush=True)
        return 1


def show_compact(args):
    """Display saved compact data in native coordinates; never run upstream jobs."""
    try:
        from vascular_processing.compact_display_adapter import CompactBraVaViewer, load_compact
        if args.no_gui:
            manifest,original,records=load_compact(args.compact_results)
            print(json.dumps(dict(source='compact-brava',default_roi=manifest['default_print_candidate'],
                full_global_nodes=len(original),full_global_edges=original.number_of_edges(),
                rois=[dict(name=r.roi_id,nodes=r.node_count,edges=r.edge_count,branches=r.branch_count) for r in records]),indent=2))
            return 0
        display=bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY') or sys.platform in {'win32','darwin'})
        show=display and not args.off_screen
        viewer=CompactBraVaViewer(args.compact_results,initial_roi=args.roi,view=args.brava_view,
                                 show=show,output_dir=args.output_dir)
        report=viewer.run_window(show=show,smoke_seconds=args.smoke_gui_seconds)
        print(json.dumps(report,indent=2,ensure_ascii=False))
        return 2 if args.smoke_gui_seconds and not report['gui_smoke_passed'] else 0
    except Exception as exc:
        print(f'COMPACT_VIEW_FAILED: {type(exc).__name__}: {exc}',file=sys.stderr,flush=True)
        return 2


def show_manufacturing_stl(args):
    """Load the saved STL, preserving the current UI and full global context."""
    try:
        from vascular_processing.manufacturing_stl_display import ManufacturingSTLViewer, load_print_stl
        if args.no_gui:
            _,_,_,_,provenance=load_print_stl(args.stl_file)
            print(json.dumps(provenance,indent=2,ensure_ascii=False));return 0
        display=bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY') or sys.platform in {'win32','darwin'})
        show=display and not args.off_screen
        viewer=ManufacturingSTLViewer(args.stl_file,initial_roi=args.roi,show=show,output_dir=args.output_dir)
        report=viewer.run_window(show=show,smoke_seconds=args.smoke_gui_seconds)
        print(json.dumps(report,indent=2,ensure_ascii=False))
        return 2 if args.smoke_gui_seconds and not report['gui_smoke_passed'] else 0
    except Exception as exc:
        print(f'MANUFACTURING_STL_VIEW_FAILED: {type(exc).__name__}: {exc}',file=sys.stderr,flush=True)
        return 2


def show_brava(args):
    """Read saved results only; no registration or modeling in the viewer."""
    try:
        path=args.brava_results/'summary.json'
        if not path.is_file():raise ValueError('BRAVA_RESULTS_MISSING: run tools/transfer_topbrain_to_brava.py first')
        summary=json.loads(path.read_text())
        if args.no_gui:
            print(json.dumps(summary,indent=2,ensure_ascii=False))
            sides=[summary.get(args.brava_side,{})] if args.brava_side else summary.values()
            return 0 if any(r.get('roi_components') for r in sides) else 2
        from vascular_processing.brava_display_adapter import BraVaViewer
        display=bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY') or sys.platform in {'win32','darwin'})
        show=display and not args.off_screen
        viewer=BraVaViewer(args.brava_results,side=args.brava_side,initial_roi=args.roi,view=args.brava_view,show=show)
        report=viewer.run_window(show=show,smoke_seconds=args.smoke_gui_seconds)
        print(json.dumps(report,indent=2,ensure_ascii=False))
        if args.smoke_gui_seconds and not report['gui_smoke_passed']:return 2
        return 0
    except Exception as exc:
        print(f'BRAVA_VIEW_FAILED: {type(exc).__name__}: {exc}',file=sys.stderr,flush=True)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
