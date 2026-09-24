"""TopBrain anatomical MeVO entry; human/BraVa entry and renderer remain read-only."""
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
    return parser.parse_args(argv)


def main(argv=None):
    args = arguments(argv)
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


if __name__ == '__main__':
    raise SystemExit(main())
