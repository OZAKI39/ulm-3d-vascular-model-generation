"""Freeze explicit human review coverage; exports must never reread live labels."""
from pathlib import Path
import json
import shutil
from annotation_store import CLASSES, load_annotation, save_annotation, sha, export_cellpose
from annotation_regions import reviewed_boxes, reviewed_counts


def revoke_full_frame(directory, record):
    """User-requested scope correction; keep exact pixels and historical revisions."""
    meta, masks = load_annotation(directory, record)
    if not meta or not meta.get('expert_confirmed'):
        raise ValueError('Expected a currently confirmed whole frame: ' + record['id'])
    revised = save_annotation(directory, record, masks, meta['reviewer'], False, meta.get('track_ids'), meta.get('notes', ''),
                              review_scope='REGIONS', workspace_roi=meta.get('workspace_roi_yxyx'),
                              proposal_sources=meta.get('proposal_sources'))
    return meta, revised


def freeze_reviewed(root, destination, execution=None):
    root = Path(root); destination = Path(destination); destination.mkdir(parents=True, exist_ok=False)
    records = json.loads((root / 'data/candidate_manifest.json').read_text())['rows']
    selected = []; exclusions = []; frozen_root = destination / 'snapshot'; hashes = {}
    for record in records:
        meta, masks = load_annotation(root / 'data/annotations', record)
        if not meta or not reviewed_boxes(meta, masks):
            exclusions.append({'id': record['id'], 'reason': 'NO_CONFIRMED_COVERAGE'}); continue
        if record['split'] == 'QUARANTINE_UNASSIGNED' or 'SOFTWARE' in (meta['reviewer'] + ' ' + meta.get('notes', '')).upper():
            exclusions.append({'id': record['id'], 'reason': 'QUARANTINE_OR_SOFTWARE_FIXTURE'}); continue
        image = Path(record['image'])
        if image.is_absolute() or '..' in image.parts or Path(record['id']).name != record['id']:
            raise ValueError('Invalid candidate path')
        if sha(root / image) != record['image_sha256']:
            raise ValueError('Image identity changed: ' + record['id'])
        mask = root / 'data/annotations' / record['id'] / meta['mask_file']
        target_mask = frozen_root / 'data/annotations' / record['id'] / meta['mask_file']
        target_mask.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(mask, target_mask)
        target_image = frozen_root / image; target_image.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / image, target_image)
        if sha(target_mask) != meta['mask_sha256'] or sha(target_image) != record['image_sha256']:
            raise ValueError('Snapshot copy checksum mismatch')
        # Store exactly the revision loaded above, not a potentially newer live pointer.
        for path in (target_mask.with_suffix('.json'), target_mask.parent / 'latest.json'):
            path.write_text(json.dumps(meta, indent=2))
        selected.append({'record': record, 'annotation': meta,
                         'frozen_mask_file': str(target_mask.relative_to(destination)),
                         'frozen_image_file': str(target_image.relative_to(destination))})
    if not selected: raise ValueError('No reviewed frames or regions to freeze')
    (frozen_root / 'data/candidate_manifest.json').write_text(json.dumps({'rows': [r['record'] for r in selected]}, indent=2))
    for path in sorted(frozen_root.rglob('*')):
        if path.is_file(): hashes[str(path.relative_to(destination))] = sha(path)
    result = {'execution': execution, 'protocol': json.loads((root / 'protocol.json').read_text()),
              'rows': selected, 'excluded': exclusions, 'snapshot_sha256': hashes,
              'coverage_semantics': 'Only explicitly reviewed pixels; no new approvals inferred'}
    (destination / 'frozen.json').write_text(json.dumps(result, indent=2))
    return result


def export_frozen(destination):
    destination = Path(destination); freeze = json.loads((destination / 'frozen.json').read_text())
    for name, expected in freeze['snapshot_sha256'].items():
        if sha(destination / name) != expected: raise ValueError('Frozen input changed: ' + name)
    records = [row['record'] for row in freeze['rows']]
    return {split: export_cellpose(destination / 'snapshot', records, destination / ('cellpose_' + split), split)
            for split in ('train', 'validation')}


def dataset_readiness(root):
    root = Path(root); records = json.loads((root / 'data/candidate_manifest.json').read_text())['rows']
    summary = {}; frames = []; next_validation = []
    for index, record in enumerate(records):
        meta, masks = load_annotation(root / 'data/annotations', record)
        boxes = reviewed_boxes(meta, masks) if meta else []
        counts = reviewed_counts(meta, masks) if meta else {key: 0 for key in CLASSES}
        row = {'index_1based': index+1, 'id': record['id'], 'split': record['split'],
               'status': meta['status'] if meta else 'UNREVIEWED', 'reviewed_boxes_yxyx': boxes,
               'reviewed_counts': counts, 'saved_counts': meta['instance_counts'] if meta else {},
               'revision': meta['revision'] if meta else None, 'track_ids': len(meta.get('track_ids', {})) if meta else 0}
        frames.append(row)
        split = summary.setdefault(record['split'], {'saved_frames': 0, 'reviewed_frames': 0,
                                                     'reviewed_regions': 0, 'reviewed_counts': {k: 0 for k in CLASSES}})
        split['saved_frames'] += int(meta is not None); split['reviewed_frames'] += int(bool(boxes))
        split['reviewed_regions'] += len(boxes)
        for key in CLASSES: split['reviewed_counts'][key] += counts[key]
        if record['split'] == 'validation' and not boxes and len(next_validation) < 2:
            next_validation.append({'index_1based': index+1, 'id': record['id']})
    validation = summary.get('validation', {})
    return {'split_summary': summary, 'frames': frames, 'next_validation_candidates': next_validation,
            'independent_evaluation_ready': validation.get('reviewed_frames', 0) >= 2
                and validation.get('reviewed_counts', {}).get('SINGLE_BUBBLE', 0) >= 10,
            'automatic_detection_validated': False, 'group_shape_admitted': False,
            'cross_session_generalization': 'NOT_VERIFIED', 'control_permission': False,
            'note': 'Two frames / ten bubbles is a minimal assisted-SAM admission prerequisite, not proof of detection accuracy.'}
