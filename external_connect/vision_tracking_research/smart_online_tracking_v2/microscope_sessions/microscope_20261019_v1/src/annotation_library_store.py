"""Saved annotation inventory and reversible removal of active revision pointers.

Only <frame>/latest.json is consumed by annotation loaders and training. Moving
that pointer into deleted_records is atomic; immutable masks/history stay put.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import uuid
import numpy as np
from annotation_store import CLASSES, sha, validate, save_annotation
from annotation_regions import valid_regions


def frame_folder(directory, record):
    name = record['id']
    if not name or Path(name).name != name or name in ('.', '..'):
        raise ValueError('图片 ID 不是合法目录名')
    return Path(directory) / name


def read_revision(folder, pointer, record):
    meta = json.loads(pointer.read_text())
    if meta['candidate_id'] != record['id'] or meta['image_sha256'] != record['image_sha256']:
        raise ValueError('标注与原图身份不一致')
    name = meta['mask_file']
    if Path(name).name != name:
        raise ValueError('标注文件路径无效')
    path = folder / name
    if sha(path) != meta['mask_sha256']:
        raise ValueError('标注 SHA-256 校验失败')
    with np.load(path, allow_pickle=False) as z:
        masks = {k: z[k].copy() for k in CLASSES}
    validate(masks, (record['resolution_wh'][1], record['resolution_wh'][0]), meta['expert_confirmed'])
    if len(valid_regions(meta, masks)) != len(meta.get('review_regions', [])):
        raise ValueError('审核区域校验失败')
    return meta, masks


def saved_inventory(directory, records, deleted=False):
    """Metadata-only listing; full mask integrity is checked before mutations."""
    rows = []
    for index, record in enumerate(records):
        folder = frame_folder(directory, record)
        pointers = sorted((folder / 'deleted_records').glob('*.json'), reverse=True) if deleted else [folder / 'latest.json']
        for pointer in pointers:
            if not pointer.exists() or pointer.name.endswith('.restored.json'):
                continue
            row = {'index': index, 'record': record, 'pointer': pointer, 'deleted': deleted,
                   'key': (record['id'], pointer.name), 'meta': {}, 'error': ''}
            try:
                row['meta'] = json.loads(pointer.read_text())
                if row['meta']['candidate_id'] != record['id'] or row['meta']['image_sha256'] != record['image_sha256']:
                    raise ValueError('标注与原图身份不一致')
                for key in ('revision', 'utc', 'status', 'instance_counts'):
                    if key not in row['meta']:
                        raise ValueError('缺少字段：' + key)
            except Exception as exc:
                row['error'] = str(exc)
            rows.append(row)
    return rows


def delete_saved(directory, record, expected_revision):
    folder = frame_folder(directory, record)
    pointer = folder / 'latest.json'
    meta, _ = read_revision(folder, pointer, record)
    if meta['revision'] != expected_revision:
        raise ValueError('保存记录已变化，请刷新列表后重试')
    trash = folder / 'deleted_records'
    trash.mkdir(exist_ok=True)
    name = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ') + '_' + uuid.uuid4().hex[:8] + '.json'
    destination = trash / name
    pointer.replace(destination)
    return destination


def restore_saved(directory, record, archived_pointer, reviewer=''):
    folder = frame_folder(directory, record)
    pointer = Path(archived_pointer)
    if pointer.parent.resolve() != (folder / 'deleted_records').resolve() or pointer.name.endswith('.restored.json'):
        raise ValueError('已删除记录路径无效')
    if (folder / 'latest.json').exists():
        raise ValueError('此图片已有保存记录；为避免覆盖，请先编辑现有记录')
    meta, masks = read_revision(folder, pointer, record)
    # Restoring never silently reinstates expert approval or reviewed regions.
    restored = save_annotation(directory, record, masks, reviewer, False, meta.get('track_ids'), meta.get('notes', ''),
                               review_scope=meta.get('review_scope', 'REGIONS'),
                               workspace_roi=meta.get('workspace_roi_yxyx'),
                               proposal_sources=meta.get('proposal_sources'))
    pointer.rename(pointer.with_suffix('.restored.json'))
    return restored
