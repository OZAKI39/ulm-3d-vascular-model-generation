#!/usr/bin/env python3
"""Export a size-aware project snapshot; never delete source files or push Git.

Only the destination is written. Large regenerable artifacts are inventoried,
while current deliverables, source, figures and verification records survive.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
from collections import Counter

MIB = 1024 ** 2
TEXT = {'.py', '.sh', '.md', '.txt', '.json', '.yaml', '.yml', '.csv', '.tsv',
        '.log', '.xml', '.html', '.patch', '.rev', '.cpp', '.cc', '.c', '.h',
        '.hpp', '.cmake', '.in', '.toml', '.ini', '.cfg', '.conf', '.rst',
        '.js', '.css', '.m', '.swc', '.graphml', '.svg'}
FIGURES = {'.png', '.jpg', '.jpeg', '.pdf', '.gif', '.mp4', '.webm'}
GEOMETRY = {'.stl', '.vtp', '.vtk', '.3mf', '.3dm', '.step', '.stp', '.ply', '.obj'}
SKIP_PARTS = {'.git', '.venv', '__pycache__', '.pytest_cache', '.ruff_cache',
              'node_modules', 'build-wsl'}
CURRENT_ZIPS = {
    'P1S_04_ABS_PLA_unsliced_3MF_GUI_verified.zip',
    'Straight_2Top_P1S04_ABS_fine_vertical_unsliced.zip',
    'Straight_2Top_P1S04_ABS_orientation_options_unsliced.zip',
    'candidate_0.zip', 'candidate_15.zip',
}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(MIB), b''):
            h.update(block)
    return h.hexdigest()


def key_geometry(path):
    s = path.as_posix()
    return (s.startswith('outputs/straight_channel/Straight_2Top_core_')
            or 'ABScore_PLAfixture' in s
            or '/compact_manufacturing_roi/candidates/' in s
            or '/print_fixture_design_all_ports_aligned/' in s
            or '/candidate_0_and_15/' in s
            or any(p in path.parts for p in ('selected', 'final', 'source', 'vascular_core'))
            or ('/manufacturing_roi/' in s and '/orientation/' not in s and '/BambuStudio/' not in s)
            or '/refined_roi/' in s
            or ('/all_abs_sacrificial_casting_system/' in s and '/bambu/' not in s)
            or '/support_aware_roi_60x50/casting_mold' in s)


def decision(path, size, symlink=False):
    s = path.as_posix(); parts = path.parts; ext = path.suffix.lower()
    if any(p in SKIP_PARTS for p in parts):
        return False, 'runtime_or_cache'
    if s.startswith('tools/external/') or s.startswith('outputs/topbrain_brava_transfer/runtime/'):
        return False, 'downloadable_runtime'
    if s.startswith('sync/20261022/') or s.startswith('logs/github_sync_20261022/'):
        return False, 'export_audit_generated_separately'
    if s == 'cleanup_inventory.json':
        return False, 'large_historical_cleanup_inventory'
    if path.name.endswith(':Zone.Identifier'):
        return False, 'windows_download_metadata'
    if path.name in {'.env', 'server_runtime.env', 'id_rsa', 'id_ed25519', 'credentials.json'} or ext in {'.pem', '.key'}:
        return False, 'machine_local_configuration'
    if symlink:
        if s == 'data/TopBrain' or s.startswith('deliverables/'):
            return True, 'relative_project_link'
        return False, 'machine_specific_or_legacy_link'
    if s.startswith('data/external/TopBrain/'):
        if ext in {'.txt', '.json'} or path.name == 'License.txt':
            return True, 'dataset_provenance_license_labelmap'
        return False, 'redownloadable_topbrain_volume_or_archive'
    if s.startswith('vessel_model/'):
        if s.startswith('vessel_model/T - Brava/'):
            return size <= 45 * MIB, 'brava_source_and_radius_validation'
        if 'mouse brain vasculature_vmd_branch_aic/' in s and not symlink:
            return size <= 10 * MIB, 'mouse_aic_source_and_validation'
        if 'mouse brain vasculature/' in s and path.name.startswith('fMOST_0_5_6_0_0_6_0001_02_01.'):
            return True, 'configured_mouse_example'
        return False, 'unused_bulk_reference_dataset'
    if s.startswith('third_party/vascularmd/Output/'):
        return False, 'vendor_example_generated_output'
    if size > 45 * MIB:
        return False, 'oversize_regenerable_file'
    if parts[0] in {'tests', 'test_data'}:
        return True, 'regression_fixture_or_test_source'
    if ext in {'.gcode', '.bgcode', '.pkl', '.pickle', '.npy'}:
        return False, 'regenerable_toolpath_or_geometry_cache'
    if 'actual_extrusion_paths' in path.name or 'support_toolpath_segments' in path.name:
        return False, 'regenerable_toolpath_cache'
    if ext in {'.exe', '.dll', '.so', '.a', '.o', '.pyc', '.pyd', '.deb', '.appimage', '.conda', '.jar'}:
        return False, 'compiled_runtime'
    if ext in {'.zip', '.gz', '.7z', '.tar'}:
        keep = parts[0] == 'deliverables' and (path.name in CURRENT_ZIPS or 'ABScore_PLAfixture' in path.name)
        return keep, 'selected_delivery_archive' if keep else 'duplicate_or_downloadable_archive'
    if ext in TEXT or path.name in {'LICENSE', 'License.txt', '.gitignore', '.gitattributes', 'Makefile'}:
        return True, 'source_configuration_report_or_log'
    if ext in FIGURES:
        return size <= 20 * MIB, 'figure_or_review_media' if size <= 20 * MIB else 'large_review_media'
    if parts[0] in {'tests', 'test_data', 'data'}:
        return True, 'test_fixture_or_primary_input'
    if ext in GEOMETRY:
        if ext in {'.step', '.stp'} and size <= MIB:
            return True, 'compact_geometry_regression_reference'
        if (parts[0] == 'history' or s.startswith(('outputs/straight_channel/', 'outputs/y_channel/'))) and not key_geometry(path):
            return False, 'superseded_channel_geometry'
        if key_geometry(path) or size <= 2 * MIB:
            return True, 'key_model_or_small_geometry'
        return False, 'intermediate_geometry_above_2MiB'
    if ext == '.npz':
        return size <= 5 * MIB, 'compact_numeric_data' if size <= 5 * MIB else 'large_numeric_cache'
    return size <= MIB, 'small_project_auxiliary' if size <= MIB else 'large_unclassified_auxiliary'


def export(source, destination, audit):
    source = source.resolve(); destination = destination.resolve()
    if source == destination or source in destination.parents:
        raise ValueError('Destination must be a separate checkout outside the source project')
    destination.mkdir(parents=True, exist_ok=True); audit.mkdir(parents=True, exist_ok=True)
    included = []; excluded = []; by_reason = Counter(); bytes_by_reason = Counter()
    for parent, dirs, files in os.walk(source, followlinks=False):
        for d in list(dirs):
            p = Path(parent) / d; rel = p.relative_to(source)
            if p.is_symlink():
                files.append(d); dirs.remove(d)
            elif d in SKIP_PARTS:
                dirs.remove(d); excluded.append([rel.as_posix() + '/', '', 'runtime_or_cache'])
        for name in sorted(files):
            p = Path(parent) / name; rel = p.relative_to(source); st = p.lstat()
            keep, reason = decision(rel, st.st_size, p.is_symlink())
            if not keep:
                excluded.append([rel.as_posix(), st.st_size, reason])
                by_reason[reason] += 1; bytes_by_reason[reason] += st.st_size
                continue
            dest = destination / rel; dest.parent.mkdir(parents=True, exist_ok=True)
            if p.is_symlink():
                target = os.readlink(p)
                if os.path.isabs(target) or not p.resolve().is_relative_to(source):
                    raise ValueError('Non-portable project link: ' + str(rel))
                if dest.exists() or dest.is_symlink():
                    dest.unlink()
                dest.symlink_to(target)
                digest = hashlib.sha256(target.encode()).hexdigest(); kind = 'symlink'
            else:
                digest = sha(p); kind = 'file'
                shutil.copy2(p, dest)
                if sha(dest) != digest:
                    raise ValueError('Copy hash mismatch: ' + str(rel))
            included.append(dict(path=rel.as_posix(), bytes=st.st_size, sha256=digest, kind=kind, reason=reason))
    included.sort(key=lambda r: r['path']); excluded.sort()
    (audit/'included.json').write_text(json.dumps(included, ensure_ascii=False, indent=2) + '\n')
    with (audit/'excluded.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['path', 'bytes', 'reason']); w.writerows(excluded)
    summary = dict(date='2026-10-22', source=str(source), destination=str(destination),
        included_files=len(included), included_bytes=sum(r['bytes'] for r in included),
        excluded_file_records=len(excluded), excluded_bytes=sum(bytes_by_reason.values()),
        exclusion_reasons={k:dict(files=n,bytes=bytes_by_reason[k]) for k,n in by_reason.items()},
        largest_included=max(included, key=lambda r:r['bytes']),
        local_source_files_deleted=0, models_regenerated=False, slicing_performed=False)
    (audit/'selection_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--destination', type=Path, required=True)
    p.add_argument('--audit', type=Path, required=True)
    a = p.parse_args(); export(a.source, a.destination, a.audit)
