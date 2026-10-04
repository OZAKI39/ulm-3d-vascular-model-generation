"""Read-only inventory, hashes and reference analysis for the cleanup audit.

All paths in CSV files are relative to the outer Git root. Symlinks are not
followed. This program cannot delete files; deletion is a separate gated step.
"""
from __future__ import annotations
import ast
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import struct
import subprocess
import zipfile
from urllib.parse import unquote

AUDIT = Path(__file__).resolve().parent
PROJECT = AUDIT.parents[1]
ROOT = PROJECT.parent
ARRAY_EXT = {'.npz', '.npy', '.pkl', '.pickle'}
TEMP_TEST_NAMES = {'baseline_tmp', 'final_tmp', 'all_tmp', 'tmp_baseline', 'tmp_final',
                   'pytest_tmp', 'baseline_tests', 'full_tests', 'targeted_tmp', 'new_tests_first'}
CACHE_NAMES = {'__pycache__', '.pytest_cache', '.mypy_cache', '.ruff_cache',
               '.ipynb_checkpoints', 'htmlcov', 'numba_cache', 'runtime_cache'}


def dump_json(name, data):
    (AUDIT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')


def dump_csv(name, rows, fields):
    with (AUDIT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def relative(path):
    return str(path.relative_to(ROOT))


def test_temp(path):
    p = Path(path)
    return any(part in TEMP_TEST_NAMES for part in p.parts) and (
        '/tests/' in path or any(part.startswith('build_') for part in p.parts))


def inventory():
    tracked = set(subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0'))
    raw = subprocess.check_output(['git', 'status', '--porcelain=v1', '-z', '--untracked-files=all'], cwd=ROOT).decode().split('\0')
    statuses = {line[3:]: line[:2] for line in raw if len(line) >= 4}
    paths = []
    for directory, dirs, files in os.walk(ROOT, followlinks=False):
        if Path(directory) == AUDIT:
            dirs[:] = []
            continue
        for name in files + [d for d in dirs if (Path(directory) / d).is_symlink()]:
            paths.append(Path(directory) / name)
    rels = [relative(p) for p in paths]
    check = subprocess.run(['git', 'check-ignore', '--stdin', '-z'], input='\0'.join(rels).encode() + b'\0', cwd=ROOT, capture_output=True)
    ignored = set(check.stdout.decode().split('\0'))
    rows = []
    for p, name in zip(paths, rels):
        s = p.lstat()
        kind = 'symlink' if stat.S_ISLNK(s.st_mode) else 'file' if stat.S_ISREG(s.st_mode) else 'other'
        rows.append(dict(path=name, type=kind, size_bytes=s.st_size, mtime=datetime.fromtimestamp(s.st_mtime, timezone.utc).isoformat(),
                         mtime_ns=s.st_mtime_ns, git_tracked=name in tracked,
                         git_status=statuses.get(name, '!!' if name in ignored else '  ' if name in tracked else 'UNTRACKED_OR_GIT_METADATA'),
                         extension=p.suffix.lower(), top_level_group=name.split('/')[0],
                         candidate_category='UNCLASSIFIED', reason='', hash='', link_target=os.readlink(p) if kind == 'symlink' else ''))
    return sorted(rows, key=lambda r: r['path'])


def hash_file(row):
    if row['type'] == 'symlink':
        return hashlib.sha256(os.readlink(ROOT / row['path']).encode()).hexdigest()
    if row['type'] != 'file':
        return ''
    h = hashlib.sha256()
    with (ROOT / row['path']).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def npy_header(f):
    if f.read(6) != b'\x93NUMPY':
        raise ValueError('invalid NPY signature')
    major, minor = f.read(2)
    n = struct.unpack('<H' if major == 1 else '<I', f.read(2 if major == 1 else 4))[0]
    header = ast.literal_eval(f.read(n).decode('utf8' if major >= 3 else 'latin1'))
    return {'shape': list(header['shape']), 'dtype': str(header['descr']), 'fortran_order': header['fortran_order']}


def array_metadata(row):
    p = ROOT / row['path']
    item = {'path': row['path'], 'sha256': row['hash'], 'arrays': {}}
    try:
        if p.suffix == '.npz':
            with zipfile.ZipFile(p) as z:
                for name in z.namelist():
                    if name.endswith('.npy'):
                        with z.open(name) as f:
                            item['arrays'][name] = npy_header(f)
        elif p.suffix == '.npy':
            with p.open('rb') as f:
                item['arrays'][p.name] = npy_header(f)
        else:
            item['note'] = 'Pickle hashed only; not deserialized.'
    except Exception as e:
        item['metadata_error'] = repr(e)
    return item


def formal_files(rows):
    for row in rows:
        path = row['path']; p = ROOT / path
        if row['type'] != 'file' or '.git' in p.parts or test_temp(path) or any(part in CACHE_NAMES for part in p.parts):
            continue
        if any(x in p.parts for x in ('external', 'python_deps', 'julia_depot', 'pdf_tools')):
            continue
        if row['extension'] == '.md' or (row['extension'] == '.json' and p.name.startswith(('STATUS', 'SUMMARY'))):
            yield p


def json_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for k, v in value.items():
            if isinstance(k, str) and '/' in k:
                yield k
            yield from json_strings(v)
    elif isinstance(value, list):
        for v in value:
            yield from json_strings(v)


def extract_references(rows):
    known = {str(ROOT / row['path']): row['path'] for row in rows}
    refs, unresolved, nonlocal_refs = [], [], []
    seen = set()
    for source in formal_files(rows):
        content = source.read_text(errors='replace')
        tokens = []
        if source.suffix == '.md':
            tokens += [(m.group(1), 'markdown_link') for m in re.finditer(r'\]\(([^)]+)\)', content)]
            tokens += [(m.group(1), 'inline_path') for m in re.finditer(r'`([^`\n]+)`', content)
                       if '/' in m.group(1) or re.search(r'\.(?:npz|npy|csv|json|log|png|pdf|py|sh|txt)$', m.group(1))]
        else:
            try:
                tokens += [(v, 'json_path') for v in json_strings(json.loads(content))
                           if (('/' in v and not '\n' in v and len(v) < 1024) or re.search(r'\.(?:npz|npy|csv|json|log|png|pdf)$', v))]
            except Exception as e:
                unresolved.append(dict(source=relative(source), token='', kind='json_parse_error', reason=repr(e)))
        for token, kind in tokens:
            token = token.strip().strip('<>')
            if not token or token.startswith(('https://', 'http://', '#', 'mailto:', 'data:')):
                continue
            if kind == 'markdown_link':
                token = re.sub(r'\s+"[^"]*"$', '', token)
            token = unquote(token.split('#')[0])
            token = re.sub(r'(?<=\.py):\d+(?:[–-]\d+)?$', '', token)
            key = (relative(source), token, kind)
            if key in seen:
                continue
            seen.add(key)
            # Code snippets, expressions and glob templates are not literal artifact paths.
            if any(c in token for c in ('\n', '*', '{', '}', '=', '|', ';')) or token.startswith(('$', '--')):
                nonlocal_refs.append(dict(source=relative(source), token=token, kind=kind, reason='expression_or_template'))
                continue
            q = Path(token)
            ancestors = [p for p in source.parents if p == PROJECT or PROJECT in p.parents]
            choices = [q] if q.is_absolute() else [source.parent / q, PROJECT / q, ROOT / q, *[p / q for p in ancestors]]
            matches = []
            for choice in choices:
                c = Path(os.path.abspath(choice))
                if c.exists() or c.is_symlink():
                    if str(c) not in [str(a) for a in matches]:
                        matches.append(c)
            if matches:
                # Prefer the Markdown document-relative link; otherwise preserve every plausible target.
                if kind == 'markdown_link':
                    matches = matches[:1]
                for match in matches:
                    within = match == ROOT or ROOT in match.parents
                    refs.append(dict(source=relative(source), token=token, kind=kind,
                                     target=relative(match) if within else str(match),
                                     within_repository=within, is_directory=match.is_dir(), exists_before=True))
            else:
                reason = ('external_or_historical_absolute_path' if q.is_absolute() else
                          'unresolved_literal_markdown_link' if kind == 'markdown_link' else 'unresolved_contextual_path_or_label')
                unresolved.append(dict(source=relative(source), token=token, kind=kind, reason=reason))
    return refs, unresolved, nonlocal_refs


def source_audit(rows):
    sources = [r for r in rows if r['extension'] == '.py' and r['type'] == 'file'
               and len(Path(r['path']).parts) > 2 and Path(r['path']).parts[1] in ('pd_clot', 'tests', 'scripts')]
    content = {r['path']: (ROOT / r['path']).read_text(errors='replace') for r in sources}
    imports = defaultdict(set)
    syntax_errors = []
    for path, text in content.items():
        try:
            tree = ast.parse(text, filename=path)
        except SyntaxError as e:
            syntax_errors.append({'path': path, 'error': str(e)}); continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for name in node.names:
                    imports[name.name].add(path)
            elif isinstance(node, ast.ImportFrom):
                prefix = node.module or ''
                imports[prefix].add(path)
                for name in node.names:
                    imports[(prefix + '.' + name.name).strip('.')].add(path)
    output = []
    for row in sources:
        path = row['path']; p = Path(path); module = '.'.join(p.with_suffix('').parts[1:])
        consumers = sorted({v for key, values in imports.items() if key == module or key.startswith(module + '.') for v in values})
        textual = [q for q, text in content.items() if q != path and (p.name in text or module in text)]
        signal = bool(re.search(r'deprecated|legacy|obsolete|unused|backup|prototype|scratch|\bv0\b', path + '\n' + content[path], re.I))
        category = 'ACTIVE_CODE' if consumers or textual or p.parts[1] == 'tests' else 'HISTORICAL_CODE_KEEP' if signal else 'AMBIGUOUS'
        output.append(dict(path=path, classification=category, import_consumers=';'.join(consumers),
                           textual_consumers=';'.join(textual), legacy_keyword_signal=signal,
                           deletion_decision='KEEP', reason='Production/tests and historical reproduction sources are protected; absence of static import is not evidence of dead code.'))
    dump_csv('SOURCE_DEPENDENCY_AUDIT.csv', output, ['path', 'classification', 'import_consumers', 'textual_consumers', 'legacy_keyword_signal', 'deletion_decision', 'reason'])
    dump_json('SOURCE_AST_AUDIT.json', {'parsed_files': len(sources), 'syntax_errors': syntax_errors})


def main():
    rows = inventory()
    print('Enumerated', len(rows), 'entries; hashing with 4 threads.', flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for row, sha in zip(rows, pool.map(hash_file, rows)):
            row['hash'] = sha
    dump_json('INVENTORY_RAW.json', rows)
    dirs = defaultdict(lambda: [0, 0])
    for row in rows:
        parent = Path(row['path']).parent
        for directory in [parent, *parent.parents]:
            dirs[str(directory)][0] += 1
            dirs[str(directory)][1] += row['size_bytes']
    dump_csv('DIRECTORY_SIZE_REPORT.csv', [dict(directory=k, file_count=v[0], total_size=v[1]) for k, v in sorted(dirs.items(), key=lambda a: a[1][1], reverse=True)], ['directory', 'file_count', 'total_size'])
    refs, unresolved, expressions = extract_references(rows)
    dump_json('REFERENCE_SCAN_BEFORE.json', {'resolved': refs, 'unresolved': unresolved, 'expressions': expressions,
                                            'scope': 'First-party Markdown plus STATUS*.json and SUMMARY*.json; vendor documentation and synthetic pytest outputs excluded.'})
    with (AUDIT / 'REFERENCED_ARTIFACTS.txt').open('w') as f:
        for r in refs:
            f.write(r['target'] + '\t' + r['source'] + '\t' + r['kind'] + '\n')
    arrays = [r for r in rows if r['extension'] in ARRAY_EXT and r['type'] == 'file']
    with ThreadPoolExecutor(max_workers=4) as pool, (AUDIT / 'ARRAY_METADATA.jsonl').open('w') as f:
        for item in pool.map(array_metadata, arrays):
            f.write(json.dumps(item, ensure_ascii=False) + '\n')
    source_audit(rows)
    print(json.dumps({'files': len(rows), 'bytes': sum(r['size_bytes'] for r in rows),
                      'tracked_bytes': sum(r['size_bytes'] for r in rows if r['git_tracked']),
                      'arrays_hashed': len(arrays), 'resolved_references': len(refs),
                      'unresolved_references_for_review': len(unresolved)}), flush=True)


if __name__ == '__main__':
    main()
