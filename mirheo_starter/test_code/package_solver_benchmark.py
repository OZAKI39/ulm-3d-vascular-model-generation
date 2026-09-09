"""Package this completed campaign only; no solver, network or build imports."""
from pathlib import Path
from hashlib import sha256
import json
import zipfile


ROOT = Path(__file__).resolve().parents[1]
HEMO = ROOT.parent / 'hemocell_starter'
CAMPAIGN = 'solver_benchmark_20260909'


def read(path):
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    return json.loads(path.read_text(encoding='utf-8'), parse_constant=invalid)


def digest(path):
    with path.open('rb') as stream:
        result = sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def write_new(path, value):
    data = json.dumps(value, ensure_ascii=False, allow_nan=False, indent=2) + '\n'
    if path.exists():
        if path.read_text(encoding='utf-8') != data:
            raise ValueError('Preserve existing delivery: ' + str(path))
    else:
        with path.open('x', encoding='utf-8') as stream:
            stream.write(data)


def main():
    latest = read(ROOT / 'data/solver_benchmark/LATEST.json')
    page = Path(latest['html'])
    dest = page.parent
    results = read(dest / 'results.json')
    verification = read(dest / 'delivery_verification.json')
    if verification['status'] != 'PASS' or digest(page) != latest['html_sha256']:
        raise ValueError('Delivery verification or HTML hash failed')
    browser = read(dest / 'browser_local_profile/browser_check.json')
    if browser['status'] != 'PASS' or browser['html_sha256'] != digest(page):
        raise ValueError('Browser evidence must match the exact HTML')
    files = {}

    def add(path, archive_name):
        if not path.is_file() or path.is_symlink():
            raise ValueError('Expected a regular file: ' + str(path))
        if archive_name in files or '..' in Path(archive_name).parts:
            raise ValueError('Duplicate or unsafe archive path: ' + archive_name)
        files[archive_name] = path

    def add_tree(path, archive_base):
        for entry in sorted(path.rglob('*')):
            if entry.is_file() and '__pycache__' not in entry.parts:
                add(entry, archive_base + '/' + entry.relative_to(path).as_posix())

    for name in ['benchmark_review.html', 'results.json', 'artifact_sha256.json',
                 '结果说明.md', 'delivery_verification.json']:
        add(dest / name, name)
    for folder in ['browser_local_profile', 'browser_final']:
        for entry in sorted((dest / folder).glob('*')):
            if entry.is_file():
                add(entry, 'browser/' + folder + '/' + entry.name)
    for folder in ['runs/solver_benchmark/' + CAMPAIGN,
                   'data/solver_benchmark/' + CAMPAIGN, 'py_scripts/solver_benchmark']:
        add_tree(ROOT / folder, 'mirheo_starter/' + folder)
    for name in ['py_scripts/benchmark_mirheo_hemocell.py', 'py_scripts/solver_benchmark.yaml',
                 'test_code/test_solver_benchmark.py', 'test_code/check_solver_benchmark_mpi.py',
                 'test_code/check_solver_benchmark_browser.cjs', 'test_code/review_solver_benchmark.py',
                 'test_code/package_solver_benchmark.py', 'test_code/README_solver_benchmark.md']:
        add(ROOT / name, 'mirheo_starter/' + name)
    for folder in ['cases/pure_fluid_benchmark', 'scripts', 'metadata', 'logs']:
        add_tree(HEMO / folder, 'hemocell_starter/' + folder)
    add(ROOT / '.venv/lib/python3.12/site-packages/plotly-6.5.2.dist-info/licenses/LICENSE.txt',
        'licenses/plotly_python_LICENSE.txt')
    add(ROOT / '.venv/lib/python3.12/site-packages/Mirheo-1.6.2.dist-info/LICENSE',
        'licenses/Mirheo_LICENSE.txt')

    entries = {}
    json_count = 0
    for name, path in sorted(files.items()):
        if path.suffix == '.json':
            read(path)
            json_count += 1
        entries[name] = {'sha256': digest(path), 'bytes': path.stat().st_size}
    manifest = {
        'campaign_id': CAMPAIGN, 'file_count': len(entries), 'files': entries,
        'scope': 'Only this campaign, new code, installation evidence and offline review; no old result trees, binaries or virtual environment.',
        'reproduction': 'Offline HTML is self-contained. Solver reruns use the existing workspace modules and pinned installed dependencies; the ZIP is evidence, not a standalone execution environment.',
        'human_review': 'PENDING',
    }
    manifest_path = dest / 'bundle_manifest.json'
    write_new(manifest_path, manifest)
    archive = dest / 'solver_benchmark_results.zip'
    if not archive.exists():
        with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
            for name, path in sorted(files.items()):
                bundle.write(path, name)
            bundle.write(manifest_path, 'bundle_manifest.json')
    with zipfile.ZipFile(archive) as bundle:
        if bundle.testzip() is not None or len(bundle.namelist()) != len(entries) + 1:
            raise ValueError('ZIP CRC or entry count failed')
        for name, entry in entries.items():
            if sha256(bundle.read(name)).hexdigest() != entry['sha256']:
                raise ValueError('ZIP content hash failed: ' + name)
        if bundle.read('bundle_manifest.json') != manifest_path.read_bytes():
            raise ValueError('ZIP manifest differs')
    record = {
        'status': 'DELIVERED_WITH_SCIENTIFIC_LIMITATIONS', 'campaign_id': CAMPAIGN,
        'html': str(page), 'html_sha256': digest(page),
        'archive': str(archive), 'archive_sha256': digest(archive), 'archive_bytes': archive.stat().st_size,
        'manifest': str(manifest_path), 'manifest_sha256': digest(manifest_path),
        'archive_entries': len(entries) + 1, 'strict_json_files_checked': json_count,
        'zip_crc_and_every_file_sha256': 'PASS', 'verification': str(dest / 'delivery_verification.json'),
        'scientific_status': results['comparison']['status'],
        'qualified_speedup': None, 'human_review': 'PENDING',
        'cpu_charged_s': results['cpu_charged_s'], 'gpu_charged_s': results['gpu_charged_s'],
        'new_gpu_600s_remaining_s': 600 - results['gpu_charged_s'],
        'view_command': 'cd /home/lzy/projects/mirheo_starter && .venv/bin/python -B -m test_code.review_solver_benchmark --open',
    }
    write_new(dest / 'delivery_record.json', record)
    print(json.dumps(record, ensure_ascii=False, allow_nan=False, indent=2))


if __name__ == '__main__':
    main()
