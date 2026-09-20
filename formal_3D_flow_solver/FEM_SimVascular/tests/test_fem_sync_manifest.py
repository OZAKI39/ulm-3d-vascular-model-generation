import csv
import importlib.util
from fem_sync_support import ROOT, read, digest


def test_every_published_file_and_hash():
    spec = importlib.util.spec_from_file_location('verify_manifest',ROOT/'scripts/fem_freeze_sync/verify_manifest.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.verify() > 1000


def test_all_copied_source_and_tests_preserve_original_hashes():
    copied=read('sync_metadata/copied_from_wsl.json')
    for row in copied['files']:
        assert digest(ROOT/row['relative_path']) == row['sha256'], row['relative_path']
    import hashlib
    for row in read('sync_metadata/upstream_lfs_pointer_inventory.json')['files']:
        assert hashlib.sha256(row['original_pointer_text'].encode()).hexdigest()==row['sha256']
    assert any(row['source_path'].startswith('tests/test_sv13q_') for row in copied['files'])
    assert read('sync_metadata/patch_reconstruction.json')['status']=='PASS'


def test_no_secret_or_oversized_payload():
    assert read('sync_metadata/secret_scan.json')['status']=='PASS'
    assert read('sync_metadata/secret_scan.json')['findings']==[]
    with (ROOT/'sync_metadata/source_manifest.csv').open() as f:
        for row in csv.DictReader(f):
            assert int(row['size_bytes']) < 50*2**20
            assert row['category'] in {'SOURCE','CONFIG','TEST','SCRIPT','REPORT','REFERENCE','LOG','RESULT','DOC'}
            assert not set(row['relative_path'].split('/')) & {'external','build','CMakeFiles','.venv','__pycache__'}
