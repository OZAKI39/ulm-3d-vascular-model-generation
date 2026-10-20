"""Verify every file in the server delivery manifest after the WSL transfer."""
from common import *
import datetime


def main():
    assert str(R).startswith('/home/lzy/projects/')
    manifest_path = O/'DELIVERY_MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    mismatches = []
    for relative, expected in manifest['files'].items():
        path = R/relative
        if not path.is_file():
            mismatches.append(dict(path=relative, reason='missing'))
            continue
        if path.stat().st_size != expected['bytes'] or sha(path) != expected['sha256']:
            mismatches.append(dict(path=relative, reason='size_or_hash'))
        if path.is_symlink() != ('symlink' in expected):
            mismatches.append(dict(path=relative, reason='symlink_identity'))
        elif path.is_symlink() and os.readlink(path) != expected['symlink']:
            mismatches.append(dict(path=relative, reason='symlink_target'))
    library = R/manifest['native_library_external_dependency']
    external = dict(path=str(library), included_in_this_delivery=False,
                    execution_environment='Isolated server environment /workspace/y_clot/env.sh',
                    expected_sha256=manifest['native_library_sha256'])
    if library.is_file():
        external.update(local_file_present=True, local_sha256=sha(library))
        external['local_identity_matches'] = external['local_sha256'] == external['expected_sha256']
    else:
        external['local_file_present'] = False
    result = dict(status='PASS' if not mismatches else 'FAIL',
                  utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  files=len(manifest['files']), total_bytes=manifest['total_bytes'],
                  manifest_sha256=sha(manifest_path), mismatches=mismatches,
                  native_external_dependency=external,
                  physics_recomputed=False, GitHub_pushed=False)
    dump(O/'LOCAL_DELIVERY_READBACK.json', result)
    assert not mismatches, mismatches
    print('LOCAL_DELIVERY_READBACK', result, flush=True)


if __name__ == '__main__':
    main()
