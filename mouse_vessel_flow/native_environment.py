"""Check the explicitly installed local native stack without launching a solve."""
import json
from pathlib import Path
from .paths import project_path
from .provenance import sha256


def verified_installation(build=None):
    config = json.loads(project_path('configs/environment.json').read_text())
    path = Path(config['installation_manifest'])
    if sha256(path) != config['installation_manifest_sha256']:
        raise ValueError('Native installation manifest changed; review and version it explicitly')
    manifest = json.loads(path.read_text())
    for name, expected in manifest['pinned_files'].items():
        if sha256(name) != expected:
            raise ValueError('Installed native file changed: ' + name)
    installed_build = json.loads(Path(manifest['build_json']).read_text())
    if build is not None and build != installed_build:
        raise ValueError('Requested build differs from the audited local installation')
    return manifest


def native_status(component):
    manifest = verified_installation()
    receipt = Path(manifest['checks'][component])
    result = json.loads(receipt.read_text())
    if result['status'] != 'PASS':
        raise ValueError('Native small test did not pass: ' + component)
    return {'status': 'PASS', 'scope': 'Installed runtime and recorded small execution only',
            'receipt': str(receipt), 'new_full_mouse_run': False,
            'known_limits': manifest['known_limits']}
