#!/usr/bin/env python3
"""Read and SHA-check the minimal frozen reference for the native calculation."""
import json
import shutil
import sys
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from sv_validation.provenance import sha256, write_json, now, git_state
from sv_validation.validation import reference_condition, SI_UNITS

old = ROOT.parent / 'FEM'
snapshot = json.loads((old / 'inputs/stage03/reference_condition.json').read_text())
source_config = old / 'configs/stage03_reference_vascular.yaml'
assert sha256(source_config) == snapshot['yaml_sha256']
selected = Path(snapshot['config']['mesh']['source'])
assert selected.as_posix() == 'outputs/stage01_7/selected'
sources = {
    'exterior_surface.npz': old / selected / 'surface/tagged_surface_si.npz',
    'volume_mesh.npz': old / selected / 'mesh/volume_mesh.npz',
    'port_contract.json': old / selected / 'planar_port_contract_v2.json',
    'stage03_reference_vascular.yaml': source_config,
    'stage03_reference_condition.json': old / 'inputs/stage03/reference_condition.json',
}
records = {}
for name, source in sources.items():
    target = ROOT / 'inputs/fem_reference' / name
    digest = sha256(source)
    if target.exists():
        assert sha256(target) == digest
    else:
        shutil.copy2(source, target)
    records[name] = {'source_path': str(source), 'source_sha256': digest,
                     'copied_path': str(target.relative_to(ROOT)), 'copied_sha256': sha256(target),
                     'size_bytes': target.stat().st_size, 'verified_at': now()}
assert records['volume_mesh.npz']['source_sha256'] == snapshot['config']['mesh']['volume_mesh_sha256']
write_json(ROOT / 'inputs/MANIFEST.json', {'status': 'PASS', 'source_git': git_state(old), 'files': records})
frozen = yaml.safe_load((ROOT / 'inputs/fem_reference/stage03_reference_vascular.yaml').read_text())
assert frozen == snapshot['config']
previous = {key: frozen[key] for key in ('schema_version', 'condition_type', 'experimental', 'warning', 'physics')}
previous['units'] = dict(SI_UNITS)
previous['source'] = {'artifact': 'inputs/fem_reference/stage03_reference_vascular.yaml',
                      'original_sha256': snapshot['yaml_sha256']}
reference_condition(previous)
previous['source']['manifest'] = 'inputs/MANIFEST.json'
previous['current_lbm_configuration_read'] = False
previous['condition'] = 'REFERENCE_NUMERICAL_CONDITION'
(ROOT / 'configs/sv_reference.yaml').write_text(yaml.safe_dump(previous, sort_keys=False))
(ROOT / 'configs/reference_condition.yaml').unlink(missing_ok=True)
write_json(ROOT / 'reports/sv1/reference_physics.json', previous)
print('Formal inputs PASS:', len(records), 'files; frozen reference:', previous['physics'])
