from pathlib import Path
import json,shutil
import pytest
from particle_3d.audit import sha256
from particle_3d.particle6_checkpoint import read_checkpoint
from particle_3d.particle6_cases import mixed_provider

def test_checkpoint_is_binary_state_with_global_only_sidecar(restart_cases):
    for d in restart_cases:
        p=Path(d['checkpoint_path']);m=json.loads((p/'checkpoint_manifest.json').read_text());s=json.loads((p/'particle_sidecar.json').read_text())
        for name,digest in m['files'].items():assert sha256(p/name)==digest
        assert len(m['source_commit'])==40 and m['source_sha256']
        assert 'particles' not in s and s['physical_time_s']==.04 and s['particle_step_index']==40
        assert s['rng_states'] is None and len(s['provenance']['frozen_sha256'])==31
        assert (p/'state.restart').read_bytes().startswith(b'LammpS RestartT')

@pytest.mark.parametrize('damage',['binary','sidecar','schema','provenance','viscosity'])
def test_checkpoint_rejects_corruption(restart_cases,tmp_path,provenance,mu,damage):
    src=Path(restart_cases[1]['checkpoint_path']);p=tmp_path/'copy';shutil.copytree(src,p)
    if damage in ['binary','sidecar']:
        f=p/('state.restart' if damage=='binary' else 'particle_sidecar.json');f.write_bytes(f.read_bytes()+b'x')
    elif damage=='schema':
        f=p/'checkpoint_manifest.json';m=json.loads(f.read_text());m['schema_version']='WRONG';f.write_text(json.dumps(m))
    elif damage=='viscosity':mu*=2
    else:provenance=dict(provenance,dependency_commits={})
    with pytest.raises(ValueError):read_checkpoint(p,provenance,mixed_provider,mu)
