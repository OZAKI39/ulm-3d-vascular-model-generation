from pathlib import Path
import json
import pytest
from fem3d.audit import sha256
from fem3d.mesh_input import read_contract,adapt_surface,PORT_TAGS

ROOT=Path(__file__).resolve().parents[1]

def test_frozen_contract_and_both_surface_hashes():
    c,lock=read_contract(ROOT)
    assert sha256(ROOT/"reports/stage00/source_contract.json")==lock["source_contract_sha256"]
    assert {p['name']:p['surface_entity_id'] for p in c['inlets']+c['outlets']}==PORT_TAGS

def test_tagged_surface_and_once_only_si_conversion():
    import numpy as np
    import pyvista as pv
    p,t,tags,m,c=adapt_surface(ROOT)
    raw=pv.read(c['geometry_path'])
    source=np.asarray(raw.points,dtype=np.float64)[raw.faces.reshape(-1,4)[:,1:]]*1e-6
    np.testing.assert_array_equal(p[t],source)
    assert m['input_units']=='um' and m['solver_units']=='m' and m['conversion_factor']==1e-6
    assert set(tags)=={1,2,3,4,5} and len(t)==67262

def test_contract_tamper_fails_immediately(tmp_path):
    (tmp_path/'configs').mkdir()
    path=tmp_path/'tampered.json';path.write_text('{}')
    lock=json.loads((ROOT/'configs/stage01_source_lock.json').read_text())
    lock['source_contract_path']=str(path)
    (tmp_path/'configs/stage01_source_lock.json').write_text(json.dumps(lock))
    with pytest.raises(ValueError,match='hash mismatch'):read_contract(tmp_path)
