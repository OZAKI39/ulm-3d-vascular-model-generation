from sv13l_support import *
def test_frozen_references_match():
    d=load('reference_manifest')
    assert d['CUDA_winner']=='12.3.2'
    for f in d['files']:assert sha256(ROOT/f['path'])==f['sha256']
def test_stage_j_delivery_remains_frozen():
    d=json.loads((ROOT/'reports/sv1_3j/delivery_manifest.json').read_text())
    for f in d['files']:assert sha256(ROOT/f['path'])==f['sha256']
