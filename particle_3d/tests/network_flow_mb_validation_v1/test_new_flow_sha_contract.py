import hashlib,json
from pathlib import Path
def test_new_flow_sha_contract(report):
    c=json.loads((report/'data/old_new_flow_contract.json').read_text())
    approved=json.loads(Path(c['authoritative_new_source']).read_text())
    for label in ['OLD','NEW']:
        path=Path(c['paths'][label]);h=hashlib.sha256(path.read_bytes()).hexdigest()
        assert h==c['sha256'][label]
        assert h==hashlib.sha256((report/'server_bundle/inputs'/f'{label}.vtu').read_bytes()).hexdigest()
    assert approved['export']['files'][Path(c['paths']['NEW']).name]['sha256']==c['sha256']['NEW']
