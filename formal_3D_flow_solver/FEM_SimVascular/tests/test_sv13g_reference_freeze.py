from sv13g_support import *
from sv_validation.provenance import sha256
def test_frozen_reference_bytes():
    d=load('reference_manifest')
    assert d['status']=='PASS' and d['CPU_production_read_only']
    for f in d['files']:assert sha256(ROOT/f['path'])==f['sha256'],f['path']
