from sv13n_support import *
def test_reference_hashes_unchanged():
 from sv_validation.provenance import sha256
 for f in accepted('reference_manifest')['files']:assert sha256(ROOT/f['path'])==f['sha256']
 assert read('reference_manifest')['CPU_production_read_only']
