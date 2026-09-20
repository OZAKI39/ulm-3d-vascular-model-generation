from sv13m_support import *
from sv_validation.sv13m import *
def test_references_unchanged():
 from sv_validation.provenance import sha256
 d=accepted('reference_manifest')
 assert d['stack']['svmp_gpu_build']['source_unmodified']
 for f in d['files']:assert sha256(ROOT/f['path'])==f['sha256']
 assert d['MPI_prefix']==d['stack']['compatibility_winner']['MPI_prefix']
