from sv13m_support import *
from sv_validation.sv13m import *
def test_official_history_is_hashed():
 from sv_validation.provenance import sha256
 d=accepted('upstream_audit');assert d['official_repository']=='https://gitlab.com/petsc/petsc.git'
 assert d['first_fixed_commit']=='8ba8e283c7d3ba8c2bd73cc674b79590a073ab5a' and d['first_fixed_tag']=='v3.25.0'
 for c in d['commands']:assert c['exit_code']==0 and sha256(ROOT/c['artifact'])==c['sha256']
def test_local_adaptation_not_misrepresented_as_upstream():
 d=read('repair_03');assert d['upstream_exact_commit'] is None and 'NOT claimed' in d['provenance']
