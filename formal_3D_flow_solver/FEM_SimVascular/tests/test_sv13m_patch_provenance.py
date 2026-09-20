from sv13m_support import *
from sv_validation.sv13m import *
def test_provenance_complete_and_bounded():
 d=accepted('repair_iterations');assert len(d['iterations'])<=3
 for r in d['iterations']:
  for k in ('trigger','root_cause','modified_files','modified_functions','added_lines','deleted_lines','tests_added','before_result','after_result'):assert k in r
 text=(ROOT/'patches/sv1_3m/PATCH_PROVENANCE.md').read_text();assert 'NONE — compatibility only' in text and 'Revert:' in text
