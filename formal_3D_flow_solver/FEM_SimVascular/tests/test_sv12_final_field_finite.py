from sv12_support import accepted
from sv_validation.sv12 import ROOT
from sv_validation.validation import parse_result_vtu
from sv_validation.provenance import sha256

def test_final_native_fields_finite():
    final=accepted();grid,u,p=parse_result_vtu(ROOT/final['path'])
    assert sha256(ROOT/final['path'])==final['sha256']
    assert final['velocity_finite'] and final['pressure_finite']
