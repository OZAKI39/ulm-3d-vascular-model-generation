from sv_validation.sv11 import load,linear_gate,ROOT
from sv_validation.validation import parse_result_vtu
from sv_validation.provenance import sha256

def test_official_petsc_smoke_actual_fields_and_reasons():
    data=load('petsc_smoke')
    assert data['status']=='PASS'
    assert linear_gate(data['run'])
    for record in data['result_files']:
        path=ROOT/record['path'];assert sha256(path)==record['sha256']
        mesh,u,p=parse_result_vtu(path)
        assert mesh.n_points==record['points'] and mesh.n_cells==record['cells']
