from sv_validation.sv11 import load,ROOT
from sv_validation.provenance import sha256

def test_figures_have_actual_evidence_and_honest_scope():
    data=load('visuals')
    assert data['files']
    for f in data['files']:
        assert sha256(ROOT/f['path'])==f['sha256']
    if not data['accepted_solution_available']:
        for name in ('velocity_global.png','velocity_slices.png','pressure_global.png','pressure_sections.png','outlet_flow_split.png'):
            assert not (ROOT/'reports/sv1_1'/name).exists()
