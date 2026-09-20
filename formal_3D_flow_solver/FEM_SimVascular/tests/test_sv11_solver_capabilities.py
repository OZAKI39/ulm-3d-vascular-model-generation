from sv_validation.sv11 import load,ROOT
from sv_validation.provenance import sha256

def test_pinned_capabilities_and_actual_source():
    data=load('linear_solver_capabilities')
    assert data['FSILS_supported'] and data['PETSc_supported'] and data['Trilinos_supported']
    assert not data['baseline_build_PETSc']
    for path,evidence in data['sources'].items():
        assert sha256(ROOT/'external/svMultiPhysics'/path)==evidence['sha256']
    for key in ('relative_tolerance','absolute_tolerance','restart'):
        assert data['xml'][key]['effective_for_PETSc'] is False
    assert data['existing_runtime_interface']['available']

def test_runtime_critical_options_are_observed_not_just_requested():
    from sv_validation.sv11_runtime import OPTIONS
    settings=load('petsc_settings')
    assert settings['options']==OPTIONS
    for name in ('petsc_smoke','petsc_short'):
        path=ROOT/'logs/sv1_1'/(name+'.log')
        if not path.exists():continue
        text=path.read_text()
        for item in settings['runtime_confirmed'].values():assert item['text'] in text
    baseline=load('baseline_solver_diagnosis')['configuration']['settings']
    assert float(OPTIONS['-ksp_rtol'])<=float(baseline['Tolerance'])
    assert float(OPTIONS['-ksp_atol'])<=float(baseline['Absolute_tolerance'])
    assert not any(settings['ambient_options_audit']['paths_checked'].values())
