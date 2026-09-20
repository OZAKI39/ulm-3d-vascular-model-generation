from sv_validation.sv11 import load,parse_solver_log,ROOT

def test_complete_historical_failure_history():
    recorded=load('baseline_solver_diagnosis')
    run=load('flow_execution','sv1')['runs'][0]
    import json
    dt=json.loads((ROOT/'configs/time_policy.json').read_text())['dt_s']
    actual=parse_solver_log((ROOT/run['log']).read_text(),dt)
    assert len(actual['linear_solves'])==30
    assert actual['failed_linear_solves']==recorded['failed_linear_solves']==28
    assert actual['ill_conditioned_warnings']==recorded['ill_conditioned_warnings']==12
    assert not actual['unparsed_rows'] and not actual['unassigned_warnings']
    assert actual['first_failed_step']==1
    assert actual['linear_solves'][0]['residual_trustworthy'] is False
