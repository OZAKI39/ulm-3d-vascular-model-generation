from sv_validation.provenance import sha256

def test_native_output_provenance(root, report):
    q=report("flow_qc"); run=report("flow_execution")
    assert run["runs"][-1]["exit_code"] == 0
    assert sha256(root/q["path"]) == q["sha256"]
    assert q["native_vs_vtu_max_difference_over_Q"] < 1e-6

def test_real_solver_converged(report):
    q=report("flow_qc")
    assert q["solver_success"], f"FLOW_SOLVE_FAIL: {q['linear_nonconvergence_warnings']} linear warnings"

def test_required_steady_intervals_reached(report):
    s=report("steady_state")
    assert s["status"] == "PASS", f"NOT_REACHED: {s['saved_states']} saved state(s), five consecutive intervals required"
