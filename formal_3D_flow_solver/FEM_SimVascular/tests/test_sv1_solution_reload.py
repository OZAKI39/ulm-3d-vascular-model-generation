import os

def test_new_process_reload_agrees(report, measured):
    r=report("solution_reload")
    assert r["status"] == "PASS" and r["fresh_process"] and r["process_id"] != os.getpid()
    q=measured[3]
    assert abs(r["measurements"]["epsilon_mass"]-q["epsilon_mass"]) < 1e-12
    assert r["result_sha256"] == report("flow_qc")["sha256"]
