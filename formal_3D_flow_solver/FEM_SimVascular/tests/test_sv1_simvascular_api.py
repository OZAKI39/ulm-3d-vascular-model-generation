def test_official_api_actual_probe(report):
    a=report("simvascular_api")
    assert a["status"] == "PASS" and a["exit_code"] == 0
    assert a["result"]["model_type"] == "<class 'modeling.PolyData'>"
    assert a["result"]["actual_mesher_type"] == "<class 'meshing.TetGen'>"
    assert a["result"]["options_surface_mesh_flag"] and a["result"]["options_volume_mesh_flag"]
