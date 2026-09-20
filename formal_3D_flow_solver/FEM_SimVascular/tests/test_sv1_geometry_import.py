def test_import_and_four_ports(report):
    assert report("geometry_import")["status"] == "PASS"
    g=report("geometry_qc")
    assert g["boundary_roles_complete"]
    assert set(g["ports"]) == {"INLET","OUTLET_01","OUTLET_02","OUTLET_03"}
    for p in g["ports"].values():
        assert p["connected_regions"] == 1 and p["topological_disk"]
        assert not p["orientation_flip"] and p["normal_dot"] > .999
    assert g["wall"]["forward_sample_count"] > 0 and g["wall"]["backward_sample_count"] > 0
