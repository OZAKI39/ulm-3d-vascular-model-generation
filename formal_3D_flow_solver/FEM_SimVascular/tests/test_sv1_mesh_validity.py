from sv_validation.provenance import sha256

def test_actual_mesh_validity(root, report):
    m=report("mesh_validity")
    assert m["status"] == "PASS" and m["connected_volumes"] == 1
    assert m["inverted"] == m["degenerate"] == m["nonfinite"] == 0
    assert m["exterior_closed"] and m["boundary_roles_complete"]
    assert sha256(root/"outputs/sv1/SV_MESH/mesh-complete.mesh.vtu") == m["solver_mesh_sha256"]
