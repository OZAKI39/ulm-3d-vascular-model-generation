from PIL import Image
from sv_validation.provenance import sha256

NAMES="source_geometry simvascular_surface mesh_cutaway mesh_quality real_geometry_and_bc velocity_global velocity_slices pressure_global pressure_sections flux_balance outlet_flow_split steady_convergence solver_resource_usage".split()
def test_actual_figure_artifacts(root, report):
    provenance=report("visual_provenance")
    assert provenance["flow_result_sha256"] == report("flow_qc")["sha256"]
    hashes=[]
    for name in NAMES:
        path=root/"reports/sv1"/(name+".png")
        with Image.open(path) as im:
            assert im.width >= 800 and im.height >= 500
            im.verify()
        hashes.append(sha256(path))
        assert hashes[-1] == provenance["figures"][name+".png"]
    assert len(set(hashes)) == len(NAMES)
