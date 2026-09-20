from sv13_support import artifact
def test_cpu_equivalent_to_accepted_step400():
    e=artifact('cpu_reference_equivalence')
    assert e['status']=='PASS' and all(e['checks'].values())
    assert e['reference_sha256']==artifact('reference_freeze')['accepted_solution']['sha256']
def test_fresh_process_reload():
    r=artifact('solution_reload')
    assert r['status']=='PASS' and r['fresh_process'] and all(r['checks'].values())
def test_perturbed_native_field_fails_equivalence(tmp_path):
    import pyvista as pv
    from sv13_support import ROOT
    from sv_validation.sv13 import equivalence,frozen_policy
    from sv_validation.postprocess import SolutionMeasurements
    reference=artifact('reference_freeze')['accepted_solution'];policy=frozen_policy()
    source=ROOT/reference['path'];grid=pv.read(source)
    velocity=next(k for k in grid.point_data if k.lower()=='velocity')
    grid.point_data[velocity]=grid.point_data[velocity]*1.01
    path=tmp_path/'perturbed.vtu';grid.save(path)
    measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',reference['Q_target_m3_s'],policy['Umean_m_s'])
    comparison=equivalence(measure,path,source,policy)
    assert comparison['status']=='FAIL'
    assert not comparison['checks']['velocity_relative_volume_L2']
    assert not comparison['checks']['relative_Qin']
