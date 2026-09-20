from sv13_support import artifact
from sv_validation.sv13 import equivalence,frozen_policy
from sv_validation.postprocess import SolutionMeasurements
from sv13_support import ROOT
def test_selected_field_independently_equivalent():
    d=artifact('production_selection');ref=artifact('reference_freeze');p=frozen_policy()
    m=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',ref['accepted_solution']['Q_target_m3_s'],p['Umean_m_s'])
    e=equivalence(m,ROOT/d['accepted_field'],ROOT/ref['accepted_solution']['path'],p)
    assert e['status']=='PASS'
