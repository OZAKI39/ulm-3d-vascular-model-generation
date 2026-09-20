from sv12_support import artifact
from sv_validation.sv12 import ROOT,REPORT
from sv_validation.provenance import sha256
def test_figures_provenance_and_gating():
    data=artifact('visuals')
    for row in data['figures']:assert sha256(ROOT/row['path'])==row['sha256']
    fields=['velocity_global.png','velocity_slices.png','pressure_global.png','pressure_sections.png','flux_balance.png','outlet_flow_split.png']
    for name in fields:assert (REPORT/name).exists()==data['accepted_solution_available']
