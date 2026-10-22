import csv
import numpy as np
from mouse_vessel_flow.paths import project_path
DATA=project_path('inputs/network/reports/a_network_1d0d_boundary_v2_idealized/data')
from mouse_vessel_flow.network_models.idealized_h0 import solve_operating_point

def test_all_saved_leave_one_out_variants_reproduce(domain):
    rows=list(csv.DictReader((DATA/'terminal_leave_one_out_sensitivity.csv').open()))
    assert len(rows)==122
    assert {int(r['removed_terminal_id']) for r in rows}==set(domain.ids[domain.terminals])
    for row in rows:
        removed=int(np.flatnonzero(domain.ids==int(row['removed_terminal_id']))[0])
        solved=solve_operating_point(domain,removed_terminal=removed)
        assert solved['status']==row['status']=='PASS'
        assert removed not in solved['terminals'] and len(solved['terminals'])==121
        assert abs(solved['node_outflow'][removed])<1e-10*solved['mass']['A_total_source_inflow']
        np.testing.assert_allclose(solved['fractions'],[float(row[n+'_fraction']) for n in ['O1','O2','O3']],rtol=0,atol=1e-10)
