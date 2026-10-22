"""Original SI 1D/0D resistance equations and lightweight frozen-network regressions."""
import json
import numpy as np
from .paths import project_path
from .network_models.roi_only_hydraulics import load_roi_cache, solve_roi_mixed_bc
from .network_models.roi_boundary_design import regress_h0, design
from .network_models.hydraulic_resistance import ROI_TARGET_Q_M3_S

def roi_cache():
    return load_roi_cache(
        project_path('inputs/network/reports/a_network_1d0d_boundary_v2_idealized/data/analysis_A_H0_graph_si.npz'),
        project_path('inputs/network/reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json'))

def validate_design():
    """Recompute only the small resistance network, never CFD or particle dynamics."""
    cache=roi_cache()
    h0=json.loads(project_path('inputs/network/reports/a_network_1d0d_boundary_v2_idealized/data/roi_operating_point_H0.json').read_text())
    regression=regress_h0(cache,h0)
    cert=json.loads(project_path('inputs/network/reports/balanced_three_outlet_forward_v1/feasibility_certificate.json').read_text())
    solution,basis,trace=design(cache,regression,cert)
    original=json.loads(project_path('results/roi_only/network/roi_equal_split_solution.json').read_text())
    comparisons={}
    for key in ('relative_realcut_pressure_pa','gauge_realcut_pressure_pa','outlet_flow_m3_s','outlet_fraction'):
        a=np.array(list(solution[key].values()));b=np.array(list(original[key].values()))
        np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-24)
        comparisons[key]={'max_abs_difference':float(np.max(abs(a-b)))}
    assert len(cache.node_ids)==109 and len(cache.edge_nodes)==108
    return dict(status='PASS',nodes=len(cache.node_ids),edges=len(cache.edge_nodes),cache=cache.summary(),
        h0_regression=regression,comparisons=comparisons,mass_audit=solution['mass_audit'],
        pressure_basis=basis,CFD_feedback_used=False)
