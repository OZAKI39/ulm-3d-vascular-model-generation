import numpy as np
from network_1d0d.hydraulic_resistance import ROI_TARGET_Q_M3_S

def test_positive_signed_scaling_and_separate_total_source(domain,baseline):
    assert baseline['unit_roi_q'][0]>0
    scale=ROI_TARGET_Q_M3_S/baseline['unit_roi_q'][0]
    assert baseline['scaling_lambda']==scale
    assert baseline['ports'][0]['signed_Q_m3s']==ROI_TARGET_Q_M3_S
    np.testing.assert_allclose(baseline['pressure'],baseline['unit'].pressure*scale,rtol=0,atol=0)
    np.testing.assert_allclose(baseline['flow'],baseline['unit'].flow*scale,rtol=0,atol=0)
    assert baseline['mass']['A_total_source_inflow']>10*ROI_TARGET_Q_M3_S
    assert baseline['mass']['global_relative_residual']<1e-10
    assert baseline['mass']['ROI_relative_residual']<1e-10
