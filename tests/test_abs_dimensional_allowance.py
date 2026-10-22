"""Finished minimum vs. actual CAD gap, plus duplicate slicer compensation."""
import copy
import numpy as np
import pytest
import trimesh
from vascular_processing import y_channel_casting as y
from vascular_processing.abs_dimensional_allowance import dimension_budget, minimum_cad, slicer_scale_qc


def config():
    return y.load_config(y.ROOT/'config/history/before_bottomless_dovetail_20261012/y_channel_p1s_abs.yaml')


def test_nominal_nine_mm_fails_finished_requirement_from_actual_mesh():
    cfg = config()
    mold = trimesh.boolean.difference([y.box([[-10,-6.5,0],[10,6.5,13.4]]),
        y.box([[-8,-4.5,2],[8,4.5,14.4]])], engine='manifold')
    inner = [[-8,-4.5],[8,4.5]]
    for result in [y.pdms_width_qc(mold,inner,cfg,11), y.pdms_height_qc(mold,inner,cfg,11)]:
        assert not result['passed']
        assert result['allowance']['budgeted_finished_lower_bound_mm'] < 9


def test_reserve_is_total_dimension_and_larger_targets_round_up():
    cfg = config()
    for dimension in ['width', 'height']:
        assert minimum_cad(cfg,dimension) == 9.6
        result = dimension_budget(cfg,dimension,9.6)
        assert result['budgeted_finished_lower_bound_mm'] == pytest.approx(9.108)
        assert not result['physical_finished_dimensions_validated']
        larger=copy.deepcopy(cfg);larger['minimum_pdms_'+dimension+'_mm']=12
        assert minimum_cad(larger,dimension) == 12.6
    cfg['abs_dimensional_allowance']['width_shrink_fraction'] = .01
    assert minimum_cad(cfg,'width') == 9.4
    assert minimum_cad(cfg,'height') == 9.6
    bounds=y.enforce_pdms_width([[-10,-12],[10,12]],cfg)
    assert np.ptp(bounds[:,1]) == 24


@pytest.mark.parametrize('key,value', [('width_shrink_fraction',float('nan')),
    ('height_shrink_fraction',1.), ('width_error_reserve_mm',-.3), ('cad_round_up_mm',0)])
def test_invalid_budget_rejected(key,value):
    cfg=config();cfg['abs_dimensional_allowance'][key]=value
    with pytest.raises(ValueError,match='INVALID_ABS_DIMENSIONAL_ALLOWANCE'):
        y.validate_pdms_width_requirement(cfg)


def test_double_slicer_compensation_cannot_pass():
    cfg=config()
    settings=dict(filament_shrink=['100%'],xy_contour_compensation='0',
        xy_hole_compensation='0',enable_circle_compensation='0')
    assert slicer_scale_qc(cfg,settings)['passed']
    for key,value in [('filament_shrink',['98%']),('xy_hole_compensation','.2'),
                      ('xy_contour_compensation','.2'),('enable_circle_compensation','1')]:
        assert not slicer_scale_qc(cfg,{**settings,key:value})['passed']
    assert not slicer_scale_qc(cfg,{})['passed']
