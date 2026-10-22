"""Local cavity allowances, not a global scale or a measured material law.

Finished PDMS minima stay separate from CAD dimensions. ABS cooling shrinkage
and a total dimensional error reserve are budgeted independently for Y and Z.
Archived configurations without this policy retain their original geometry.
"""
import math

SOURCES = [
    'https://help.prusa3d.com/article/abs_2058',
    'https://store.bblcdn.com/s7/default/23b4cf2b83d5470bb96d19970b5f3ae8/Bambu_ABS_Technical_Data_Sheet_V3.pdf',
]


def dimension_budget(cfg, dimension, measured_mm=None):
    if dimension not in ('width', 'height'):
        raise ValueError('UNKNOWN_PDMS_DIMENSION')
    target = float(cfg.get('minimum_pdms_'+dimension+'_mm', 0.))
    if not math.isfinite(target) or target < 0:
        raise ValueError('INVALID_PDMS_MINIMUM')
    if reserve_policy := cfg.get('cad_dimensional_allowance'):
        if cfg.get('abs_dimensional_allowance') or reserve_policy.get('calibration_status')!='UNCALIBRATED_DESIGN_RESERVE':
            raise ValueError('AMBIGUOUS_DIMENSIONAL_ALLOWANCE')
        required=float(reserve_policy[dimension+'_cad_minimum_mm'])
        reserve=float(reserve_policy[dimension+'_error_reserve_mm'])
        if not all(math.isfinite(v) for v in (required,reserve)) or reserve<0 or required-reserve<target:
            raise ValueError('INVALID_PRESERVED_CAD_RESERVE')
        actual=required if measured_mm is None else float(measured_mm)
        return dict(dimension=dimension,finished_minimum_mm=target,required_cad_minimum_mm=required,
            measured_cad_mm=actual,total_dimension_error_reserve_mm=reserve,
            budgeted_finished_lower_bound_mm=actual-reserve,material=reserve_policy['material'],
            assumed_material_shrink_fraction=None,calibration_status=reserve_policy['calibration_status'],
            provenance=reserve_policy['provenance'],sources=[],physical_finished_dimensions_validated=False,
            pdms_cure_shrinkage_characterized=False,
            passed=math.isfinite(actual) and actual>=required-1e-5 and actual-reserve>=target-1e-5)
    policy = cfg.get('abs_dimensional_allowance')
    shrink = reserve = 0.
    step = None
    required = target
    if policy is not None:
        if policy.get('material') != 'ABS' or policy.get('calibration_status') != 'UNCALIBRATED_DESIGN_BUDGET':
            raise ValueError('ABS_ALLOWANCE_POLICY_REQUIRED')
        if target < 9:
            raise ValueError('ABS_ALLOWANCE_REQUIRES_FINISHED_MINIMUM_9_MM')
        shrink = float(policy[dimension+'_shrink_fraction'])
        reserve = float(policy[dimension+'_error_reserve_mm'])
        step = float(policy['cad_round_up_mm'])
        if not all(math.isfinite(v) for v in (shrink, reserve, step)) or not (0 < shrink < .1 and reserve > 0 and step > 0):
            raise ValueError('INVALID_ABS_DIMENSIONAL_ALLOWANCE')
        required = round(math.ceil((target+reserve)/(1-shrink)/step-1e-12)*step, 10)
    actual = required if measured_mm is None else float(measured_mm)
    lower = actual*(1-shrink)-reserve
    return dict(dimension=dimension, finished_minimum_mm=target,
        required_cad_minimum_mm=required, measured_cad_mm=actual,
        assumed_abs_shrink_fraction=shrink, total_dimension_error_reserve_mm=reserve,
        cad_round_up_mm=step, budgeted_finished_lower_bound_mm=lower,
        calibration_status='UNCALIBRATED_DESIGN_BUDGET' if policy else 'LEGACY_NO_ALLOWANCE',
        physical_finished_dimensions_validated=False, pdms_cure_shrinkage_characterized=False,
        sources=SOURCES if policy else [],
        passed=math.isfinite(actual) and actual >= required-1e-5 and lower >= target-1e-5)


def minimum_cad(cfg, dimension):
    return dimension_budget(cfg, dimension)['required_cad_minimum_mm']


def slicer_scale_qc(cfg, settings):
    """Reject a second compensation applied to already dimensioned geometry."""
    if not (cfg.get('abs_dimensional_allowance') or cfg.get('cad_dimensional_allowance')):
        return dict(applicable=False, passed=True)
    shrink = settings.get('filament_shrink')
    shrink = shrink if isinstance(shrink, list) else [shrink]
    try:
        scale_ok = bool(shrink) and all(float(str(v).rstrip('%')) == 100 for v in shrink)
        xy_ok = all(float(settings[k]) == 0 for k in ['xy_contour_compensation', 'xy_hole_compensation'])
        circle_ok = str(settings['enable_circle_compensation']) == '0'
    except (TypeError, ValueError, KeyError):
        scale_ok = xy_ok = circle_ok = False
    return dict(applicable=True, filament_shrink=shrink,
        xy_contour_compensation=settings.get('xy_contour_compensation'),
        xy_hole_compensation=settings.get('xy_hole_compensation'),
        enable_circle_compensation=settings.get('enable_circle_compensation'),
        passed=scale_ok and xy_ok and circle_ok)
