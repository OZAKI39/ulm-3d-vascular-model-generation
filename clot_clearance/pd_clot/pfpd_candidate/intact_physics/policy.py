"""Executable, purpose-scoped CONTROL-ORIENTED SUFFICIENCY V2 policy.

The 30% interval is a project screening assumption, never an experimental
confidence interval or an estimated biological error from these six solves.
"""
import copy
import math
import numpy as np

POLICY_ID = 'CONTROL_ORIENTED_SUFFICIENCY_V2'
EPSILON = 0.30
LEVELS = {
    0: 'LEVEL_0_REDUCED_ANALYTICAL_LOAD',
    1: 'LEVEL_1_FAST_FEM',
    2: 'LEVEL_2_FEM_FRACTURE',
    3: 'LEVEL_3_PD2_SPOT_CHECK',
    4: 'LEVEL_4_EXPLICIT_DYNAMIC_FRAGMENTATION',
    5: 'LEVEL_5_FLUID_FRAGMENT_TRANSPORT',
}
REQUIRED_METADATA = (
    'model_level', 'method_id', 'material_id', 'geometry_id', 'load_id',
    'quantity_definition', 'reference_scope', 'intended_use', 'accuracy_tier',
    'estimated_model_uncertainty', 'validation_status', 'ranking_status',
    'threshold_status', 'escalation_status',
)

def uncertainty(epsilon=EPSILON):
    if not math.isfinite(epsilon) or not 0 <= epsilon < 1:
        raise ValueError('Uncertainty must be finite in [0,1)')
    return dict(relative=epsilon, source=POLICY_ID,
                interpretation='PROJECT_SCREENING_POLICY_NOT_EXPERIMENTAL_UNCERTAINTY',
                empirically_estimated=False, biological_validation='NOT_ESTABLISHED')

def compare_reference(q, reference, floor, direction=None, reference_direction=None,
                      epsilon=EPSILON, reference_is_credible=False):
    """Assess a particular continuous output, not an entire model."""
    uncertainty(epsilon)
    if not all(math.isfinite(x) for x in [q, reference, floor]) or floor <= 0:
        raise ValueError('Finite values and a positive quantity-specific floor required')
    relative = abs(q-reference)/max(abs(reference), floor)
    magnitude = max(abs(q), floor)/max(abs(reference), floor)
    signs_ok = q*reference >= 0
    angle = None
    direction_status = 'NOT_APPLICABLE'
    if direction is not None and reference_direction is not None:
        a, b = np.asarray(direction, float), np.asarray(reference_direction, float)
        if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError('Matching finite direction vectors required')
        if np.linalg.norm(a) <= floor or np.linalg.norm(b) <= floor:
            direction_status = 'UNDEFINED_NEAR_ZERO'
        else:
            cosine = float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))
            angle = math.degrees(math.acos(min(1., max(-1., cosine))))
            direction_status = 'PASS' if angle <= 30. else 'SPOT_CHECK_REQUIRED'
    numerical_pass = (relative <= epsilon and 1/3 <= magnitude <= 3 and signs_ok
                      and direction_status not in ['SPOT_CHECK_REQUIRED', 'UNDEFINED_NEAR_ZERO'])
    status = ('PASS_FOR_SCREENING' if numerical_pass else 'SPOT_CHECK_REQUIRED') if reference_is_credible else 'EXPLORATORY_ONLY'
    return dict(relative_difference=relative, magnitude_ratio=magnitude,
                sign_gate=signs_ok, direction_angle_deg=angle, direction_status=direction_status,
                magnitude_gate=1/3 <= magnitude <= 3, relative_gate=relative <= epsilon,
                validation_status=status, credible_same_quantity_reference=reference_is_credible)

def ranking(a, b, epsilon=EPSILON):
    uncertainty(epsilon)
    if not all(math.isfinite(x) and x >= 0 for x in [a, b]):
        raise ValueError('Ranking policy requires finite nonnegative metrics')
    if a*(1-epsilon) > b*(1+epsilon):
        relation, status = 'A_GREATER_THAN_B', 'PASS_FOR_RANKING'
    elif b*(1-epsilon) > a*(1+epsilon):
        relation, status = 'B_GREATER_THAN_A', 'PASS_FOR_RANKING'
    else:
        relation, status = 'UNRESOLVED', 'RANKING_UNCERTAIN'
    return dict(status=status, relation=relation, epsilon=epsilon,
                interval_A=[a*(1-epsilon), a*(1+epsilon)],
                interval_B=[b*(1-epsilon), b*(1+epsilon)],
                required_positive_ratio=(1+epsilon)/(1-epsilon))

def risk_classification(q, threshold=None, epsilon=EPSILON):
    uncertainty(epsilon)
    if not math.isfinite(q) or q < 0:
        raise ValueError('Finite nonnegative risk source quantity required')
    if threshold is None:
        return dict(status='RISK_PROXY_RELATIVE_ONLY', ratio=None,
                    ABSOLUTE_DAMAGE_THRESHOLD='NOT_ESTABLISHED')
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError('Positive threshold with compatible evidence required')
    status = ('ROBUST_HIGH_RISK' if q*(1-epsilon) > threshold else
              'ROBUST_LOW_RISK' if q*(1+epsilon) < threshold else
              'NEAR_THRESHOLD_UNCERTAIN')
    return dict(status=status, ratio=q/threshold, threshold=threshold,
                relative_margin=(q-threshold)/threshold, epsilon=epsilon)

def escalation(*, threshold_status='RISK_PROXY_RELATIVE_ONLY', direction_anomaly=False,
               magnitude_anomaly=False, ranking_conflict=False, nonlocal_central=False,
               high_value_surrogate_uncertain=False, requires_fracture_output=False,
               impact=False, t_load=None, t_clot=None, fragment_transport=False,
               release_state_available=False):
    if t_load is not None and (not math.isfinite(t_load) or t_load <= 0):
        raise ValueError('t_load must be positive')
    if t_clot is not None and (not math.isfinite(t_clot) or t_clot <= 0):
        raise ValueError('t_clot must be positive')
    if impact or (t_load is not None and t_clot is not None and t_load <= t_clot):
        return 'DYNAMIC_MODEL_REQUIRED'
    if fragment_transport and release_state_available:
        return 'LEVEL_5_TRANSPORT_REQUIRED'
    if any([direction_anomaly, magnitude_anomaly, ranking_conflict,
            nonlocal_central, high_value_surrogate_uncertain]):
        return 'LEVEL_3_PD2_SPOT_CHECK'
    if requires_fracture_output or threshold_status in ['ROBUST_HIGH_RISK', 'NEAR_THRESHOLD_UNCERTAIN']:
        return 'LEVEL_2_FRACTURE_CANDIDATE'
    if threshold_status == 'RISK_PROXY_RELATIVE_ONLY':
        return 'NEEDS_EXPERIMENTAL_THRESHOLD'
    return 'NO_ESCALATION'

def decorate_result(result, request):
    """Attach current use context without approving reuse or old evidence."""
    r = request if isinstance(request, dict) else {}
    z = result
    tier = r.get('accuracy_tier', 'TIER_A')
    z.update(model_level=LEVELS[3] if 'PD2' in str(z.get('method_id')) else LEVELS[1],
             intended_use=r.get('intended_use'), accuracy_tier=tier,
             source=r.get('source'), use_requirement=copy.deepcopy(r.get('allowed_error')),
             estimated_model_uncertainty=uncertainty(),
             policy=dict(id=POLICY_ID, accuracy_tier=tier, intended_use=r.get('intended_use'),
                         allowed_error=copy.deepcopy(r.get('allowed_error')), source=r.get('source'),
                         default_solver='FEM', PD2_role='SPOT_CHECK_REFERENCE',
                         biological_accuracy='NOT_ESTABLISHED'),
             validation_status='EXPLORATORY_ONLY' if z.get('success') else 'NUMERICAL_OR_REQUEST_FAILURE',
             ranking_status='NOT_ASSESSED', threshold_status='RISK_PROXY_RELATIVE_ONLY',
             escalation_status='NEEDS_EXPERIMENTAL_THRESHOLD' if z.get('success') else 'SPOT_CHECK_REQUIRED')
    for k in ['geometry_id','load_id','quantity_definition','reference_scope']:
        z.setdefault(k, r.get(k, 'NOT_ASSESSED'))
    return z

def validate_metadata(result):
    missing = [k for k in REQUIRED_METADATA if k not in result]
    if missing:
        raise ValueError('Missing policy metadata: '+', '.join(missing))
    if result['model_level'] not in LEVELS.values():
        raise ValueError('Unknown model level')
    if result['accuracy_tier'] not in ['TIER_A','TIER_B','TIER_C']:
        raise ValueError('Unknown accuracy tier')
    u = result['estimated_model_uncertainty']
    if u['source'] != POLICY_ID or u['relative'] != EPSILON or u['empirically_estimated']:
        raise ValueError('Default screening policy must not be relabelled as measured error')
    if not result.get('success') and result.get('response') is not None:
        raise ValueError('Failed results cannot expose dataset labels')
    return True
