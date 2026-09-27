"""Performance-stage gates; all measurements use the frozen physical mesh."""
import json
import math
from pathlib import Path
import numpy as np
from sv_validation.provenance import sha256
from sv_validation.validation import require
from .transient_checks import interval_metrics

ROOT = Path(__file__).resolve().parents[3]
REPORT = ROOT/'reports/sv1_3'
CONFIG = ROOT/'configs/sv1_3'
OUTPUT = ROOT/'outputs/sv1_3'
LOG = ROOT/'logs/sv1_3'

def load(name, stage='sv1_3'):
    return json.loads((ROOT/'reports'/stage/(name+'.json')).read_text())

def frozen_policy():
    return json.loads((CONFIG/'policy.json').read_text())

def check_reference(files=None):
    files = files if files is not None else load('reference_freeze')['files']
    for item in files:
        path = ROOT/item['path']
        require(path.is_file() and sha256(path)==item['sha256'], 'REFERENCE_CHANGED: '+str(path))
    return True

class SteadyStopMonitor:
    def __init__(self, measure, policy):
        self.measure, self.policy = measure, policy
        self.states, self.intervals = [], []
        self.previous_u = None

    def observe(self, state, velocity, linear_failures=0, nonlinear_failures=0):
        if self.states:
            self.intervals.append(interval_metrics(self.measure, velocity, self.previous_u,
                state, self.states[-1], self.policy['save_interval_steps']))
        self.states.append(state)
        self.previous_u = velocity.copy()
        return stop_gate(self.intervals, state, linear_failures, nonlinear_failures, self.policy)

def stop_gate(intervals, state, linear_failures, nonlinear_failures, policy):
    n = policy['steady_last_intervals']
    consecutive = len(intervals)>=n and all(
        math.isfinite(r[k]) and r[k]<=policy[limit]
        for r in intervals[-n:]
        for k,limit in [('E_u','velocity_change_limit'),('E_Q','flow_change_limit')])
    mass = all(math.isfinite(state[k]) and state[k]<=policy['mass_limit']
               for k in ('epsilon_Q','epsilon_mass'))
    return bool(consecutive and mass and state['velocity_finite'] and state['pressure_finite']
                and state['wall_noslip_pass'] and not linear_failures and not nonlinear_failures)

def equivalence(measure, candidate, reference, policy):
    u,p = measure.read(candidate)
    ur,pr = measure.read(reference)
    m,mr = measure.measure(u,p),measure.measure(ur,pr)
    l2 = lambda a: measure.velocity_l2(a[:,None])
    pressure_scale = max(abs(v) for v in mr['pressure_range_pa'])
    e = {
        'velocity_relative_volume_L2': measure.velocity_l2(u-ur)/measure.velocity_l2(ur),
        'pressure_relative_volume_L2': l2(p-pr)/l2(pr),
        'relative_Qin': abs(m['Q_in_m3_s']-mr['Q_in_m3_s'])/abs(mr['Q_in_m3_s']),
        'relative_Qout': {r:abs(v-mr['outlet_flows_m3_s'][r])/abs(mr['outlet_flows_m3_s'][r]) for r,v in m['outlet_flows_m3_s'].items()},
        'absolute_fraction_difference': {r:abs(v-mr['outlet_fractions'][r]) for r,v in m['outlet_fractions'].items()},
        'port_pressure_scale_normalized': {r:abs(v-mr['area_average_pressure_pa'][r])/pressure_scale for r,v in m['area_average_pressure_pa'].items()},
        'port_pressure_relative': {r:abs(v-mr['area_average_pressure_pa'][r])/abs(mr['area_average_pressure_pa'][r]) for r,v in m['area_average_pressure_pa'].items()},
        'relative_max_velocity': abs(m['velocity_max_m_s']-mr['velocity_max_m_s'])/mr['velocity_max_m_s']}
    thresholds=policy['equivalence']
    checks={k:all(math.isfinite(x) and x<=thresholds[k] for x in (v.values() if isinstance(v,dict) else [v]))
            for k,v in e.items() if k in thresholds}
    checks['finite_and_wall']=m['velocity_finite'] and m['pressure_finite'] and m['wall_noslip_pass']
    checks['mass']=max(m['epsilon_Q'],m['epsilon_mass'])<=policy['mass_limit']
    return {'status':'PASS' if all(checks.values()) else 'FAIL','errors':e,'checks':checks,
            'candidate':m,'reference':mr,'candidate_path':str(candidate),'reference_path':str(reference),
            'candidate_sha256':sha256(candidate),'reference_sha256':sha256(reference),
            'pressure_normalization':'one frozen reference max(abs(p)) scale for all ports; raw pressure, no offset'}

def gpu_gate(data):
    if data.get('mat_type') not in ('seqaijcusparse','mpiaijcusparse','aijcusparse') or data.get('vec_type') not in ('seqcuda','mpicuda','cuda'):
        return 'GPU_BACKEND_NOT_REACHABLE'
    if data.get('oom') or data.get('peak_memory_fraction',1)>=.9:
        return 'DEVICE_MEMORY_MARGIN_TOO_SMALL'
    if not data.get('scientific_equivalence') or not data.get('mass_pass'):
        return 'GPU_REJECTED_NUMERICAL_DIFFERENCE'
    if not data.get('linear_convergence'):
        return 'LINEAR_SOLVER_FAILURE'
    if data.get('transfer_scales_with')=='KSP_iterations' or not data.get('matrix_resident') or not data.get('vectors_resident'):
        return 'EXCESSIVE_HOST_DEVICE_TRANSFER'
    speed=data.get('speedup',0)
    return 'PASS' if speed>=1.25 else ('GPU_WORKS_BUT_NOT_WORTH_COMPLEXITY' if speed>1 else 'GPU_SLOWER_THAN_CPU')

def repeated_timing(times):
    require(len(times)>=2 and all(math.isfinite(t) and t>0 for t in times),'At least two measured runs required')
    third=abs(times[1]-times[0])/times[0]>.10
    require(not third or len(times)>=3,'Third run required by >10% difference')
    return float(np.median(times[:3])) if third else times[1]

def select_production(cpu, gpu=None):
    require(cpu.get('science_equivalent') and cpu.get('mass_pass') and cpu.get('steady_pass')
        and cpu.get('checkpoint_complete') and cpu.get('linear_failures')==0
        and cpu.get('nonlinear_failures')==0,'EARLY_STOP_INVALID')
    if gpu and gpu_gate(gpu)=='PASS' and gpu.get('full_run_pass') and gpu.get('full_run_speedup',0)>=1.25:
        return gpu['name']+'_PRODUCTION'
    return 'CPU_EARLY_STOP_PRODUCTION'
