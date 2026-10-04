"""Read-only coupling diagnostics. No production solver or phase-law repair.

Phi reproduces the existing solver's moving inner lower bound. The separate
monolithic objective uses the PREVIOUS ACCEPTED LOAD bound. These two bounds
are deliberately not conflated.
"""
from dataclasses import dataclass
from pathlib import Path
import copy
import json
import numpy as np
from .adr_adapter_at1 import fixture, boundary_view, pair_weights, equilibrate_tensile
from .at1_evolution import critical_driver, at1_constrained_update, HighDamageGate
from .baqp import evaluate


ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT/'results/pfpd_at1_evolution_v1'
OUT = ROOT/'verification/pfpd_at1_coupling_stability_audit_v1'


def read_json(path): return json.loads(Path(path).read_text())


def save_json(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def audit_config(): return read_json(OUT/'provenance/PREDECLARED_AUDIT_CONFIG.json')


def physics_config(): return read_json(ROOT/audit_config()['physics_config_source'])


def read_npz(path):
    with np.load(path,allow_pickle=False) as data:
        return {k:data[k].copy() for k in data.files}


def green_minimum(state):
    green = .5*(state.Ftilde.transpose(0,2,1)@state.Ftilde-np.eye(3))
    return float(np.linalg.eigvalsh(green).min())


def projected_kkt(s, q_scaled, lower, upper):
    """Natural projected residual, including both bounds, with no active-set fudge."""
    return s-np.clip(s-q_scaled,lower,upper)


class AuditProblem:
    def __init__(self, config, strain, previous_load_phase):
        self.config = copy.deepcopy(config)
        self.cloud,self.family,self.material,self.top = fixture(config)
        self.strain = float(strain)
        self.boundary = boundary_view(self.cloud,self.top,strain).X
        self.free = ~self.cloud.fixed
        self.nx = 3*int(self.free.sum()); self.ns = len(self.family.pairs)
        self.lower = np.array(previous_load_phase,float,copy=True)
        self.sc = float(config['phase']['sc'])
        if self.lower.shape != (self.ns,) or not np.isfinite(self.lower).all() or np.any((self.lower<0)|(self.lower>=self.sc)):
            raise ValueError('Valid previous-load lower bound required')
        self.lower.setflags(write=False)
        self.weights = pair_weights(self.family)
        self.yc = critical_driver(config['phase']['Gc_demo_J_per_m2'],self.family.horizon,config['phase']['kernel'])
        self.length_scale = self.cloud.spacing*config['phase']['threshold_engineering_strain']
        self.energy_scale = float(self.family.volume.sum()*self.yc)

    def prescribe(self,x):
        x = np.array(x,float,copy=True); x[self.cloud.fixed] = self.boundary[self.cloud.fixed]
        return x

    def pack(self,x,s):
        return np.r_[((np.asarray(x)[self.free]-self.family.X[self.free])/self.length_scale).ravel(),np.asarray(s)]

    def unpack(self,z):
        z = np.asarray(z,float)
        if z.shape != (self.nx+self.ns,) or not np.isfinite(z).all(): raise ValueError('Invalid joint state')
        x = self.boundary.copy()
        x[self.free] = self.family.X[self.free]+self.length_scale*z[:self.nx].reshape((-1,3))
        return x,z[self.nx:].copy()

    def state(self,x,s):
        return evaluate(self.family,x,self.material,g=(1-np.asarray(s))**2,h=np.ones(self.ns))

    def parts(self,state,s):
        crack = float(self.weights@(2*self.yc*s))
        return dict(elastic_J=state.energy,crack_J=crack,total_J=state.energy+crack)

    def objective_gradient(self,z):
        """Analytic gradient of identical BAQP pair-weighted candidate energy.

        Domain checks are separate: derivative probes can sample the algebraic
        energy extension, but no optimizer trial there is accepted as physical.
        """
        x,s = self.unpack(z); state = self.state(x,s)
        q = -2*(1-s)*state.W0+2*self.yc
        gradient = np.r_[-state.force[self.free].ravel()*self.length_scale,self.weights*q]/self.energy_scale
        parts = self.parts(state,s)
        if not np.isfinite(gradient).all() or not np.isfinite(parts['total_J']): raise FloatingPointError('Nonfinite objective')
        return parts['total_J']/self.energy_scale,gradient

    def metrics(self,x,s):
        state = self.state(x,s); q = -2*(1-s)*state.W0+2*self.yc
        scaled_q = q/np.maximum(state.W0,self.yc)
        kkt = projected_kkt(s,scaled_q,self.lower,self.sc)
        return dict(**self.parts(state,s),force_residual=float(np.linalg.norm(state.force[self.free])/
            self.config['mechanical_equilibrium']['force_floor_N']),projected_phase_KKT=float(np.max(abs(kkt))),
            maximum_s=float(np.max(s)),mean_s=float(np.average(s,weights=self.weights)),
            minimum_Green_strain=green_minimum(state),minimum_lower_margin=float(np.min(s-self.lower)),
            minimum_upper_margin=float(np.min(self.sc-s)))

    def classify(self,x,s,settings):
        state = self.state(x,s); q = (-2*(1-s)*state.W0+2*self.yc)/np.maximum(state.W0,self.yc)
        near_lower = s-self.lower<=settings['active_s_tolerance']
        strong_lower = near_lower & (q>settings['active_q_scaled_tolerance'])
        interior = (s-self.lower>settings['active_s_tolerance']) & (self.sc-s>settings['near_sc_tolerance'])
        return dict(LOWER_BOUND_ACTIVE=strong_lower,INTERIOR_ACTIVE=interior,
            WEAK_LOWER_OR_UNCLASSIFIED=~(strong_lower|interior),
            NEAR_THRESHOLD=abs(state.W0/self.yc-1)<=settings['near_threshold_relative_tolerance'],
            NEAR_SC=self.sc-s<=settings['near_sc_tolerance'],q_scaled=q)

    def phi(self,s_input,x_start):
        """Exactly the existing one-step mechanics -> irreversible phase map."""
        s_input = np.array(s_input,float,copy=True)
        if s_input.shape != (self.ns,) or np.any(s_input<self.lower) or np.any(s_input>=self.sc):
            raise ValueError('Phi input outside physical lower/sc bounds')
        result = equilibrate_tensile(self.cloud,self.family,self.material,self.top,
                                    self.prescribe(x_start),s_input,self.strain,self.config)
        x = result.x.copy(); state = self.state(x,s_input)
        metrics = self.metrics(x,s_input)
        answer = dict(x_equilibrium=x,Y=state.W0.copy(),mechanical_summary=result.summary,
                      ADR_iterations=result.summary['iterations'],before=metrics,s_input=s_input.copy())
        if not result.converged:
            answer['status']='MECHANICAL_EQUILIBRIUM_FAILED';return answer
        if metrics['minimum_Green_strain']<self.config['tolerances']['tensile_min_green_strain']:
            answer['status']='PURE_TENSION_DOMAIN_VIOLATION';return answer
        try: update = at1_constrained_update(s_input,state.W0,self.yc,self.sc)
        except HighDamageGate as exc:
            answer.update(status='HIGH_DAMAGE_GATE_REACHED',unaccepted_proposal=exc.update.s_new.copy());return answer
        answer.update(status='OK',s_output=update.s_new.copy(),s_trial=update.s_trial.copy(),
            selector_gap=update.s_trial-s_input,active_set=update.active_set_status.copy(),
            q=update.derivative_at_solution.copy(),after=self.metrics(x,update.s_new))
        return answer


def load_saved_iteration(run,step,iteration):
    config = read_json(OLD/run/'CONFIG.json'); folder = OLD/run/f'load_{step:03d}'
    row = read_json(folder/f'iteration_{iteration:03d}_phase.json')
    phase = read_npz(folder/f'iteration_{iteration:03d}_phase.npz')
    previous = read_npz(OLD/run/f'load_{step-1:03d}/accepted_state.npz') if step else dict(x=fixture(config)[1].X,s_pf_bond=np.zeros(len(phase['s_old'])))
    initial = (read_npz(folder/f'iteration_{iteration-1:03d}_mechanical.npz')['x']
               if iteration else previous['x'])
    problem = AuditProblem(config,row['strain'],previous['s_pf_bond'])
    return problem,dict(x=phase['x'],s=phase['s_old'],s_next=phase['s_new'],x_start=problem.prescribe(initial),
                        previous=previous,row=row,path=str(folder/f'iteration_{iteration:03d}_phase.npz'))


def load_audit_state(key):
    selection = audit_config()['states'][key]
    if selection['kind']=='iteration_input':
        return load_saved_iteration(selection['run'],selection['load_step'],selection['iteration'])
    folder = OLD/selection['run']/f"load_{selection['load_step']:03d}"
    data = read_npz(folder/'accepted_state.npz')
    prev = read_npz(OLD/selection['run']/f"load_{selection['load_step']-1:03d}/accepted_state.npz")
    config = read_json(OLD/selection['run']/'CONFIG.json')
    problem = AuditProblem(config,selection['strain'],prev['s_pf_bond'])
    return problem,dict(x=data['x'],s=data['s_pf_bond'],x_start=data['x'].copy(),previous=prev,
                        path=str(folder/'accepted_state.npz'))


def deterministic_starts():
    """All starts at target load, with the same physical previous-load bound."""
    settings = audit_config(); prev = read_npz(OLD/'loading_unloading/load_009/accepted_state.npz')
    problem = AuditProblem(physics_config(),settings['subdivision_endpoint'],prev['s_pf_bond'])
    lifted = prev['x'].copy()
    increment = (settings['subdivision_endpoint']-settings['subdivision_start'])*np.ptp(problem.family.X[:,2])
    lifted[:,2] += increment*(problem.family.X[:,2]-problem.family.X[:,2].min())/np.ptp(problem.family.X[:,2])
    lifted = problem.prescribe(lifted)
    first = read_npz(OLD/'loading_unloading/load_010/iteration_000_phase.npz')
    later = read_npz(OLD/'loading_unloading/load_010/iteration_023_phase.npz')
    rng = np.random.default_rng(settings['monolithic_perturbation_seed'])
    perturb = prev['s_pf_bond']+settings['monolithic_perturbation_phase_amplitude']*rng.uniform(0,1,problem.ns)
    pairs = [(lifted,prev['s_pf_bond']),(first['x'],first['s_new']),(later['x'],later['s_new']),(lifted,perturb)]
    return problem,dict(zip(settings['monolithic_starts'],pairs))
