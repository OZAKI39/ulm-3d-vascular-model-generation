"""P9-A.2 size-first source and finite-sphere conditional inlet sampling.

This module receives inlet velocity and solid geometry only. Trajectory results
are neither parameters nor inputs. Legacy sampling and integration are untouched.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from .sonovue_adapter import read_sonovue
from .injection_admission import FiniteSizeAdmission
from .injection_population import FluxClock, LinearProfile, ConstantMBConcentrationV0, C_MB
from .inlet_size_capacity import InletClearanceTree

METHOD = 'METHOD_C_SIZE_FIRST_FLUX_CONDITIONAL_V1'
LEGACY_METHOD = 'LEGACY_METHOD_B_POSITION_ANCHORED_SIZE_CONDITIONING'
DIAMETER_ROLE, POSITION_ROLE, ORIENTATION_ROLE = 920, 921, 922
DEFAULT_SEED = 2026092492


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def counter_rng(seed, birth_id, role, attempt):
    return np.random.default_rng([int(seed), int(birth_id), int(role), int(attempt)])


class TruncatedSonoVue:
    def __init__(self, root):
        self.original_contract, self.original, _ = read_sonovue(root)
        c, d = self.original_contract, self.original
        if c['within_bin_assumption'] != 'uniform' or c['pdf_type'] != 'piecewise constant':
            raise ValueError('UNRESOLVED_SOURCE_WITHIN_BIN_DEFINITION')
        self.mass = float(d.cdf(4.))
        width = np.maximum(0., np.minimum(d.high, 4.)-d.low)
        retained = d.probability*width/d.width
        self.contract = dict(contract_name='SONOVUE_DIAMETER_MAX4UM_V1',
            method=METHOD, D_max_m=4e-6, source_contract_sha256=sha(Path(root)/'contracts/SONOVUE_SAMPLER_CONTRACT_V0.json'),
            original_sonovue_source_sha256=c['sampler_source_sha256'], original_histogram_sha256=c['histogram_sha256'],
            original_distribution_mass_above_4um=1-self.mass, retained_mass=self.mass,
            rule='CONDITIONAL_DISTRIBUTION_D_LE_4UM; INVERSE_ORIGINAL_CDF(U*F_ORIGINAL(4UM))',
            normalization='DIVIDE_RETAINED_BIN_MASS_BY_F_ORIGINAL_4UM',
            within_bin_definition='UNIFORM_ALREADY_DEFINED_BY_FROZEN_SOURCE_CONTRACT',
            clipping=False, source_files_modified=False,
            generation_code_sha256=sha(__file__),
            bins=[dict(original_bin=int(i),low_um=float(d.low[i]),high_um=float(min(d.high[i],4.)),
                       original_probability=float(d.probability[i]), retained_probability=float(retained[i]/self.mass))
                  for i in range(len(d.low)) if width[i]>0])
        self.contract_sha256 = hashlib.sha256(canonical_bytes(self.contract)).hexdigest()

    def sample(self, rng):
        u = float(rng.random())
        diameter = float(self.original.inverse_cdf(u*self.mass))*1e-6
        assert 0 < diameter <= 4e-6
        return diameter, dict(uniform=u, original_cdf_probability=u*self.mass)

    def cdf(self, diameter_m):
        return np.minimum(self.original.cdf(np.asarray(diameter_m)*1e6)/self.mass, 1.)


def sample_truncated_sonovue_diameter(distribution, seed, particle_id, attempt):
    d, record = distribution.sample(counter_rng(seed, particle_id, DIAMETER_ROLE, attempt))
    return d, dict(record, diameter_draw_id=[int(seed), int(particle_id), DIAMETER_ROLE, int(attempt)], diameter_m=d)


class PositionSamplingProgressFailure(RuntimeError):
    def __init__(self, record):
        self.record = record
        super().__init__('POSITION_SAMPLING_PROGRESS_FAILURE')


def sample_flux_weighted_feasible_position(event, *, sampler, mapping, checker, seed, guard, bounds):
    """A rejected center cannot call, or otherwise advance, a diameter stream."""
    fixed = np.float64(event['diameter_m']).tobytes()
    attempts = []
    for attempt in range(guard):
        key = [int(seed), int(event['particle_id']), POSITION_ROLE, attempt]
        xyz, triangle = sampler.sample(counter_rng(*key))
        position = xyz[0]
        particle, status, detail = checker.check(event, position, {})
        assert np.float64(event['diameter_m']).tobytes() == fixed
        attempts.append(dict(position_draw_id=key, position_m=position.tolist(),
            inlet_triangle_id=int(mapping[triangle[0]]), admission_status=status,
            diameter_draw_id=event['diameter_draw_id'], diameter_float64_hex=fixed.hex(), **detail))
        if particle is not None:
            return position, int(mapping[triangle[0]]), attempts
    raise PositionSamplingProgressFailure(dict(particle_id=event['particle_id'], diameter_m=event['diameter_m'],
        diameter_draw_id=event['diameter_draw_id'], diameter_fixed_during_position_sampling=True,
        position_attempts=attempts, bounds=bounds, size_resampled_after_position_failure=False))


class MethodCSource:
    def __init__(self, *, sampler, wall, field, distribution, flow_sha256, geometry_sha256, seed=DEFAULT_SEED):
        self.sampler, self.wall, self.field = sampler, wall, field
        self.distribution, self.seed = distribution, int(seed)
        self.flow_sha256, self.geometry_sha256 = flow_sha256, geometry_sha256
        self.tree = InletClearanceTree(sampler, wall)
        self.capacity = self.tree.capacity()
        if not self.capacity['positive_flux_cap_covers_entire_cap']:
            raise ValueError('CAPACITY_OVER_WHOLE_CAP_REQUIRES_SEPARATE_ZERO_OR_NEGATIVE_FLUX_REGION_AUDIT')
        self.checker = FiniteSizeAdmission(None, None, wall=wall, field=field)
        self.clock = FluxClock(LinearProfile([0],[sampler.Q_m3_s]), ConstantMBConcentrationV0())

    def size_is_feasible(self, radius):
        for _ in range(8):
            lo, hi = self.capacity['radius_max_bracket_m']
            if radius <= lo:
                return True
            if radius > hi:
                return False
            self.capacity = self.tree.capacity(tolerance_m=self.capacity['tolerance_m']/8)
        raise RuntimeError('GLOBAL_SIZE_CAPACITY_BOUNDARY_UNRESOLVED')

    def event(self, particle_id, max_source_draws=100000):
        draws = []
        for index in range(max_source_draws):
            diameter, draw = sample_truncated_sonovue_diameter(self.distribution, self.seed, particle_id, index)
            feasible = self.size_is_feasible(diameter/2)
            draw['status'] = 'SIZE_FROZEN_FOR_POSITION_SAMPLING' if feasible else 'NO_FEASIBLE_INLET_POSITION_FOR_SIZE'
            draw['capacity_bracket_m'] = self.capacity['radius_max_bracket_m']
            draws.append(draw)
            if feasible:
                break
        else:
            raise RuntimeError('SOURCE_SIZE_PROGRESS_FAILURE')
        q = counter_rng(self.seed, particle_id, ORIENTATION_ROLE, 0).standard_normal(4)
        q /= np.linalg.norm(q)
        event = dict(particle_id=int(particle_id), species='MB', birth_time_s=self.clock.time_at(particle_id),
            scheduled_time_s=self.clock.time_at(particle_id), diameter_m=diameter, diameter_um=diameter*1e6,
            radius_m=diameter/2, volume_m3=float(4*np.pi*(diameter/2)**3/3), q=q.tolist(),
            diameter_draw_id=draw['diameter_draw_id'], diameter_source_draw=draw,
            diameter_global_rejections=len(draws)-1, source_diameter_draws=draws,
            diameter_fixed_during_position_sampling=True, seed=self.seed, admission_strategy=METHOD,
            source_distribution_contract_sha256=self.distribution.contract_sha256,
            flow_sha256=self.flow_sha256, geometry_sha256=self.geometry_sha256,
            entering_distribution='4UM_TRUNCATED_SONOVUE_CONDITIONED_ON_CURRENT_RIGID_SPHERE_INLET_PASSABILITY',
            birth_rate_semantics='ACTUAL_ENTERING_EVENT_RATE_C_MB_TIMES_Q', C_MB_m3=C_MB,
            orientation_role='ISOTROPIC_RANDOM_ORIENTATION_V0_MODEL_ASSUMPTION_NOT_MEASURED')
        proposal, mapping, bounds = self.tree.feasible_proposal(diameter/2)
        position, triangle, attempts = sample_flux_weighted_feasible_position(event,
            sampler=proposal, mapping=mapping, checker=self.checker, seed=self.seed,
            guard=bounds['position_guard'], bounds=bounds)
        clearance = self.wall.nearest_center_triangle(position)[1]-diameter/2
        xyz = self.sampler.triangles[triangle]
        bary12 = np.linalg.lstsq((xyz[1:]-xyz[0]).T, position-xyz[0], rcond=None)[0]
        bary = np.r_[1-bary12.sum(), bary12]
        event.update(position_xyz=position.tolist(), birth_center_m=position.tolist(), anchor_m=position.tolist(),
            anchor_triangle=triangle, inlet_triangle_id=triangle, position_proposal_count=len(attempts),
            position_attempts=attempts, flux_bounds=bounds, clearance_m=float(clearance), admission_status='ACCEPTED',
            local_flux_weight_m_s=float(max(0.,bary@self.sampler.q[triangle])),
            original_triangle_flux_m3_s=float(self.sampler.weights[triangle]),
            entry_transition_time_s=0., source_method=METHOD)
        return event


def generate_method_c_birth_event(source, particle_id):
    return source.event(particle_id)
