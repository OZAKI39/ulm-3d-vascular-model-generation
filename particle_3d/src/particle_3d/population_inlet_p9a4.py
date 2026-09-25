"""Authoritative Network-H0 input adapter and Linux CPU proposal orchestration."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace
import hashlib
import importlib.util
import json
import multiprocessing
import numpy as np
import pyvista as pv
from .continuous_infusion import (ContinuousInfusionSource, PopulationLedger,
    NEW_FLOW_SHA256, MODEL, canonical_bytes, sha256)
from .injection_admission import FiniteSizeAdmission
from .injection_method_c import TruncatedSonoVue
from .inlet_flux import frozen_boundary_flux

MASTER_SEED = 2026092594
REPORT_RELATIVE = 'particle_3d/reports/particle9a4_population_inlet'
_WORKER_SOURCE = None


def load_new_environment(root):
    root = Path(root).resolve()
    reference = root/'particle_3d/reports/network_derived_flow_mb_validation_v1'
    authority = json.loads((reference/'data/old_new_flow_contract.json').read_text())
    if authority['sha256']['NEW'] != NEW_FLOW_SHA256:
        raise ValueError('Authoritative NEW metadata changed')
    bundle = reference/'server_bundle'
    if sha256(bundle/'inputs/NEW.vtu') != NEW_FLOW_SHA256:
        raise ValueError('NEW file SHA mismatch: OLD substitution prohibited')
    # Reuse the audited input adapter, including old/new canonical mesh checks.
    path = reference/'scripts/runner.py'
    spec = importlib.util.spec_from_file_location('p9a4_audited_network_input_adapter', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    env = module.make_environment(bundle, 'NEW')
    audit, samplers = frozen_boundary_flux(env.mesh, pv.read(bundle/'inputs/NEW.vtu'), env.boundaries)
    env.sampler, env.flux_audit = samplers['INLET'], audit
    env.audit = audit
    env.new_flow_path, env.new_flow_sha256 = bundle/'inputs/NEW.vtu', NEW_FLOW_SHA256
    env.repo_root = root
    env.distribution = TruncatedSonoVue(root/'sonovue_size_distribution_v0')
    env.checker = FiniteSizeAdmission(None, None, wall=env.wall, field=env.field)
    return env


def make_contract(root, env, *, lambda_enter_estimate=None):
    """Nominal full-population anchor; continuous infusion is an assumption."""
    root = Path(root)
    old_path = root/'particle_3d/contracts/PARTICLE7_INLET_POPULATION_V0.json'
    old = json.loads(old_path.read_text())
    if old['MB_NOMINAL_CONCENTRATION_V0'] != 'NOMINAL_WELL_MIXED_POST_BOLUS_ANCHOR':
        raise ValueError('Legacy concentration provenance changed')
    total = float(old['MB_NOMINAL_CONCENTRATION_M3'])
    modeled = total*env.distribution.mass
    original = env.distribution.original_contract
    return dict(contract_name='P9A4_CONTINUOUS_INFUSION_V1', model=MODEL,
        source_distribution='SONOVUE_D_LE_4UM_CONDITIONAL',
        D_min_m=float(env.distribution.original.low[0])*1e-6, D_max_m=4e-6,
        source_original_distribution_sha256=original['histogram_sha256'],
        source_original_sampler_sha256=original['sampler_source_sha256'],
        source_original_contract_sha256=sha256(root/'sonovue_size_distribution_v0/contracts/SONOVUE_SAMPLER_CONTRACT_V0.json'),
        F_original_4um=env.distribution.mass, within_bin='UNIFORM_FROM_ORIGINAL_FROZEN_CONTRACT',
        C_total_m3=total, C_modeled_D_le_4um_m3=modeled,
        concentration_semantics='NOMINAL_TOTAL_ANCHOR_TIMES_F4; STEADY_LEVEL_MODEL_ASSUMPTION',
        concentration_provenance=dict(legacy_contract_sha256=sha256(old_path),
            number_dose=20000000, body_mass_kg=.0326, blood_volume_ml_per_kg=72,
            blood_volume_ml=.0326*72, unrounded_nominal_total_m3=20000000/(.0326*72)*1e6,
            inherited_rounded_total_m3=total, measured=False,
            steady_infusion_level_is_model_assumption=True,
            limitation='POST_BOLUS_DOSE_ANCHOR_IS_NOT_MEASURED_CONTINUOUS_INFUSION_CONCENTRATION',
            references=old['sources'][1:]),
        flow_SHA=NEW_FLOW_SHA256, Q_in_m3_s=env.sampler.Q_m3_s,
        lambda_source_s_inv=modeled*env.sampler.Q_m3_s,
        lambda_enter=dict(formula='C_modeled * integral f_source(D) Q_acc(D) dD',
            estimate=lambda_enter_estimate, estimate_role='INDEPENDENT_INLET_DIAGNOSTIC_NOT_ACCEPTANCE_GATE'),
        arrival_model='HOMOGENEOUS_POISSON_STEADY_CONTINUOUS_INFUSION',
        position_model='FEM_POSITIVE_FLUX_WEIGHTED', finite_size_model='SINGLE_SHOT_WALL_CLEARANCE_THINNING',
        open_inlet='CENTER_ON_CAP; SPHERE_MAY_STRADDLE_UPSTREAM', dynamics='UNCHANGED_P9_A_1',
        active={}, rejected_event_policy='PERSIST_NO_RETRY_NO_DELAY_NO_PARTICLE_ID',
        master_seed=MASTER_SEED, production_500_authorized=False,
        flow_limitations=dict(model='P1_P1_VMS',root_section_max_residual_percent=3.454212,
            root_section_rms_residual_percent=1.789924,scope='KNOWN_ACCEPTED_LIMITATION_NO_CFD_CHANGES'))


def make_source(env, contract, master_seed=MASTER_SEED):
    if contract['flow_SHA'] != NEW_FLOW_SHA256 or contract['Q_in_m3_s'] != env.sampler.Q_m3_s:
        raise ValueError('P9A4 flow or flux contract mismatch')
    if contract['source_distribution'] != 'SONOVUE_D_LE_4UM_CONDITIONAL':
        raise ValueError('P9A4 requires the unique conditional source')
    if contract['F_original_4um'] != env.distribution.mass:
        raise ValueError('Source CDF differs from contract')
    if contract['C_modeled_D_le_4um_m3'] != contract['C_total_m3']*env.distribution.mass:
        raise ValueError('Conditional concentration must include the original retained mass')
    snapshot_path = env.repo_root/'particle_3d/reports/network_derived_flow_mb_validation_v1/data/code_snapshot_manifest.json'
    original = json.loads(snapshot_path.read_text())
    for relative, expected in original.items():
        if sha256(env.repo_root/relative) != expected:
            raise ValueError('Original protected science changed: '+relative)
    scientific_identity = dict(original_snapshot_sha256=sha256(snapshot_path),
        continuous_infusion_sha256=sha256(Path(__file__).with_name('continuous_infusion.py')),
        population_adapter_sha256=sha256(__file__))
    return ContinuousInfusionSource(sampler=env.sampler, distribution=env.distribution,
        checker=env.checker, concentration_m3=contract['C_modeled_D_le_4um_m3'],
        flow_sha256=NEW_FLOW_SHA256, master_seed=master_seed,
        source_contract_sha256=hashlib.sha256(canonical_bytes(contract)).hexdigest(),
        scientific_identity=scientific_identity)


def _proposal_shard(ids):
    return [_WORKER_SOURCE.proposal(sid) for sid in ids]


def generate(source, count, *, workers=1, ledger=None, chunk_size=128, progress=None):
    """Linux fork workers share immutable geometry; returned rows merge in ID order."""
    global _WORKER_SOURCE
    if type(count) is not int or count < 0 or type(workers) is not int or not 1 <= workers <= 6:
        raise ValueError('Nonnegative integer count and workers in [1,6] required')
    if type(chunk_size) is not int or chunk_size < 1:
        raise ValueError('Positive integer chunk size required')
    ledger = PopulationLedger(source.identity) if ledger is None else ledger
    if ledger.identity != source.identity:
        raise ValueError('Resume identity changed')
    start = ledger.next_source_event_id
    chunks = [list(range(a, min(a+chunk_size, start+count))) for a in range(start, start+count, chunk_size)]
    if source.rate == 0:
        if count:
            raise StopIteration('ZERO_SOURCE_RATE')
        return ledger
    _WORKER_SOURCE = source
    if workers == 1:
        for ids in chunks:
            ledger.merge(_proposal_shard(ids))
            if progress: progress(len(ledger.rows))
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('fork')) as pool:
            # map yields in source order even when execution completes out of order.
            for rows in pool.map(_proposal_shard, chunks):
                ledger.merge(rows)
                if progress: progress(len(ledger.rows))
    return ledger
