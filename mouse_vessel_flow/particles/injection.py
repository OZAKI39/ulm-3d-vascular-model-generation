"""Fixed source-order cohorts and original single-shot finite-size admission."""
import json
from ..paths import case_config,project_path
from .continuous_infusion import ContinuousInfusionSource,PopulationLedger
from .injection_method_c import TruncatedSonoVue
from .injection_admission import FiniteSizeAdmission
from .population_inlet_p9a4 import generate

def saved_events(name='roi_only'):
    return json.loads(project_path(case_config(name)['particles']['cohort']).read_text())['events']

def source_for_case(name='roi_only', *, environment=None):
    """Original source law and seed, with a separate relocated-code identity.

    Creating this object generates no events. A caller explicitly requests
    proposals or invokes generate; the original frozen ledger is never modified.
    """
    from ..flow_field import load_case
    from ..provenance import verify_source_integrity,sha256
    from .particle6_stepper import bind_query_dependency
    cfg=case_config(name);spec=cfg['particles'];env=environment or load_case(name)
    if env.new_flow_sha256!=cfg['flow']['sha256']:raise ValueError('Source/field case mismatch')
    cohort=json.loads(project_path(spec['cohort']).read_text())
    old=cohort.get('identity',cohort.get('source_identity'))
    class CaseSource(ContinuousInfusionSource):
        __init__=bind_query_dependency(ContinuousInfusionSource.__init__,{'NEW_FLOW_SHA256':cfg['flow']['sha256']})
    return CaseSource(sampler=env.sampler,distribution=env.distribution,checker=env.checker,
        concentration_m3=old['concentration_m3'],flow_sha256=cfg['flow']['sha256'],
        master_seed=old['master_seed'],source_contract_sha256=sha256(project_path(spec['contract'])),
        scientific_identity={'schema':'RELOCATED_SOURCE_V1','runtime_sources':verify_source_integrity(),
                             'parent_cohort_sha256':sha256(project_path(spec['cohort']))})
