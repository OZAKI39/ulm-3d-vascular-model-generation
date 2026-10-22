"""Opt-in new run from a retained cohort; frozen campaigns are never resumed/overwritten."""
import json
from ..paths import case_config,new_output,project_path
from ..provenance import verify_source_integrity,sha256
from ..flow_field import load_case
from .injection import saved_events
from .dynamics import bind_job

def run_saved_cohort(name,output,particle_ids):
    """Explicit IDs only, serial CPU execution. New identity does not impersonate an old run."""
    cfg=case_config(name);spec=cfg['particles'];events=saved_events(name)
    ids=set(particle_ids)
    if not ids or len(ids)!=len(particle_ids):raise ValueError('Supply distinct explicit particle IDs')
    selected=[e for e in events if e['particle_id'] in ids]
    if len(selected)!=len(ids):raise ValueError('Unknown particle ID')
    source=verify_source_integrity();env=load_case(name)
    out=new_output(output)
    identity=dict(schema='MOUSE_RELOCATED_EXPLICIT_RUN_V1',flow_sha256=cfg['flow']['sha256'],dt_s=spec['dt_s'],
        source_cohort_sha256=sha256(project_path(spec['cohort'])),runtime_sources=source,
        parent_case=name,original_campaign_resumed=False)
    (out/'identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    job=bind_job(env,identity,spec['dt_s'],accelerated=spec['accelerated'])
    rows=[job((event,str(out/f"mb_{event['particle_id']:06d}"))) for event in selected]
    (out/'metrics.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows
