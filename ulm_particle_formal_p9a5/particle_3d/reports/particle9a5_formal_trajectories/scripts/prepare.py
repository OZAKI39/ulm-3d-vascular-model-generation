"""Freeze outcome-blind benchmark cohorts and audit prior formal time limits."""
from pathlib import Path
import sys,json,csv,hashlib,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import *
R=ROOT/REL;data=R/'data'

def main():
 population=load_births(ROOT)
 contract=ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json'
 assert digest(contract)==SOURCE_CONTRACT_SHA
 c=json.loads(contract.read_text())
 flow=ROOT/'particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/inputs/NEW.vtu'
 assert digest(flow)==FLOW_SHA
 cohort_sha=write_new(data/'CORE500_COHORT.json',cohort(population,500,c['master_seed']))
 write_new(data/'BENCHMARK24_COHORT.json',dict(selection='FIRST_24_OF_FROZEN_CORE500',core_cohort_sha256=cohort_sha,events=population['events'][:24]))
 hist=ROOT/'particle_3d/outputs/particle8_2a_ppt';rows=list(csv.DictReader((hist/'data/trajectory_catalog.csv').open()))
 limits=[r for r in rows if r['end_reason']=='PHYSICAL_RESIDENCE_HORIZON_REACHED']
 ages=sorted(set(float(r['last_age_s']) for r in limits));assert ages==[1.5]
 docs=[hist/'data/trajectory_catalog.csv',hist/'PPT_REVIEW.md',hist/'PPT_VALIDATION.json',ROOT/'particle_3d/reports/particle8_2a_ppt/PPT_REVIEW.md',ROOT/'particle_3d/reports/particle8_2a_ppt/PPT_VALIDATION.json']
 validation=json.loads((hist/'PPT_VALIDATION.json').read_text());assert validation['INTEGRATED_TRAJECTORIES']==len(rows)==500
 assert validation['summary']['end_reasons']['PHYSICAL_RESIDENCE_HORIZON_REACHED']==len(limits)
 write_new(data/'historical_horizon_audit.json',dict(previous_formal_horizon_s=1.5,formal_rows=len(rows),horizon_reached_count=len(limits),horizon_reached_ids=[int(r['event_id']) for r in limits],observed_horizon_ages_s=ages,document_note='PPT review/validation identify the 500-cohort; actual horizon is verified from all 21 censored catalog ages.',evidence={str(p.relative_to(ROOT)):digest(p) for p in docs}))
 baseline=json.loads((R/'logs/baseline/summary.json').read_text());assert baseline['tests']==165 and baseline['failures']==baseline['skipped']==0
 combined=ET.Element('testsuites')
 for path in sorted((R/'logs/baseline').rglob('*.xml')):
  tree=ET.parse(path).getroot()
  for suite in ([tree] if tree.tag=='testsuite' else list(tree)):combined.append(suite)
 ET.ElementTree(combined).write(R/'logs/baseline_tests.xml',encoding='utf-8',xml_declaration=True)
 policy=dict(schema='P9A5_FORMAL_PRODUCTION_V1',flow_sha256=FLOW_SHA,source_contract_sha256=SOURCE_CONTRACT_SHA,source_births_sha256=BIRTHS_SHA,core_N=500,core_cohort_sha256=cohort_sha,initial_horizon_s=3.,horizon_schedule_s=list(HORIZONS),dt_s=DT,dt_authorization='Explicit subsequent user instruction: dt = 1.0 ms; supersedes original 0.25 ms request',increment_before_2000=250,increment_after_2000=500,preferred_maximum_N=5000,population_selection='CONTIGUOUS_PREFIX_ONLY',stop_population='ALL_O1_O2_O3_NATURALLY_COMPLETED_OR_PREDECLARED_RESOURCE_MAXIMUM',point_role='POST_COHORT_DESCRIPTIVE_ONLY',stationary_role='EXISTING_CONTACT_RANK3_KKT_SUPPORT_REQUIRED',production_authorization='EXPLICIT_P9A5_USER_REQUEST; P9A4 source contract remains historical and unmodified',maximum_provider_call_budget='16000 per historical 1.5 s, scaled only with authorized horizon; all rejection/progress guards remain unchanged',GPU_TRAJECTORY_KERNEL_AVAILABLE=False)
 write_new(ROOT/'particle_3d/contracts/P9A5_FORMAL_PRODUCTION_V1.json',policy)
 write_new(data/'input_audit.json',dict(core_cohort_sha256=cohort_sha,source_events=len(population['events']),source_contract_sha256=digest(contract),actual_flow_sha256=digest(flow),baseline_tests=baseline,old_reports_modified=False,actual_previous_horizon_s=1.5))
 print(json.dumps(dict(core500_sha256=cohort_sha,flow_sha256=digest(flow),previous_formal_horizon_s=1.5,baseline_tests=baseline['tests']),indent=2))
if __name__=='__main__':main()
