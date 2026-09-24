#!/usr/bin/env python3
"""Final evidence gates: completion is a scientific audit status, never production PASS."""
from pathlib import Path
import subprocess,xml.etree.ElementTree as ET
import numpy as np
from particle_3d.routing_stationary_audit import *
from particle_3d.particle8_replay import REPO
R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';D=R/'data'

def main():
 s=read(D/'audit_summary.json');raw=read(D/'point_completed.json')['results'];p65=read(D/'p65_completed.json')['results'];st=read(D/'stationary_audit.json');pair=read(D/'paired_routing.json');manifest=read(D/'figure_manifest.json')
 suites=ET.parse(R/'logs/audit_tests.xml').getroot();suites=list(suites.iter('testsuite'));tests={k:sum(int(e.get(k,0)) for e in suites) for k in ['tests','errors','failures','skipped']}
 protected=read(D/'protected_before.json');changed=[p for p,h in protected.items() if not (REPO/p).is_file() or sha(REPO/p)!=h]
 before=read(R/'logs/git_before.json');sameindex=subprocess.check_output(['git','write-tree'],cwd=REPO,text=True)==before['git write-tree'];protection=dict(file_count=len(protected),mismatches=changed,index_unchanged=sameindex);dump(D/'protected_verification.json',protection)
 for name,cmd in [('git_diff_stat.txt',['git','diff','--stat']),('git_status_final.txt',['git','status','--short'])]:(R/'logs'/name).write_text(subprocess.check_output(cmd,cwd=REPO,text=True))
 new=[str(p.relative_to(REPO)) for parent in ['particle_3d/src/particle_3d','particle_3d/scripts','particle_3d/tests/particle9a1_audit'] for p in (REPO/parent).glob('*.py') if str(p.relative_to(REPO)) not in protected]
 (R/'logs/new_audit_files.txt').write_text('\n'.join(sorted(new))+'\n');dump(D/'audit_source_manifest.json',{p:sha(REPO/p) for p in new})
 figures_ok=all((R/'figures'/(f['figure']+'.png')).is_file() and sha(R/'figures'/(f['figure']+'.png'))==f['png_sha256'] and all((D/x).is_file() for x in f['sources']) for f in manifest)
 remote=read(D/'remote_protection_verification.json')
 executed_radius_sha=sha(R/'reference/computation_stationary_radius_audit_v1.py')
 resolution=read(D/'point_resolution.json')['results'];controls=[r for r in resolution if r['baseline_outlet']!='NO_EXIT'];resolution_ok=all(all(t['point_outlet']==r['baseline_outlet'] for t in r['trials']) for r in controls)
 gates=dict(all_original_candidates_computed=len(raw)==s['candidate_count']==read(D/'event_identity.json')['candidate_count'],accepted500_computed=sum(r['accepted'] for r in raw)==500,
  P65_same500_computed=len(p65)==500,P9_saved500_aligned=len(pair)==500 and len({r['particle_id'] for r in pair})==500,
  outlet_bias_localized=bool(s.get('causal_decomposition',{}).get('dominant_observed_enrichment_stage')),outlet01_count_chain_present=s['outlet01_loss_chain']['raw_O1']==sum(r['point_outlet']=='OUTLET_01' for r in raw),
  stationary_rank_audited=len(st)==s['stationary_count'] and all('contact_redundancy' in r for r in st),stationary_downstream_audited=all(r['downstream']['status']=='RESOLVED' for r in st),
  stationary_virtual_radius_bounded=all(r['reference_jam_reproduced'] and all(t['criteria']['local_only'] for t in r['trials']) for r in st),stationary_all_classified=all(r['classification'] in ['RIGID_GEOMETRIC_JAM','RIGID_MODEL_HIGH_SENSITIVITY','NUMERICAL_OR_GEOMETRIC_UNRESOLVED'] for r in st),
  production_files_unchanged=not changed and sameindex,remote_formal_reference_matches_local=not remote['missing_in_remote_reference'] and not remote['mismatches_against_local_frozen_reference'],executed_radius_source_archived=executed_radius_sha==read(D/'radius_config.json')['diagnostic_source_sha256'],all_eight_figures_have_sources=figures_ok and len(manifest)==8,
  point_no_exit_causes_analytically_audited=not s['point_unresolved_cause_ids'],point_resolution_controls_consistent=resolution_ok,
  audit_tests_passed=tests['tests']>=9 and tests['errors']==tests['failures']==tests['skipped']==0)
 status='SCIENTIFIC_AUDIT_COMPLETE' if all(gates.values()) and not s['unresolved_ids'] else 'SCIENTIFIC_AUDIT_PARTIAL'
 result=dict(status=status,role='SCIENTIFIC_AUDIT_NOT_PRODUCTION_PASS',gates=gates,tests=tests,production_protection=protection,point_no_exit_count=len(s['point_no_exit_ids']),point_no_exit_interpretation=s['point_no_exit_explanation'],stationary_unresolved_ids=s['unresolved_ids'])
 dump(D/'gates.json',result);s['status']=status;s['production_protection']=protection;dump(D/'audit_summary.json',s);print(result)
if __name__=='__main__':main()
