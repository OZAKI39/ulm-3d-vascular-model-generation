"""Real-source replay without admission or trajectories; historical data read-only."""
from pathlib import Path
import json,gzip,hashlib
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.source import PerfusionSource
from particle_3d.population_v1.fixed_window import validate_contract
from particle_3d.population_v1.provenance import canonical,write_json
root=Path(__file__).resolve().parent;pop=root.parent
c=json.loads((root/'FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json').read_text());env=load_environment(json.loads((pop/'data/environment.json').read_text()))
s=PerfusionSource(env,concentration_m3=c['concentration_m3'],seed=c['master_seed'],realization_id=c['realization_id']);validate_contract(c,s)
rows=[]
while (r:=s.peek()) is not None and r['birth_time_s']<c['T_infusion_s']:
 rows.append(dict(r));s.previous_time=r['birth_time_s'];s.next_sid+=1;s.pending=None
s.stop(c['T_infusion_s'])
old=list(map(json.loads,gzip.open(pop/'runs/realization_0000_contactfix/source.jsonl.gz','rt')))
matching=len(rows)<=len(old) and all(canonical({k:v for k,v in r.items() if k!='particle_id'})==canonical({k:previous[k] for k in r if k!='particle_id'}) for r,previous in zip(rows,old))
assert matching
out=root/'validation/real_source_only';out.mkdir(parents=True,exist_ok=True)
with gzip.open(out/'source_proposals.jsonl.gz','wb') as f:
 for row in rows:f.write(canonical(row))
write_json(out/'IDENTITY.json',dict(experiment_contract=c,flow_sha256=env.flow_sha256,scope='SOURCE_ONLY_NO_ADMISSION_NO_TRAJECTORY',current_protocol=c['schema'],previous_protocol='OUTCOME_DEPENDENT_COVERAGE_STOP'))
write_json(out/'SUMMARY.json',dict(PASS=True,actual_source_count=len(rows),expected_source_count=s.rate*c['T_infusion_s'],
 first_birth_s=rows[0]['birth_time_s'],last_birth_s=rows[-1]['birth_time_s'],old_source_count=len(old),
 same_seed_marks_match_historical_prefix=matching,half_open_interval_pass=all(r['birth_time_s']<c['T_infusion_s'] for r in rows),
 accepted_birth_count=None,source_rejection_count=None,reason='No admission or particle trajectories executed'))
print((out/'SUMMARY.json').read_text())
