from pathlib import Path
import json,sys,math,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.continuous_infusion import canonical_bytes,sha256
from particle_3d.injection_method_c import TruncatedSonoVue
from particle_3d.particle82a_geometry import maximum_handoff_radius
R=ROOT/'particle_3d/reports/particle9a4_population_inlet';D=R/'data'
# Fixed broad audit tolerances, recorded before evaluating the distribution gates.
policy=dict(source_cdf_absolute_tolerance=.01,enter_cdf_absolute_tolerance=.03,
 joint_cell_probability_absolute_tolerance=.02,source_position_cdf_absolute_tolerance=.02,
 entering_position_cdf_absolute_tolerance=.035,rate_combined_standard_errors=4.,
 source_interarrival_mean_relative_tolerance=.02,source_interarrival_variance_relative_tolerance=.04,
 N_min=100000,fano_is_descriptive=True)
if (D/'audit_tolerances.json').exists():
 assert (D/'audit_tolerances.json').read_bytes()==canonical_bytes(policy)
else:
 (D/'audit_tolerances.json').write_bytes(canonical_bytes(policy))
rows=[json.loads(line) for line in (D/'inlet100k/proposal_ledger.jsonl').read_text().splitlines()]
s=json.loads((D/'inlet100k/summary.json').read_text());q=json.loads((D/'qacc_diagnostic.json').read_text());c=json.loads((ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json').read_text());dist=TruncatedSonoVue(ROOT/'sonovue_size_distribution_v0')
a=np.array([r['particle_id'] is not None for r in rows]);diam=np.array([r['diameter_m'] for r in rows]);xyz=np.array([r['position_m'] for r in rows]);g=np.load(D/'inlet_geometry.npz');z=(xyz-g['origin'])@g['basis'].T
control=np.load(D/'qacc_positions.npz');control_z=(control['position_m']-g['origin'])@g['basis'].T;maxD=2*maximum_handoff_radius(control['wall_distance_m']+q['wall_roundoff_m']);prob=dist.cdf(maxD)
grid=np.linspace(.5e-6,4e-6,71)
source_error=max(abs(np.mean(diam<=d)-dist.cdf(d)) for d in grid)
enter_error=max(abs(np.mean(diam[a]<=d)-np.mean(dist.cdf(np.minimum(d,maxD)))/prob.mean()) for d in grid)
spatial=[]
for axis in range(2):
 cuts=np.quantile(control_z[:,axis],np.linspace(.1,.9,9))
 spatial.append(dict(axis=axis,source_max_cdf_error=max(abs(np.mean(z[:,axis]<=cut)-np.mean(control_z[:,axis]<=cut)) for cut in cuts),
  entering_max_cdf_error=max(abs(np.mean(z[a,axis]<=cut)-np.mean(prob*(control_z[:,axis]<=cut))/prob.mean()) for cut in cuts)))
joint=[];d_edges=np.array([0.,1.2,1.8,2.4,4.])*1e-6
for quadrant in range(4):
 obsmask=(z[:,0]>0)==bool(quadrant&1);obsmask&=(z[:,1]>0)==bool(quadrant&2)
 predmask=(control_z[:,0]>0)==bool(quadrant&1);predmask&=(control_z[:,1]>0)==bool(quadrant&2)
 for lo,hi in zip(d_edges[:-1],d_edges[1:]):
  observed=float(np.mean((obsmask&(diam>lo)&(diam<=hi))[a]))
  w=np.maximum(0.,dist.cdf(np.minimum(maxD,hi))-dist.cdf(np.minimum(maxD,lo)))
  expected=float(np.mean(w*predmask)/prob.mean());joint.append(dict(quadrant=quadrant,D_lo_um=lo*1e6,D_hi_um=hi*1e6,observed=observed,expected=expected,absolute_error=abs(observed-expected)))
p_expected=q['rate_estimate']['acceptance_probability'];se_diagnostic=(q['rate_estimate']['ci95_probability'][1]-q['rate_estimate']['ci95_probability'][0])/(2*1.96);se_sample=math.sqrt(p_expected*(1-p_expected)/len(rows));zrate=abs(s['acceptance_fraction']-p_expected)/math.hypot(se_diagnostic,se_sample)
checks=dict(count=s['N_source_proposals']>=policy['N_min'],all_events=sum([s['N_accepted'],s['N_rejected']])==len(rows),
 source_size=source_error<policy['source_cdf_absolute_tolerance'],enter_size=enter_error<policy['enter_cdf_absolute_tolerance'],
 joint=max(r['absolute_error'] for r in joint)<policy['joint_cell_probability_absolute_tolerance'],
 source_position=all(r['source_max_cdf_error']<policy['source_position_cdf_absolute_tolerance'] for r in spatial),
 entering_position=all(r['entering_max_cdf_error']<policy['entering_position_cdf_absolute_tolerance'] for r in spatial),
 rate=zrate<policy['rate_combined_standard_errors'],
 timing_mean=abs(s['source_mean_interarrival_s']*s['lambda_source_s_inv']-1)<policy['source_interarrival_mean_relative_tolerance'],
 timing_variance=abs(s['source_variance_interarrival_s2']*s['lambda_source_s_inv']**2-1)<policy['source_interarrival_variance_relative_tolerance'],
 only_wall_rejections=set(s['rejection_reasons'])<={'WALL_REJECTED','WALL_NEARFIELD_REJECTED'},
 no_retry=all(r['position_draw_count']==r['diameter_draw_count']==r['admission_check_count']==1 for r in rows),
 monotone_Qacc=all(x['fraction']>=y['fraction'] for x,y in zip(q['rows'][:-1],q['rows'][1:])))
checks={k:bool(v) for k,v in checks.items()}
evidence=dict(checks=checks,source_size_cdf_error=float(source_error),entering_size_cdf_error=float(enter_error),spatial=spatial,joint_cells=joint,rate_discrepancy_in_combined_standard_errors=zrate)
with (D/'inlet_audit_gate_evidence.json').open('xb') as f:f.write(canonical_bytes(evidence))
assert all(checks.values()),checks
reg=json.loads((D/'presmoke_regression.json').read_text())
def suite(name):
 x=ET.parse(R/'logs'/name).getroot().find('testsuite');return {k:int(x.attrib[k]) for k in ['tests','failures','errors','skipped']}
new=suite('new_pass1.xml');legacy=suite('legacy_b_c.xml')
gate=dict(inlet_audit_pass=True,portable_regression_pass=reg['tests']==115 and reg['failures']==0,new_tests_pass=new['failures']==new['errors']==0,legacy_b_c_pass=legacy['failures']==legacy['errors']==0,
 smoke30_births_sha256=sha256(D/'inlet100k/smoke30_births.json'),contract_sha256=sha256(ROOT/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json'),portable_tests=115,new_tests=new,legacy_tests=legacy,production_500_authorized=False)
assert all(gate[k] for k in ['portable_regression_pass','new_tests_pass','legacy_b_c_pass'])
with (D/'smoke_gate.json').open('xb') as f:f.write(canonical_bytes(gate))
print(json.dumps(dict(gate=gate,evidence={k:v for k,v in evidence.items() if k!='joint_cells'}),indent=2))
