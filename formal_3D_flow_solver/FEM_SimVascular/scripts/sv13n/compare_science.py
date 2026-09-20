"""Compare native (unshifted pressure) P1 fields with exact tetra mass integration."""
raise SystemExit('NOT RUN BY USER DECISION: scientific equivalence deferred')
import json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv13n import science_gate,SCIENCE_LIMITS
R=ROOT/'reports/sv1_3n';O=ROOT/'outputs/sv1_3n'
label,left,right=sys.argv[1:]
source=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
m=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',source['accepted_solution']['Q_target_m3_s'],source['policy']['Umean_m_s'])
records=[json.loads((R/(n+'_acceptance.json')).read_text()) for n in (left,right)]
N=json.loads((ROOT/'configs/sv1_3n/convergence_strategy.json').read_text())['SHORT_REGRESSION_STEP']
assert all(r['accepted'] and r['steps_completed']==N and r['initial_state']=='t=0' for r in records)
import xml.etree.ElementTree as ET
def science_xml(record):
 root=ET.parse(O/record.get('origin_case',record['name'])/'solver.xml').getroot()
 root.find('.//Number_of_time_steps').text='MATCHED_SHORT'
 def canonical(e):return (e.tag,tuple(sorted(e.attrib.items())),(e.text or '').strip(),tuple(canonical(c) for c in e))
 return canonical(root)
assert science_xml(records[0])==science_xml(records[1]),'SCIENTIFIC_INPUT_MISMATCH'
paths=[O/Path(r['results'][-1]['path']).relative_to('outputs') for r in records]
(u0,p0),(u1,p1)=[m.read(p) for p in paths]
f0,f1=[m.measure(u,p) for u,p in ((u0,p0),(u1,p1))]
relative=lambda a,b:abs(b-a)/max(abs(a),np.finfo(float).tiny)
out={k:relative(a,f1['outlet_flows_m3_s'][k]) for k,a in f0['outlet_flows_m3_s'].items()}
errors=dict(velocity_relative_volume_L2=m.velocity_l2(u1-u0)/m.velocity_l2(u0),pressure_relative_volume_L2=m.velocity_l2((p1-p0)[:,None])/m.velocity_l2(p0[:,None]),relative_Qin=relative(f0['Q_in_m3_s'],f1['Q_in_m3_s']),max_relative_Qout=max(out.values()),mass_error_difference=abs(f1['epsilon_mass']-f0['epsilon_mass']),relative_max_velocity=relative(float(np.linalg.norm(u0,axis=1).max()),float(np.linalg.norm(u1,axis=1).max())))
d=dict(status='PENDING',cases=[left,right],pressure_shift_applied=False,errors=errors,each_Qout_relative=out,thresholds=SCIENCE_LIMITS,measurements=[f0,f1],field_paths=[str(p.relative_to(ROOT)) for p in paths],field_sha256=[r['results'][-1]['sha256'] for r in records])
d['matched_step']=N
try:science_gate(d);d['status']='PASS'
except ValueError as e:d.update(status='FAIL',reason=str(e))
(R/(label+'.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');print(json.dumps(d,indent=2));raise SystemExit(d['status']!='PASS')
