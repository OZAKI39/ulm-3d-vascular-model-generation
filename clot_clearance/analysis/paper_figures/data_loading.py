"""Load existing evidence only. No imports of solvers or physical kernels."""
import csv,json,hashlib
from functools import lru_cache
from pathlib import Path
import numpy as np
from figure_specs import V,R,B,A,C5,C7,ONSET,S5,S7,BS
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
USED=set()
def record(path):
 p=ROOT/path
 if not p.is_file():raise FileNotFoundError(f'Required real data missing: {p}')
 USED.add(path);return p
@lru_cache(maxsize=None)
def _json(path):return json.loads((ROOT/path).read_text())
def read_json(path):
 record(path)
 return _json(path)
def read_csv(path):
 with record(path).open(newline='') as f:return list(csv.DictReader(f))
def numeric(rows,key):return np.asarray([float(r[key]) for r in rows])
def table(name,rows):
 p=HERE/'output/tables'/name
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def onset():
 c=read_json(ONSET);rows=[]
 for label,casepath,bracket in [('C5',BS[0],V+'pfpd_at1_c5_activation_stationarity_audit_v1/activation/BRACKET.json'),('C7',BS[1],V+C7+'activation/BRACKET.json')]:
  d=c[label];geo=read_json(casepath);interval=read_json(bracket)['interval'];n=d['count'];total=geo['bonds']
  assert np.isclose(d['fraction'],n/total,rtol=1e-13)
  assert interval[0]<=d['strain']<=interval[1]
  depth=d['nearest_clamp_percentiles_m'][2]
  rows.append(dict(configuration=label,strain=d['strain'],bracket_lower=interval[0],bracket_upper=interval[1],active_count=n,total_bonds=total,active_fraction=n/total,median_depth_m=depth,depth_over_delta=depth/geo['delta_m'],depth_over_h=depth/geo['h_m'],h_m=geo['h_m'],source=ONSET))
 table('activation_onset.csv',rows);return rows

def branches():
 c5=read_json(S5)['states'];c7=read_json(S7)['states'];rows=[]
 for level,states in [('C5',c5),('C7',c7)]:
  for i,d in enumerate(states):
   met=d['metrics'];positive=met['maximum_s']>1e-8
   if level=='C5':
    assert d['status']=='STATIONARY_STABLE'
    stress=None;src=S5
    if positive:
     lab=d['label'];src=V+C5+(f'continuation/{lab}/SUMMARY.json' if lab.startswith('E') else f'states/{lab}/INPUT_VERIFICATION.json')
     x=read_json(src);assert abs(x['metrics']['maximum_s']-met['maximum_s'])<1e-12;stress=x['reaction_stress_Pa']
    count=d['counts']['INTERIOR'];stable='AUDITED_STABLE'
   else:
    assert d['scientifically_accepted'];stress=d['nominal_stress_Pa'];count=d['active_count'];src=S7;stable='NOT_ASSESSED'
   rows.append(dict(configuration=level,index=i,label=d['label'],strain=d['strain'],maximum_s=met['maximum_s'],weighted_mean_s=met['mean_s'],nominal_stress_Pa=stress,force_residual=met['force_residual'],projected_phase_KKT=met['projected_phase_KKT'],active_count=count,positive=positive,stability=stable,source=src))
 table('branch_states.csv',rows)
 return rows

def boundary():return {d['case']:d for d in [read_json(p) for p in BS]}
def layer(c):return read_csv(R+B+f'tables/{c}_layers.csv')
def cumulative(c):return [r for r in read_csv(R+B+f'tables/{c}_cumulative.csv') if r['side']=='top']
def particles(c):return read_csv(R+B+f'tables/{c}_particles.csv')
def gaps(c):return read_csv(R+A+f'tables/{c}_FULL_CLAMP_BASELINE_0.0030_CENTERLINE_GAPS.csv')
def clamp():return read_json(V+A+'SUMMARY.json')['clamp_comparison']
def source_metadata(paths):
 ans=[]
 for rel in sorted(set(paths)):
  p=record(rel);h=hashlib.sha256()
  with p.open('rb') as f:
   for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
  ans.append(dict(path=rel,bytes=p.stat().st_size,sha256=h.hexdigest()))
 return ans
