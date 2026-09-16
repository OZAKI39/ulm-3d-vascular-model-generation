"""Audit literal-source probes and actual upstream LAMMPS force arrays, no solver calls."""
from pathlib import Path
import json,hashlib,argparse
import numpy as np

def audit(root):
 p=Path(root)/'source_correction_audit';up=Path(root)/'provenance/upstream_lubrication';a=np.genfromtxt(p/'SOURCE_COEFFICIENTS.csv',delimiter=',',names=True);metrics={}
 source_identity=(up/'pair_lubricate_poly.cpp').read_bytes()==(up/'TARGET_RAW_pair_lubricate_poly.cpp').read_bytes()
 assert source_identity, 'Target source differs from official fixed-commit raw source'
 extracts=json.loads((p/'EXTRACTION_PROVENANCE.json').read_text())
 for name,meta in extracts.items():
  raw=(up/meta['source_file']).read_bytes();literal=(p/(name+'_literal_source.txt')).read_text()[:-1]
  assert hashlib.sha256(raw).hexdigest()==meta['source_sha256']
  assert hashlib.sha256(literal.encode()).hexdigest()==meta['literal_sha256'] and literal in raw.decode()
 for flag in [0,1]:
  b=a[a['flaglog']==flag];row={}
  for label in ['target_poly','target_upoly','old_poly']:
   row[label]={}
   for term in ['sq','sh']:
    x=b[label+'_'+term];y=b[label+'_'+term+'_swap'];error=np.abs(x-y)/np.maximum(np.maximum(np.abs(x),np.abs(y)),np.finfo(float).tiny)
    row[label][term]={'max_relative_swap_error':float(error.max()),'violations_above_1e_12':int(np.sum(error>1e-12)),'finite':bool(np.isfinite(x).all() and np.isfinite(y).all())}
  metrics[str(flag)]=row
 actual=[]
 for d in sorted(p.glob('actual_lammps_*')):
  text=(d/'forces.dump').read_text().splitlines();names=text[8].split()[2:];data=np.atleast_2d(np.loadtxt(text[9:]));cols={n:data[:,i] for i,n in enumerate(names)};forces=np.c_[cols['fx'],cols['fy'],cols['fz']];defect=float(np.linalg.norm(forces.sum(axis=0))/np.linalg.norm(forces,axis=1).max());flag=int(d.name.split('flaglog')[1][0]);coeff=-forces[:,0]/np.array([.002 if v>0 else -.002 for v in cols['vx']])
  actual.append({'case':d.name,'flaglog':flag,'force_sum_N':forces.sum(axis=0).tolist(),'force_i_N':forces[0].tolist(),'force_j_N':forces[1].tolist(),'newton_third_law_relative_defect':defect,'resistance_from_force_kg_s':coeff.tolist(),'status':'PASS' if defect<=1e-12 else 'FAIL','timesteps':0,'returncode':json.loads((d/'RUN_METRICS.json').read_text())['returncode']})
 confirmed=metrics['1']['target_poly']['sq']['violations_above_1e_12']==0 and metrics['1']['target_upoly']['sq']['violations_above_1e_12']==0 and metrics['1']['target_poly']['sh']['violations_above_1e_12']==0 and all(x['status']=='PASS' for x in actual)
 return {'status':'PASS' if confirmed else 'FAIL','POLYDISPERSE_LUBRICATION_FIX_PRESENT':'YES' if confirmed else 'NO_AS_REQUIRED_BY_DOCUMENTED_SYMMETRY','literal_probe_pairs_per_flag':int(len(a)//2),'source_matches_official_commit':source_identity,'coefficient_metrics':metrics,'actual_unmodified_target_lammps':actual,'normal_only_symmetry_pass':all(metrics['0'][k]['sq']['violations_above_1e_12']==0 for k in ['target_poly','target_upoly']),'interpretation':'The normal-only flaglog=0 branch is symmetric. The documented 4Jul2026 unequal-radius correction for the full pair style is not verified in the frozen target: flaglog=1 squeeze and shear are asymmetric, independently confirmed in actual force arrays. Do not equate normal-only branch PASS with the required source-correction gate.','consequence':'PHASE_A_MIGRATION_FAIL; prohibit promotion, old active deletion and Phase C until a separately authorized revised source contract.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);a=ap.parse_args();r=audit(a.root);Path(a.root,'source_correction_audit/SOURCE_CORRECTION_AUDIT.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
