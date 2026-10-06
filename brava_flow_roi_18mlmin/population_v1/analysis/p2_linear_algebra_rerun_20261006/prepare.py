from common import *
import shutil, xml.etree.ElementTree as ET
rt=HERE/'runtime';rt.mkdir(exist_ok=True)
for name,source in [('p1_candidate','p1_schur'),('p2_system','p2_system')]:
 case=rt/name;run=case/'run';run.mkdir(parents=True,exist_ok=True)
 if not (case/'SV_MESH').exists(): (case/'SV_MESH').symlink_to(PRIOR/'runtime/p2/SV_MESH' if name.startswith('p2') else P1/'SV_MESH',target_is_directory=True)
 for n in ['solver.xml','inlet_exact_velocity.dat','outlet_exact_traction.vtp']:
  shutil.copy2(PRIOR/'runtime'/source/'run'/n,run/n)
 (run/'PETSC_OPTIONS.txt').write_text(read(HERE/'LINEAR_SOLVER_CANDIDATE_PLAN.json')['candidates'][0]['options']+'\n')
put('PREPARED_INPUT_IDENTITY.json',{str(p.relative_to(HERE)):sha(p) for p in rt.rglob('*') if p.is_file() and p.suffix in {'.xml','.dat','.vtp','.txt'}})
for n in ['P2_DOF_STORAGE_AUDIT.json','COMPACT_MIXED_GRAPH_STATS.json','P2_PETSC_GRAPH_MEASUREMENTS.json','ASM_LU_MEMORY_AUDIT.json','HIGH_ORDER_OUTPUT_GATE.json','P1_RESTART_CORRESPONDENCE.json','P2_NUMERICAL_MATRIX_COUNTS.json']:
 v=read(PRIOR/n);v=v if isinstance(v,dict) else {'data':v};v['inheritance']={'path':str((PRIOR/n).relative_to(BASE)),'sha256':sha(PRIOR/n),'kind':'FROZEN_PRIOR_MEASUREMENT_NOT_RERUN'};put(n,v)
put('INHERITED_EVIDENCE.json',{str(p.relative_to(BASE)):sha(p) for p in PRIOR.iterdir() if p.is_file() and p.suffix in {'.json','.npz','.md','.py','.h'}})
