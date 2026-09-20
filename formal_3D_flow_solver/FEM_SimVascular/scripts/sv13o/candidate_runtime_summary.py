"""Report actual PCView semantics even for rejected linear solves; no CFD."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13o import parse_runtime_semantics,semantics_gate
R=ROOT/'reports/sv1_3o';cases={}
for k in 'ABCD':
 name='PERF_'+k;a=json.loads((R/(name+'_acceptance.json')).read_text())
 log=ROOT/'logs/sv1_3o/remote'/(name+'.log');rows=parse_runtime_semantics(log.read_text())
 assert len(rows)==a['statistics']['KSP_solves']
 for s in rows:semantics_gate(s,a['PETSC_OPTIONS'])
 cases[name]=dict(solver_status=a['status'],runtime_configuration_verified=True,views=rows,
     statistics=a['statistics'],sampled_device_memory_max_MiB=a['sampled_device_memory_max_MiB'],
     memory_scope='Device-wide samples every15 s, not exact process peak; short failed runs can finish before meaningful sampling',
     KSP_reason_values=a['KSP_reason_values'],factor_packages=a['factor_packages'],log=str(log.relative_to(ROOT)))
(R/'candidate_runtime_summary.json').write_text(json.dumps(dict(status='PASS',cases=cases,
    scope='Configuration verification only. C/D solver failures remain FAIL and excluded from ranking.'),indent=2)+'\n')
print('A-D actual runtime configurations verified; unhealthy candidates remain rejected.')
