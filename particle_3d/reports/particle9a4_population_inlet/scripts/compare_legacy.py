"""Inlet-only legacy B/C comparison on NEW flow; never used for acceptance."""
from pathlib import Path
from types import SimpleNamespace
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.population_inlet_p9a4 import load_new_environment,REPORT_RELATIVE
from particle_3d.injection_method_c import MethodCSource,canonical_bytes
from particle_3d.particle82a_admission import common_event,method_b
from particle_3d.particle82a_geometry import InletGeometryAudit
R=ROOT/REPORT_RELATIVE;t=time.perf_counter();env=load_new_environment(ROOT)
ctx=SimpleNamespace(env=env,geometry=InletGeometryAudit(env),distribution=env.distribution.original)
B=[]
for sid in range(1,501):
 event=common_event(sid,master_seed=2026092596,ctx=ctx);result=method_b(event,ctx)
 B.append(dict(event=event,result=result))
 print('B',sid,flush=True) if sid%100==0 else None
c=MethodCSource(sampler=env.sampler,wall=env.wall,field=env.field,distribution=env.distribution,
 flow_sha256=env.new_flow_sha256,geometry_sha256=env.wall.provenance['sha256'],seed=2026092596)
C=[]
for sid in range(1,501):
 C.append(c.event(sid))
 print('C',sid,flush=True) if sid%100==0 else None
for name,rows in [('legacy_b_inlet_500.json',B),('legacy_c_inlet_500.json',C)]:
 with (R/'data'/name).open('xb') as f:f.write(canonical_bytes(rows))
summary=dict(N_legacy_requested=500,flow_sha256=env.new_flow_sha256,runtime_s=time.perf_counter()-t,
 B=dict(accepted=sum(x['result']['accepted'] for x in B),draws=sum(x['result']['attempt_count'] for x in B),source_distribution='ORIGINAL_FULL_SONOVUE',retry_size=True,retry_position=False,physical_rejected_events_preserved=False,rate_semantics='LEGACY_DETERMINISTIC_ENTERING_CLOCK; RETRIES_NOT_PHYSICAL_ARRIVALS'),
 C=dict(accepted=len(C),source_draws=sum(len(x['source_diameter_draws']) for x in C),position_draws=sum(x['position_proposal_count'] for x in C),source_distribution='SONOVUE_D_LE_4UM_CONDITIONAL',retry_size='WHEN_GLOBALLY_IMPOSSIBLE',retry_position=True,physical_rejected_events_preserved=False,rate_semantics='DETERMINISTIC_C_MB_TIMES_Q_ENTERING_CLOCK'),
 warning='Different source-size contracts are shown explicitly; B/C counts are requested events, not physical Poisson proposals. No outlet quota or split gate.')
with (R/'data/legacy_comparison_summary.json').open('xb') as f:f.write(canonical_bytes(summary))
print(json.dumps(summary),flush=True)
