#!/usr/bin/env python3
"""Compare distributed integrals and aggregate new-process restart evidence."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import write_json,timestamp
base=ROOT/'outputs/stage02/cases'
def qc(case,name): return json.loads((base/case/'qc'/f'{name}.json').read_text())
one,two='pipe_medium_natural','mpi2_medium'
a,b=qc(one,'analytic_comparison'),qc(two,'analytic_comparison')
keys=['lambda_pa','velocity_L2_norm_m_pow_2p5_s','pressure_integral_pa_m3','velocity_L2_relative_error','velocity_global_L2_relative_error','pressure_profile_error','lambda_relative_error']
diffs={key:abs(a[key]-b[key])/abs(a[key]) for key in keys}
fa,fb=qc(one,'flux'),qc(two,'flux')
for key in ['actual_Q_in_m3_s','actual_Q_out_m3_s']:
    diffs[key]=abs(fa[key]-fb[key])/abs(fa['target_Q_m3_s'])
assert fa['status']==fb['status']=='PASS'
sa,sb=qc(one,'solver'),qc(two,'solver')
assert sa['mpi_ranks']==1 and sb['mpi_ranks']==2
assert sa['real_global_dofs']==sb['real_global_dofs']==sum(sb['real_owned_dofs_per_rank'])==1
assert max(diffs.values())<=1e-9
write_json(ROOT/'reports/stage02/mpi_reproducibility.json',{'timestamp':timestamp(),'status':'PASS','ranks':[1,2],'cases':[one,two],'real_global_dofs':1,'real_owned_dofs_per_rank':sb['real_owned_dofs_per_rank'],'tolerance':1e-9,'relative_differences':diffs,'comparison':'Physical integrals and analytic diagnostics; file hashes are not an MPI equality criterion','one_rank':{k:a[k] for k in keys},'two_ranks':{k:b[k] for k in keys}})
results={p.name:qc(p.name,'result_reload') for p in base.iterdir() if (p/'qc/result_reload.json').is_file()}
assert len(results)==12 and all(d['status']=='PASS' for d in results.values())
write_json(ROOT/'reports/stage02/result_roundtrip.json',{'timestamp':timestamp(),'status':'PASS','case_count':len(results),'format':'Exact native P2/P1/Real restart from nodal NPZ + profile XDMF/HDF5 geometry; derived DG0 NPZ checked against recomputation. Visualization XDMF/HDF5 mesh independently reloaded; no DOLFINx read_function API is assumed.','cases':results})
print('MPI and 12-case new-process reload: PASS',diffs)
