"""Exact H0 mesh/BC: CPU8 recheck and separate dt/2 numerical-uncertainty check."""
import shutil,json
from case_common import *
base=V/'stage3/vessel_baseline';assert (base/'reports/baseline_reuse.json').exists()
for name,scale in [('vessel_baseline_cpu_mpi8',1.),('vessel_baseline_cpu_mpi8_halfdt',.5)]:
 case=V/'stage3'/name;assert not case.exists();case.mkdir();shutil.copytree(base/'SV_MESH',case/'SV_MESH');pol=vascular_policy(name,'H0_same_mesh_backend_recheck_CFD' if scale==1 else 'H0_fixed_mesh_time_step_CFD');pol.update(mesh_reference='vessel_baseline',dt_s=pol['dt_s']*scale,time_step_scale_relative_to_H0=scale)
 if scale!=1:pol['dt_rule']='Separate dt/2 experiment on original H0 mesh; excluded from main spatial mesh sequence'
 dump(case/'policy.json',pol);xml_case(case,pol['dt_s']);opts=case/'run/PETSC_OPTIONS.txt';opts.write_text(opts.read_text().replace('-mat_type aijcusparse -vec_type cuda','-mat_type aij -vec_type standard'));lock_case(case);print(case)
