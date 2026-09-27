from pathlib import Path
import importlib.util

def test_resource_plan_accounts_for_refinement(root):
 spec=importlib.util.spec_from_file_location('p9a5_resource_plan',root/'particle_3d/reports/particle9a5_formal_trajectories/scripts/resource_plan.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 benchmark=dict(dt_s=.001,benchmark_cost_rows=[dict(particle_id=1,bytes=1000000,accepted_steps=1000)],worst_case_5000_runtime_seconds=100.)
 result=dict(rows=[dict(particle_id=1,trajectory_age_s=.1)])
 p=m.make_plan(benchmark,result,1000*10**9)
 # 1000 accepted substeps cover 100 nominal intervals. Extrapolate bytes over
 # actual physical age, including all those substeps, never just 12000 samples.
 assert p['worst_case_5000_storage_bytes']==750000000000
 assert p['resource_limited_N_max']==5000 and p['declared_before_production']
