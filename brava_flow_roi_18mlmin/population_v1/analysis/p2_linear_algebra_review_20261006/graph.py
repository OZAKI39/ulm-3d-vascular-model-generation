from launch import *
import shutil
case=HERE/'runtime/p2_graph_observerfix';run_dir=case/'run';run_dir.mkdir(parents=True,exist_ok=True)
if not (case/'SV_MESH').exists():(case/'SV_MESH').symlink_to(HERE/'runtime/p2/SV_MESH',target_is_directory=True)
for name in ['solver.xml','PETSC_OPTIONS.txt','inlet_exact_velocity.dat','outlet_exact_traction.vtp']:
 shutil.copy2(HERE/'runtime/p2_graph/run'/name,run_dir/name)
run('p2_graph_observerfix','graph')
