
import os,json,socket
from pathlib import Path
from datetime import datetime,timezone
KEY_PATHS=['/workspace/formal_3D_flow_solver_FEM_lzy', '/workspace/formal_3D_flow_solver_FEM_lzy/remote/.env', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/outputs/REAL_VASCULAR_GPU_ILU_REUSE_WINNER', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/logs', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/reports', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/svMultiPhysics-reuse', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3n/external/petsc325/install_gpu13', '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran', '/workspace/flow_mean_2p0_mmps_20260922', '/workspace/flow_mean_2p0_mmps_20260922/flow_cases/mean-2p0-mmps', '/workspace/flow_mean_2p0_mmps_20260922/flow_cases/mean-2p0-mmps/run/solver.log', '/root/particle8_2_runs', '/root/particle8_2_runs/env', '/root/particle8_2_runs/env/bin/python', '/root/particle8_2_runs/readonly_inputs', '/root/particle8_2_runs/readonly_inputs/p81/particle_3d/outputs/particle8_1', '/root/particle8_2_runs/2ec91fe4c6d142400a352bb9d949800e0cb3b489/final_execution/repo', '/root/particle8_2_runs/cbb45d8222b3569bce4bb4e2034f1c444b18584e/full_audit', '/root/particle8_2_runs/a4820609c9ed71f09bcea968b5819701632a2a1e/diagnostic_pipeline/point_basin_100000', '/workspace/particle8_2a_20260922', '/workspace/particle8_2a_20260922/results', '/workspace/particle8_2a_20260922/results/admission', '/workspace/particle8_2a_20260922/results/trajectories/formal', '/workspace/particle8_2a_20260922/results/ppt', '/workspace/shear_lift_audit_20260923', '/workspace/shear_lift_audit_20260923/final_source/code', '/workspace/shear_lift_audit_20260923/final_source/tests', '/workspace/shear_lift_audit_20260923/input', '/workspace/shear_lift_audit_20260923/results', '/workspace/shear_lift_audit_20260923/results/states', '/home/lzy/projects/sonovue_size_distribution_v0', '/root/particle8_2_runs/readonly_inputs/sonovue_size_distribution_v0', '/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z', '/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z', '/workspace/particle9a4_population_inlet_20260924T233708Z', '/workspace/particle9a4_population_inlet_20260924T233708Z/particle_3d/src/particle_3d', '/workspace/particle9a4_population_inlet_20260924T233708Z/deployment_manifest.json', '/workspace/particle9a4_population_inlet_20260924T233708Z/smoke30.log', '/workspace/particle9a4_population_inlet_20260924T233708Z/particle_3d/reports/particle9a4_population_inlet', '/workspace/particle9a4_population_inlet_20260924T233708Z/particle_3d/reports/particle9a4_population_inlet/outputs/smoke30']
prefixes=('flow_','formal_','particle','shear_','hemocell','lammps','microbubble','reference','archives')
roots=[p for p in sorted(Path('/workspace').iterdir()) if p.is_dir() and p.name.startswith(prefixes)]
roots += [Path('/root/particle8_2_runs'),Path('/home/lzy/projects/sonovue_size_distribution_v0')]
stop={'.git','.venv','env','.env','__pycache__','node_modules','external','site-packages','build','lib','include'}
dirs=[]
for root in roots:
 if not root.exists():continue
 for base,children,files in os.walk(root,followlinks=False):
  p=Path(base);depth=len(p.relative_to(root).parts)
  dirs.append(dict(root=str(root),path=str(p),depth=depth,direct_file_count=len(files),child_directory_count=len(children)))
  children[:]=sorted(n for n in children if n not in stop and not n.startswith('.') and depth<5)
checked={p:dict(exists=Path(p).exists(),kind='directory' if Path(p).is_dir() else 'file' if Path(p).is_file() else 'missing') for p in KEY_PATHS}
rbc=[]
for root in roots:
 if root.name.startswith(('particle','hemocell','lammps','microbubble')):
  for base,children,files in os.walk(root,followlinks=False):
   p=Path(base);depth=len(p.relative_to(root).parts)
   children[:]=sorted(n for n in children if n not in stop and not n.startswith('.') and depth<8)
   if 'rbc' in p.name.lower():rbc.append(str(p))
   rbc.extend(str(p/n) for n in files if 'rbc' in n.lower() and n.endswith(('.py','.md','.json','.csv','.npz','.h5','.log')))
print(json.dumps(dict(checked_utc=datetime.now(timezone.utc).isoformat(),hostname=socket.gethostname(),roots=list(map(str,roots)),directories=dirs,key_paths=checked,rbc_related_paths=sorted(set(rbc)),scan_scope='directories depth<=5; RBC names depth<=8; no env/vendor/build traversal; path metadata only')))
