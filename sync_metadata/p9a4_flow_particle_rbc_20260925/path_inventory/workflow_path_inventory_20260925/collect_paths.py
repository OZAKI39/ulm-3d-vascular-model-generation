from pathlib import Path
from datetime import datetime,timezone
import os,json,csv,subprocess
W=Path('/home/lzy/projects');OUT=W/'workflow_path_inventory_20260925'
old=json.loads((W/'workflow_path_inventory_20260923/key_paths.json').read_text())
keys=[dict(x) for x in old]
for x in keys:
 x['category']=x['category'].replace('最新 2.0 mm/s 完整流场','历史 equal-pressure 2.0 mm/s 流场').replace('正式粒子所依赖的旧冻结 FEM','历史粒子使用的冻结 FEM')
 x['purpose']=x['purpose'].replace('最新','已有').replace('新流场','历史 equal-pressure 流场')
def add(category,purpose,p,remote=False):
 keys.append(dict(location='root@50.115.148.16:4159' if remote else 'WSL',category=category,purpose=purpose,path=str(p)))
P=W/'ulm_particle_3d_particle0/particle_3d';N=W/'ulm_3D_vascular';F=W/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular';A=W/'ulm_particle_population_inlet_p9a4';R=A/'particle_3d/reports/particle9a4_population_inlet';H=F/'flow_cases/mean-2p0-mmps-A-H0-pressure-v1'
for label,p in [('全新开发worktree',A),('新入口源码',A/'particle_3d/src/particle_3d/continuous_infusion.py'),('NEW流场适配与worker',A/'particle_3d/src/particle_3d/population_inlet_p9a4.py'),('新浓度与源合同',A/'particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json'),('新永久测试',A/'particle_3d/tests/particle9a4_population_inlet'),('P9-A.4全部报告及结果',R),('总预览',R/'OPEN_RESULTS.html'),('100k全部source ledger',R/'data/inlet100k/proposal_ledger.jsonl'),('accepted births',R/'data/inlet100k/accepted_births.json'),('审计数据与机器摘要',R/'data'),('生成、部署、审核与绘图脚本',R/'scripts'),('本地/回传日志',R/'logs'),('30泡与point完整输出',R/'outputs/smoke30'),('30泡轨迹',R/'outputs/smoke30/trajectories'),('对应point轨迹',R/'outputs/smoke30/point'),('8组PNG/PDF图件',R/'figures'),('P9-A.4使用的NEW流场副本',A/'particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/inputs/NEW.vtu')]:add('当前 P9-A.4',label,p)
for label,p in [('vascular A/ROI原始主工程',N),('血管原数据',N/'vessel_model'),('SWC/ROI/表面数据',N/'outputs'),('预处理共享代码',N/'vascular_processing'),('Network 1D/0D模块',N/'network_1d0d'),('H0求解入口',N/'scripts/run_a_network_h0_v2.py'),('H0 FEM构建',N/'network_1d0d/fem_h0_case.py'),('H0远程求解',N/'scripts/solve_a_h0_fem_remote.py'),('Network v1报告',N/'reports/a_network_1d0d_boundary_v1'),('H0 v2报告/数据/日志',N/'reports/a_network_1d0d_boundary_v2_idealized'),('NEW H0完整FEM算例',H),('NEW H0 frozen_flow',H/'frozen_flow'),('NEW H0最终场',H/'frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu'),('NEW H0网格',H/'SV_MESH'),('NEW H0原生运行/日志/检查点',H/'run'),('NEW H0求解配置',H/'run/solver.xml'),('NEW H0验收与通量报告',H/'reports'),('共享FEM测量代码',F/'src/sv_validation'),('已停止Taylor-Hood审核',F/'reports/taylor_hood_p2p1_validation_v1')]:add('当前 vascular A → Network H0 → FEM',label,p)
for q in sorted((F/'flow_cases').glob('*taylor*')):add('Taylor-Hood 已停止的历史尝试',q.name,q)
B=P/'reports/network_derived_flow_mb_validation_v1'
for label,p in [('Network OLD/NEW配对验证',B),('配对runner/analyze',B/'scripts'),('OLD/NEW输入与源码快照',B/'server_bundle'),('流场身份与配对事件',B/'data'),('配对轨迹输出',B/'outputs'),('配对日志',B/'logs'),('P9-A.1动力学',P/'src/particle_3d/particle9a_motion.py')]:add('Network-H0 配对验证与冻结动力学',label,p)
for name in ['particle9a_2mmps','particle9a_2mmps_diagnosis','particle9a1_2mmps','particle9a1_routing_stationary_audit','particle9a2_inlet_sampling','particle9a3_interior_inlet','particle9a3b_flowfield_conservation']:
 add('P9历史验证报告',name,P/'reports'/name)
for q in sorted((P/'outputs').iterdir()):
 if q.is_dir():add('Particle 各阶段原始输出',q.name,q)
for q in sorted((P/'src/particle_3d').glob('*rbc*.py')):add('RBC 源码与验证入口',q.name,q)
for name in ['generate_rbc_population.py','run_particle2_validation.py','run_particle3_real_validation.py','prepare_rbc_mb_flow_rotation.py','render_rbc_mb_flow_rotation.py','verify_rbc_reproducibility.py']:
 add('RBC 源码与验证入口',name,P/'scripts'/name)
for q in sorted((P/'tests').glob('rbc*')):add('RBC 源码与验证入口','测试 '+q.name,q)
for part in ['particle2','particle3']:
 for name in ['data','figures','logs']:add('RBC 运动、形状与展示',part+' '+name,P/'reports'/part/name)
for q in sorted((W/'archives').glob('*')):pass
for q in sorted(Path('/home/lzy/archives').iterdir()):add('清理、同步与部署归档',q.name,q)
for name in ['network_h0_particle_20260925','flow_microbubble_rbc_20260923','flow_mb_rbc_sync_20260923_audit','ulm_sonovue_20260920_111258']:add('Git 同步快照',name,W/'github_sync'/name)
G=W/'github_sync/network_h0_particle_20260925'
for name in ['vascular_network','particle_3d','formal_3D_flow_solver','sonovue_size_distribution_v0','source_dependencies','server_evidence','sync_metadata/network_h0_particle_20260925','THESIS_SUMMARY']:add('最新已同步快照内容',name,G/name)
for name in ['local_file_inventory.csv','server_file_inventory.json','original_path_map.json','SYNC_REPORT_ZH.md']:add('最新已同步逐文件清单',name,G/'sync_metadata/network_h0_particle_20260925'/name)
remote_roots=['/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z','/workspace/particle_network_flow_mb_validation_v1_20260924T210047Z','/workspace/particle9a4_population_inlet_20260924T233708Z']
for p in remote_roots:add('服务器当前主线',Path(p).name,p,True)
rr=remote_roots[-1]
for label,rel in [('P9-A.4源代码','particle_3d/src/particle_3d'),('P9-A.4部署manifest','deployment_manifest.json'),('P9-A.4终端日志','smoke30.log'),('P9-A.4审计入口、数据及结果','particle_3d/reports/particle9a4_population_inlet'),('30泡与point输出','particle_3d/reports/particle9a4_population_inlet/outputs/smoke30')]:add('服务器当前主线',label,rr+'/'+rel,True)
for x in keys:
 if x['location']=='WSL':
  p=Path(x['path']);x['exists']=p.exists();x['kind']='directory' if p.is_dir() else 'file' if p.is_file() else 'missing';x['symlink']=p.is_symlink()
remote_paths=[x['path'] for x in keys if x['location']!='WSL']
remote_code=r'''
import os,json,socket
from pathlib import Path
from datetime import datetime,timezone
KEY_PATHS=__KEY_PATHS__
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
'''.replace('__KEY_PATHS__',repr(remote_paths))
(OUT/'remote_scan.py').write_text(remote_code)
proc=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15','vast4090','python3','-'],input=remote_code,text=True,capture_output=True,timeout=120)
(OUT/'remote_scan_stderr.log').write_text(proc.stderr)
if proc.returncode:raise RuntimeError(proc.stderr)
remote=json.loads(proc.stdout);(OUT/'remote_paths.json').write_text(json.dumps(remote,indent=2))
for x in keys:
 if x['location']!='WSL':x.update(remote['key_paths'][x['path']])
for p in remote['roots']:
 if p not in remote_paths:keys.append(dict(location='root@50.115.148.16:4159',category='服务器历史运行/独立模拟工程',purpose=Path(p).name,path=p,exists=True,kind='directory'))
stop={'.git','.venv','env','.env','__pycache__','node_modules','external','site-packages','build','vessel_model','.cache'}
roots=[p for p in sorted(W.iterdir()) if p.is_dir() and (p.name.startswith(('ulm','formal','sonovue','vascular','CLEANUP')) or p.name=='github_sync')]+[Path('/home/lzy/archives')]
dirs=[]
for root in roots:
 for base,children,files in os.walk(root,followlinks=False):
  p=Path(base);depth=len(p.relative_to(root).parts)
  dirs.append(dict(root=str(root),path=str(p),depth=depth,direct_file_count=len(files),child_directory_count=len(children)))
  children[:]=sorted(n for n in children if n not in stop and not n.startswith('.') and depth<5)
local=dict(checked_utc=datetime.now(timezone.utc).isoformat(),roots=list(map(str,roots)),directories=dirs,scan_scope='directories depth<=5, no Git/env/vendor/raw dataset traversal; top-level roots and key paths retained')
(OUT/'local_paths.json').write_text(json.dumps(local,indent=2))
unique={}
for x in keys:unique.setdefault((x['location'],x['category'],x['path']),x)
keys=list(unique.values());(OUT/'key_paths.json').write_text(json.dumps(keys,indent=2,ensure_ascii=False))
for name,rows in [('local_directories',dirs),('remote_directories',remote['directories'])]:
 with (OUT/(name+'.csv')).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=['root','path','depth','direct_file_count','child_directory_count']);wr.writeheader();wr.writerows(rows)
with (OUT/'key_paths.csv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=['location','category','purpose','path','exists','kind'],extrasaction='ignore');wr.writeheader();wr.writerows(keys)
(OUT/'git_worktrees.txt').write_text(subprocess.check_output(['git','-C',str(W/'ulm-3d-vascular-model-generation'),'worktree','list','--porcelain'],text=True))
print(json.dumps(dict(key_paths=len(keys),local_directories=len(dirs),remote_directories=len(remote['directories']),remote_roots=remote['roots'],missing=[x for x in keys if not x['exists']],remote_rbc_paths=len(remote['rbc_related_paths'])),indent=2,ensure_ascii=False))
