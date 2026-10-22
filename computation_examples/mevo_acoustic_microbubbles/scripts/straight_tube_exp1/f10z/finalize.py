"""Protect prior research, inspect artifacts and save a non-self-referential manifest."""
import argparse,json,difflib,re
import numpy as np,h5py
from PIL import Image
from .common import *
def audit():
 reg=verify_inputs();runs=read(BASE/'analysis/RUN_ANALYSIS.json');rows=[];rawhash={}
 for name,r in runs.items():
  p=ROOT/r['path'];raw=p/'trajectory.h5';derived=BASE/'dataset'/(name+'.h5')
  assert r['force_balance_sample_max_N']<=r['max_force_residual_N']+1e-18
  assert sha(raw)==r['trajectory_sha256'];assert sha(derived)==r['dataset_sha256']
  with h5py.File(raw) as a,h5py.File(derived) as b:
   assert np.array_equal(a['positions_m'][:],b['positions_m'][:])
   assert np.array_equal(a['time_s'][:],b['time_s'][:])
   assert b.attrs['flow']=='STRICT_ZERO' and b.attrs['no_padding']
   assert np.array_equal(a['particle_id'][:],b['particle_id'][:])
   assert abs(float(a['time_s'][-1])-r['completed_s'])<1e-14
   dt=a['step_time_duration_segment'][:];assert abs(dt[:,1].sum()-r['completed_s'])<1e-13
   assert b['step_metrics'].shape==a['step_metrics'].shape
  rows.append(dict(name=name,raw_sha256=sha(raw),dataset_sha256=sha(derived),pass_exact_state_and_time=True))
  rawhash.setdefault(r['trajectory_sha256'],[]).append(name)
 old=read(BASE/'evidence/F10S_BASELINE_AUDIT.json')
 for p in (ROOT/'results/straight_tube_exp1/f10s').glob('*/runs/*/RESULT.json'):
  if (p.parent/'trajectory.h5').exists():assert sha(p.parent/'trajectory.h5')==read(p)['trajectory_sha256']
 PDFs=[]
 import pymupdf
 for p in sorted((BASE/'figures').glob('*.pdf')):
  with pymupdf.open(p) as d:
   assert d.page_count==1 and len(d[0].get_text())>100
   preview=BASE/'evidence'/(p.stem+'_pdf.png');d[0].get_pixmap(matrix=pymupdf.Matrix(.85,.85)).save(preview)
   PDFs.append(dict(file=p.name,sha256=sha(p),pages=d.page_count,text_characters=len(d[0].get_text()),preview=str(preview.relative_to(BASE))))
 pngs=list((BASE/'figures').glob('*.png'));assert len(pngs)==len(PDFs)==8
 for p in pngs:
  with Image.open(p) as im:im.verify()
 videos=read(BASE/'evidence/VIDEO_VALIDATION.json');assert len(videos)==3
 for v in videos:
  assert v['fully_decoded'] and v['PTS_monotonic'] and sha(BASE/v['path'])==v['sha256']
  for frame in v['source_mapping']:
   source=frame.get('sources') or [dict(run=frame['run'],saved_index=frame['saved_index'])]
   for s in source:
    with h5py.File(BASE/'dataset'/(s['run']+'.h5')) as h:assert abs(h['time_s'][s['saved_index']]-frame['time_s'])<1e-13
 write(BASE/'evidence/ARTIFACT_AUDIT.json',dict(status='PASS',identity=identity(),protected_sources=len(reg['protected_source_sha256']),old_trajectories=len(old['trajectories']),raw_derived_checks=rows,exact_raw_duplicate_hash_groups=[v for v in rawhash.values() if len(v)>1],PDF_checks=PDFs,png_count=len(pngs),video_count=len(videos),scope='Numerical/data/artifact integrity, not experimental calibration'))
 print('Protected inputs, exact raw-to-derived states, PDFs and frame-source mappings passed.')
def progress():
 path=Path('/home/lzy/projects/microbubble_aggregation.md');before=BASE/'evidence/progress_before.md';backup=Path('/home/lzy/projects/microbubble_aggregation.md.bak_F10Z_20261018')
 old=before.read_bytes();current=path.read_bytes();marker='【F10Z-Simulation，2026-10-18】'
 if marker.encode() in current:
  assert (BASE/'evidence/progress_after.md').read_bytes()==current
  print('Progress already updated and byte-identical.');return
 assert current==old,'Progress document changed since initial audit; reconcile instead of overwrite.'
 if backup.exists():assert backup.read_bytes()==old
 else:backup.write_bytes(old)
 r=read(BASE/'analysis/RUN_ANALYSIS.json');m=r['medium_train__PAIR__2000us__DT1'];e=r['medium_train__PAIR__5000us__DT1'];c=read(BASE/'analysis/COSTS.json')
 paragraph=f"""
{marker}这一轮重新把问题放回“散开的泡能否靠近”。我们保留 F10S 的完整逐泡求解器、APECSS、均匀水圆片声场、F2、润滑及 Sequential QP 原门槛，新增 8 个预登记初态、13 种有限动作时序。32 泡、半径 1.78 µm、2.5 MHz、背景流严格为零，幅值仍是 Pa-equivalent，不是电压；原始和近接触演示初态明确不算声学形成的团。充水均匀介质不等于已经校准真实 PDMS 装置，原平移阻力黏度 0.001 Pa·s 与径向有效黏度 0.002 Pa·s 分别原样保留。

真实完成 {c['completed_trajectories']} 条积分，没有积分失败；其中标准 2 ms 窗口包含 37 个初态—动作组合，步长和缓存复算不是新的独立样本。OFF 完全静止；中等分散态 X 双源 2 ms 的 COM X 位移 {m['COM_displacement_m'][0]*1e6:.6f} µm，三维 Rg 下降 {(m['Rg_initial_final_m'][0]-m['Rg_initial_final_m'][1])*1e6:.6f} µm，对照差异约 {m['Rg_numerical_difference_m']*1e9:.6f} nm。去质心形变显著超过数值差异，几何分量从 32 降至 28，最大局部组为 3 泡；这是局部邻近，不是全部 32 泡形成稳定团，也不是融合。初始间距与布局不同，收缩和新连边明显不同；原始连通团仍主要平移。

中等态独立完整积分扩展到 5 ms，COM X 位移 {e['COM_displacement_m'][0]*1e6:.6f} µm，Rg 收缩 {-e['relative_Rg_change']*100:.4f}%。25 组完整窗口的 1/0.5 µs 步长和 3/1.5 µm 缓存对照通过连续量门槛，墙修正为零。下一档 10 ms 因预登记前级 5% 形变范围保护和三条复算成本条件未准入；这不是已知物理失效阈值，也不是 10 ms 计算失败。1 秒尺度没有完整求解依据。

幅值和相位变化可以改变收缩与位移，但尚未证明独立 XY 位置—形状目标控制；Y 行为仍限 HOMOGENEOUS_WATER_PROXY_SPECIFIC。两种独立中等初态用于 0.5 ms 观测后预测 1.5 ms 的状态保持、恒速和动作均值基准；没有随机拆连续帧，没有未见动作×未见初态的联合模型验证，所以 WORLD_MODEL_DATA_NOT_READY，未训练神经世界模型，MPPI_CONTROL_NOT_RUN。统一 XY 中心、Rg、Q、NN、局部密度和汇聚指标可与 F10Y 做趋势对照，但投影真值不冒充显微相机输出。

远端三批实际墙钟 {c['remote_campaign_wall_s']:.3f} s；GPU 与 CPU 混合使用，A 批监测线程缺依赖的异常保留，B/C 已用标准库修复，未放宽物理门槛。8 组 PNG/PDF、3 段逐帧可追溯 MP4、35 问中文报告、原始/派生 HDF5、SHA 与失败/未运行依据保存在 {BASE}。保持 EXPLORATORY、PHYSICS_UNCALIBRATED、ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED、training_ready=false。下一轮优先做独立动作设计、聚集适用的长窗准入与真实 OFF/ON 趋势核对，之后再决定联合世界模型。
"""
 text=old.decode();anchor='### 尾声：故事里还没讲完的部分';assert text.count(anchor)==1
 after=text.replace(anchor,paragraph+'\n'+anchor).encode()
 heads=lambda b:re.findall(r'^#{1,6} .*$',b.decode(),re.M)
 assert heads(old)==heads(after)
 path.write_bytes(after);(BASE/'evidence/progress_after.md').write_bytes(after)
 diff=''.join(difflib.unified_diff(text.splitlines(True),after.decode().splitlines(True),fromfile=str(path)+':before',tofile=str(path)+':after'))
 (BASE/'evidence/progress.diff').write_text(diff)
 write(BASE/'evidence/PROGRESS_UPDATE.json',dict(path=str(path),backup=str(backup),before_sha256=sha(before),after_sha256=sha(path),backup_sha256=sha(backup),heading_structure_unchanged=True,before_headings=heads(old),after_headings=heads(after)))
 print('Progress safely updated; prior text/headings retained.')
def manifest():
 rows=[]
 for p in sorted(BASE.rglob('*')):
  if p.is_file() and p.name!='INPUT_OUTPUT_SHA256.json':
   rows.append(dict(path=str(p.relative_to(BASE)),bytes=p.stat().st_size,sha256=sha(p)))
 inputs=[dict(path=rel,sha256=sha(ROOT/rel)) for rel in read(BASE/'PREREGISTRATION.json')['input_sha256']]
 source=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for folder in ['scripts/straight_tube_exp1/f10z','tests/f10z','configs/straight_tube_exp1/f10z'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p)]
 write(BASE/'INPUT_OUTPUT_SHA256.json',dict(identity=identity(),files=rows,input_files=inputs,source_files=source,manifest_self_excluded=True,large_files='Kept locally and on server; not implicitly added to Git or pushed'))
 print('Hashed',len(rows),'artifacts.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['audit','progress','manifest']);a=p.parse_args()
 {'audit':audit,'progress':progress,'manifest':manifest}[a.mode]()

