from common import *
import csv,re,shutil
from html.parser import HTMLParser
from urllib.parse import unquote
def main():
 assert not (O/'PAYLOAD_MANIFEST.json').exists(),'No resealing frozen results'
 # Clarify the exact scope of recovered diagnostics; never promote keyframe
 # extrema to complete online extrema after an interrupted finalizer.
 p=O/'DYNAMIC_ENERGY_LEDGER.json';j=json.loads(p.read_text());j['independent_saved_keyframe_errors']=j.pop('independent_per_accepted_step_errors');j['online_all_accepted_steps_passed_frozen_gates']=True;j['exact_online_extrema']='NOT_SAVED_AFTER_NATIVE_EXCEPTION';dump(p,j)
 p=O/'CONTACT_AND_TIMELAYER_AUDIT.json';j=json.loads(p.read_text());j['error_scope']='Contact force/energy and loading extrema are saved full keyframes; VV/momentum extrema are every accepted scalar row. All accepted steps passed online checks.';dump(p,j)
 p=O/'RUN_RESOURCES.json';j=json.loads(p.read_text());j['compile_calls']=2;j['compile_scope']='One initial complete build plus one incremental read-only observer build';j['main_native_candidates_including_aborted']=2760;dump(p,j)
 # Validate scientific records without another dynamics case.
 rows=list(csv.DictReader((O/'full/step_scalars.csv').open()));a=np.array([[float(r[k]) for k in rows[0]] for r in rows]);cols=list(rows[0]);t={c:a[:,i] for i,c in enumerate(cols)}
 assert len(rows)==2756 and np.all(np.diff(t['time_s'])>0)
 assert np.max(abs(np.diff(t['time_s'])-t['dt_s'][1:]))<1e-15
 assert t['max_dd'][1:].max()<=.10 and t['relative_R'].max()<=.02
 assert t['time_s'][-1]<.0007 and abs(t['time_s'][-1]-.00068775)<1e-15
 assert np.all(t['commits'][1:]==t['step'][1:]) and np.all(t['token'][1:]==t['step'][1:])
 assert np.all(np.diff(t['D_J'])>=-1e-20) and np.all(np.diff(t['total_failed'])>=0)
 for case in ['local','full']:
  tx=json.loads((O/case/'TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());rb=json.loads((O/case/'ROLLBACK_STATE_CONSISTENCY.json').read_text());assert tx['rejected']==4 and len(rb['records'])==4
  assert all(x['pass_'] and x['contact_exclusions_exact'] and x['topology_exact'] for x in rb['records'])
  assert all(abs(x['covered_s']-BASE)<1e-15 for x in tx['base_intervals'])
 last=np.load(O/'full/keyframes/step_002755.npz');assert len(last['x'])==341 and np.sum(last['d']==1)==4918
 assert np.all(np.isfinite(last['x'])) and np.all(np.isfinite(last['velocity']))
 frozen=json.loads((O/'EVENT_INTEGRATOR_FREEZE.json').read_text());assert all(sha(A/n)==h for n,h in frozen['driver_sha256'].items())
 assert sha(O/'binary/liblammps.so')==frozen['library_sha256']
 doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());assert sha(O/'wsl_document/clot_rupture.after.md')==doc['readback_sha256']
 protected=json.loads((O/'PROTECTED_BASELINES.json').read_text());assert all(sha(Path(n))==h for n,h in protected.items())
 # This mutable pointer is intentionally backed up and is not an old result.
 index=f'''# 当前clot原生软化开发入口

本轮：{TASK} / {RUN}（2026-10-14 Europe/Zurich）。

- **局部损伤事件PASS，完整700µs FAILED_ATOM_LOSS。最后接受687.75µs，下一688µs试算固定域丢点。不得继续或称完整协议通过。**
- 原bond1813的0.11363469跳跃用真实VV两段0.125µs分辨，最初两段Δd=0.05736736/0.05626720；局部15接受/4拒绝，正式2755接受/4损伤拒绝，另1次未完成原生候选。
- 拒绝后x/v/F/κ/d/U/D/计数/time、实际存活键和完整special排除均精确恢复。正式35个检查点，未每步写restart；预注册最小0.03125µs，实际最小0.125µs。
- 最大接受Δd=0.0978676344，接受区间最大相对能量余额0.372426823%。全局时间收敛仍未认证。
- 卸载到650µs完成，末保持仅37.75/50µs。末态4918完全失效、7129存活；原图连通分量192/147/1/1。Gc不可用。
- 丢点诊断：原盒子fff、各轴±250µm。自由节点73最后y=249.953987866µm、vy=0.400969997m/s、F=0；下一漂移预测y=250.054230365µm。原生日志实报340vs341。预测不是接受态。
- 因原生异常后gather/restart失败，最终restart、全键稀疏事件流、精确在线误差极值和正式峰值RAM未落盘。逐步标量及18个完整FP64关键态仍在；禁止用NPZ制造restart。
- 先读verification/{TASK}/{RUN}/FINAL_STATUS.json、24项REVIEW、NATIVE_ATOM_LOSS_DIAGNOSIS.json、EVENT_INTEGRATOR_DESIGN.md。
- 活跃实现analysis/event_resolved_v1/；冻结执行代码SHA在EVENT_INTEGRATOR_FREEZE.json。run_protocol.py保留真实失败路径，禁止在封存目录重跑。
- 图表results/{TASK}/{RUN}/OPEN_RESULTS.html（五组PNG/SVG）。
- 研究库verification/{TASK}/{RUN}/build/liblammps.so，便携副本binary/liblammps.so；SHA 08ccf24110e2cf2d07dbc8c0f31c0f368368b19b5427e6e5de9cc108a5ae115d。
- 仅fix_move.cpp限定六equal规定运动变dt，新增event_state_audit.cpp只读排除表。bond/协调fix/pair/nve/Verlet源码不变；输入G/K/ρ/dx/δ/s0/sc/α、m/V/r0/taper和半余弦加载均不变。
- 动态接触最大重叠仍0.246649196dx；未作接触刚度标定，物理精度未认证。
- 生产环境source /workspace/y_clot/env.sh保持；生产库SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a。旧归档102文件在本轮入口更新前全部SHA一致，旧manifest未改；before入口保存在本轮目录。
- WSL实际文档已更新并读回，凭据WSL_DOCUMENT_UPDATE.json及wsl_document/；所有旧文本和①–⑨保留。
- 下一轮优先LPS_NATIVE_SOFTENING_DOMAIN_RETENTION_AND_FAILURE_SAFE_OUTPUT_V1，先审查自由分离点的计算域保留与异常证据持久化；本轮不自动执行。
- P4非均匀小裂纹NOT_ADMITTED/NOT_RUN；未运行V3、CFD、微泡或AI，也未部署生产或更改材料。
'''
 assert sha(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md')==sha(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
 (R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md').write_text(index)
 (O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.after.md').write_text(index)
 # All native source/binaries remain on the server; only scoped source identity,
 # actual loaded shared library and evidence are included in the local payload.
 required=['REVIEW_LPS_NATIVE_SOFTENING_EVENT_RESOLVED_INTEGRATION_V1.md','FINAL_STATUS.json','SOURCE_AND_LIBRARY_IDENTITY.json','EVENT_INTEGRATOR_DESIGN.md','TRIAL_ACCEPT_REJECT_AUDIT.json','ROLLBACK_STATE_CONSISTENCY.json','LOCAL_BOND1813_EVENT_AUDIT.json','ADAPTIVE_SUBSTEP_HISTORY.csv','FULL_PROTOCOL_DYNAMIC_AUDIT.json','DAMAGE_HISTORY_AND_TOPOLOGY.json','DYNAMIC_ENERGY_LEDGER.json','CONTACT_AND_TIMELAYER_AUDIT.json','RUN_RESOURCES.json','NATIVE_SOURCE_PATCH.diff','REPRODUCE.txt','OPEN_RESULTS.html','WSL_DOCUMENT_UPDATE.json']
 assert all((O/n).is_file() for n in required)
 assert len(list(F.glob('*.png')))==5 and len(list(F.glob('*.svg')))==5
 class Links(HTMLParser):
  def __init__(self):super().__init__();self.refs=[]
  def handle_starttag(self,tag,attrs):self.refs.extend(v for k,v in attrs if k in ['src','href'])
 checked=0
 for page in [F/'OPEN_RESULTS.html',O/'OPEN_RESULTS.html']:
  parser=Links();parser.feed(page.read_text())
  for ref in parser.refs:
   if not ref.startswith(('http','#','data:')):assert (page.parent/unquote(ref)).resolve().is_file(),ref;checked+=1
 dump(O/'DELIVERY_VALIDATION.json',dict(status='PASS_FOR_STATED_PARTIAL_SCIENTIFIC_OUTCOME',required_files=len(required),figures_png=5,figures_svg=5,html_local_links_checked=checked,accepted_rows=2755,last_accepted_time_s=float(t['time_s'][-1]),full700_pass=False,original_model_sources_unchanged=True,old_outputs_preserved=True,doc_readback_sha256=doc['readback_sha256'],missing_data_explicitly_disclosed=True,no_physical_computation_after_atom_loss=True))
 selected=[R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md']
 for folder in [A,O,F]:
  for path in folder.rglob('*'):
   rel=path.relative_to(folder)
   if not path.is_file() or '__pycache__' in rel.parts:continue
   if folder==O and (rel.parts[0] in ['source','build'] or path.name=='pre_observer_liblammps.so'):continue
   selected.append(path)
 files={str(p.relative_to(R)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(selected))}
 dump(O/'PAYLOAD_MANIFEST.json',dict(version=1,task=TASK,run_id=RUN,file_count=len(files),total_bytes=sum(v['bytes'] for v in files.values()),scope='Scoped code, modified/source-audit files, actual new binary, scientific evidence, figures, reports, WSL before/after/diff and current entrypoint. Full copied source/build caches and pre-observer library remain server-only.',files=files))
 (O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'\n')
 print('SEALED',len(files),sum(v['bytes'] for v in files.values()),sha(O/'PAYLOAD_MANIFEST.json'))
if __name__=='__main__':main()
