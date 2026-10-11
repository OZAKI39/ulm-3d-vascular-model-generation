"""Read-only provenance closure and portable review indexes; no numerical solves."""
from analysis.dual_track_v1.common import *
import ast,datetime,html,importlib.metadata as md,os,subprocess,shutil

def read(p):return json.loads(p.read_text())
def command(args):
 r=subprocess.run(args,cwd=ROOT,text=True,capture_output=True)
 return dict(returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
def href(p,base):return html.escape(os.path.relpath(p,base),quote=True)
def page(title,body):
 return '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>body{max-width:1200px;margin:32px auto;padding:0 22px;color:#253440;font:17px/1.6 system-ui,sans-serif;background:#fafbfd}h1,h2{line-height:1.3}a{color:#176182}img{width:100%;height:auto;background:white}figure{margin:28px 0;padding:16px;background:white;border:1px solid #ccd4dc}figcaption{padding:10px 0}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #cad2d9;padding:9px;text-align:left;overflow-wrap:anywhere}code{overflow-wrap:anywhere}.note{padding:16px;background:#eef3f5}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style><h1>'+html.escape(title)+'</h1>'+body+'</html>'
def table(states):return '<table><tr><th>分项</th><th>状态</th></tr>'+''.join('<tr><td>'+html.escape(k)+'</td><td>'+html.escape(str(v))+'</td></tr>' for k,v in states.items())+'</table>'

def protection():
 before=read(COMMON/'provenance/PROTECTED_START.json');checks=[]
 for name,expected in before.items():
  p=Path(name);p=p if p.is_absolute() else ROOT/p;actual=sha(p) if p.is_file() else None
  checks.append(dict(path=name,before_sha256=expected,after_sha256=actual,unchanged=actual==expected))
 failure=[x for x in checks if not x['unchanged']]
 save(COMMON/'PROTECTION_AUDIT.json',dict(status='PASS' if not failure else 'FAIL',checked_files=len(checks),changed_or_missing=failure,checks=checks,scope='Exact rehash of all pre-registered production/material/native libraries and historical accepted files; no restoration/overwrite attempted'))
 assert not failure,failure
 packages={d.metadata['Name']:d.version for d in md.distributions()};original=read(COMMON/'provenance/ENVIRONMENT_BEFORE.json')['packages'];assert packages==original
 start=read(COMMON/'provenance/SERVER_START.json');head=command(['git','rev-parse','HEAD']);assert head['stdout'].strip()==start['HEAD'].strip()
 capture=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),HEAD=head,branch=command(['git','branch','--show-current']),git_status=command(['git','status','--short']),git_tracked_diff=command(['git','diff','--stat']),packages=packages,packages_unchanged=True,python=sys.executable,python_version=sys.version,compiler=command(['g++','--version']),mpi=command(['mpirun','--version']),supervisor=command(['supervisorctl','status']),processes=command(['ps','-eo','pid,ppid,pcpu,pmem,comm']),gpu=command(['nvidia-smi','--query-gpu=name,utilization.gpu,memory.used,memory.total','--format=csv,noheader']),affinity=sorted(os.sched_getaffinity(0)),CPU_quota_us=Path('/sys/fs/cgroup/cpu/cpu.cfs_quota_us').read_text().strip(),CPU_period_us=Path('/sys/fs/cgroup/cpu/cpu.cfs_period_us').read_text().strip(),meminfo=Path('/proc/meminfo').read_text(),disk=shutil.disk_usage(ROOT)._asdict())
 save(COMMON/'provenance/ENVIRONMENT_AFTER.json',capture)
 source=[];patch=[]
 for p in sorted((ROOT/'analysis/dual_track_v1').iterdir()):
  if p.suffix not in ['.py','.cpp','.txt']:continue
  if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
  source.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size));r=command(['git','diff','--no-index','--','/dev/null',str(p.relative_to(ROOT))]);assert r['returncode'] in [0,1];patch.append(r['stdout'])
 save(COMMON/'SOURCE_MANIFEST.json',dict(files=source,python_AST='PASS',new_dependencies_installed=False,isolated_compilers_only=True))
 (COMMON/'CODE_ADDITIONS.diff').write_text(''.join(patch));shutil.copy2(ROOT/'analysis/dual_track_v1/REPRODUCE.txt',COMMON/'REPRODUCE.txt')

def indexes():
 a=read(A/'ANALYSIS_SUMMARY.json');b=read(B/'ANALYSIS_SUMMARY.json');t=b['traction'][b['selected_level']];metrics=b['flows'][b['selected_level']]['metrics'];nxt=read(B/'NEXT_PRIORITY.json')['recommendation'];doc=COMMON/'document_edit/AUDIT.json';docstate='PASS_WRITTEN_AND_READ_BACK' if doc.exists() and read(doc)['status']=='PASS' else 'PENDING_LOCAL_SERIAL_UPDATE';visual=COMMON/'VISUAL_REVIEW.json';visualstate='PASS_14_REAL_FIGURES_REVIEWED' if visual.exists() and read(visual)['status']=='PASS' else 'PENDING_FINAL_VISUAL_REVIEW';protection_audit=read(COMMON/'PROTECTION_AUDIT.json')
 cost=dict(run_id=RUN,formal_A_equilibrium_attempts=5,formal_A_budget=8,formal_B_CFD_attempts=2,formal_B_budget=2,A_known_static_wall_s=a['measured_low_high_wall_s'],A_recovery_wall_upper_bound_s=a['recovery_wall_upper_bound_s'],A_known_CPU_s=a['measured_low_high_CPU_s'],A_recovery_CPU_RSS_unknown=True,A_peak_recorded_RSS_MiB=a['peak_observed_RSS_MiB'],B_native_runtime_s=b['total_accepted_native_runtime_s'],B_telemetry=b['telemetry'],known_formal_wall_sum_excluding_recovery_s=a['measured_low_high_wall_s']+b['total_accepted_native_runtime_s'],time_accounting='Sum of formal solve wall times is not end-to-end duration. A recovery wall is an upper bound, CPU/RSS missing. B CPU is sampled live-process lower bound; GPU samples are device-wide.',report_generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_since_run_id_s=(datetime.datetime.now(datetime.timezone.utc)-datetime.datetime.strptime(RUN,'%Y%m%dT%H%M%SZ').replace(tzinfo=datetime.timezone.utc)).total_seconds())
 save(COMMON/'RESOURCE_AND_COST_SUMMARY.json',cost)
 states=dict(PRODUCTION_MATERIAL_UNCHANGED=protection_audit['status'],NATIVE_LAMMPS_UNCHANGED=protection_audit['status'],HISTORICAL_RESULTS_PROTECTED=protection_audit['status'],WSL_DOCUMENT_UPDATE=docstate,VISUALIZATION_COMPLETE=visualstate)
 save(COMMON/'FINAL_STATUS.json',dict(run_id=RUN,common_states=states,A_status_file=str((A/'FINAL_STATUS.json').relative_to(ROOT)),B_status_file=str((B/'FINAL_STATUS.json').relative_to(ROOT)),no_global_scientific_PASS=True,task_numerical_work_stopped=True,next_round_automatically_started=False))
 # Machine-readable figure identity after final rendering.
 for verification,figdir in [(A,FA),(B,FB)]:
  figs=sorted(figdir.glob('Fig*.png'));assert len(figs)==7
  save(verification/'FINAL_FIGURE_MANIFEST.json',dict(figures=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p),svg_sha256=sha(p.with_suffix('.svg'))) for p in figs],field_source='refined accepted native CFD' if verification==B else 'actual saved topology/equilibria and unaccepted high-candidate diagnostic',deformation_scale=1,no_AI_or_medical_shell=True))
 for verification,figdir,title,lead in [(A,FA,'A｜宽薄周期试样与 PVCE','低载三态通过；旧高载恢复通过；宽薄高载候选未验收，已停止。G 是试样驱动力，不是材料 Gc。'),(B,FB,'B｜原生 Y 血管中的 V3 含栓 CFD','基线与一次局部加密均通过；牵引已提取，15% 工程网格筛查通过。刚性背景流，不含 streaming 或双向 FSI。')]:
  review=next(verification.glob('REVIEW_*.md'));stat=read(verification/'FINAL_STATUS.json')['states']
  for target in [verification/'OPEN_RESULTS.html',figdir/'OPEN_RESULTS.html']:
   base=target.parent;body='<p class="note">'+lead+'</p><p>run_id：<code>'+RUN+'</code></p><p><a href="'+href(COMMON/'OPEN_RESULTS.html',base)+'">总索引</a> · <a href="'+href(review,base)+'">中文科研报告（12 项回答）</a> · <a href="'+href(verification/'FINAL_STATUS.json',base)+'">状态 JSON</a> · <a href="'+href(verification/'DESIGN_FREEZE.json',base)+'">冻结方案</a></p>'+table(stat)
   for p in sorted(figdir.glob('Fig*.png')):body+='<figure><a href="'+href(p,base)+'"><img loading="lazy" src="'+href(p,base)+'" alt="'+html.escape(p.stem)+'"></a><figcaption>'+html.escape(p.stem)+' · <a href="'+href(p.with_suffix('.svg'),base)+'">SVG</a></figcaption></figure>'
   body+='<p>完整原始结果留在服务器：<code>'+str(verification)+'</code>；本地回收范围见总索引的 COMPACT_MANIFEST.json。</p>';target.write_text(page(title,body))
 report=f'''# Y_CLOT_CFD_AND_WIDE_PVCE_DUAL_TRACK_V1 共同审查

run_id：`{RUN}`。A、B 独立验收，没有单一总科学 PASS。A 的 12 项回答见相邻 A 报告，B 的 12 项回答见相邻 B 报告；下列覆盖共同问题 13–19。

13. **正式算例和资源。** A 共 5/8：低载三态、旧保存态恢复一次、宽薄高载候选一次；B 共 2/2：含栓基线和一次局部体网格加密。没有新动态断裂、无栓重算或第三个网格。A 已记录正式墙钟 {a['measured_low_high_wall_s']:.3f} s、CPU {a['measured_low_high_CPU_s']:.3f} s，另有恢复墙钟上界 {a['recovery_wall_upper_bound_s']:.3f} s；恢复 CPU/RSS/迭代遥测因报告序列化异常不可恢复，不能以零代替。B 两次原生墙钟合计 {b['total_accepted_native_runtime_s']:.3f} s（{b['total_accepted_native_runtime_s']/60:.3f} min）。已记录正式墙钟相加 {cost['known_formal_wall_sum_excluding_recovery_s']:.3f} s，另加上述恢复上界；不把并行任务墙钟相加冒充端到端时长。完整时长、CPU 下界、内存/GPU 样本与当前报告时刻见 `RESOURCE_AND_COST_SUMMARY.json`。A 所记录峰值 RSS {a['peak_observed_RSS_MiB']:.2f} MiB；B 峰值采样 RSS 分别 {b['telemetry']['baseline']['process_tree_peak_sampled_RSS_MiB']:.2f}/{b['telemetry']['refined']['process_tree_peak_sampled_RSS_MiB']:.2f} MiB。15.36 CPU 配额下采用 A 单进程、B 8 MPI，优先 B。4090 沿用原生 PETSc GPU 路径；无压力测试或依赖升级。

14. **实际审查图。** 共 14 组 PNG＋SVG。A1 点阵/控制自由度/尺寸比，A2 周期影像邻居和原生一致性，A3 三个参考切口，A4 低载收敛/能差/G，A5 旧态恢复及宽态停止，A6 未接受高载候选真实键，A7 三个实际低载平衡。B1 原生 Y 与 V3 材料点，B2 真实体网格与开放绕流截面，B3 公共位置速度/压力，B4 实际暴露面压力/剪切，B5 守恒/稳态/压降/力/加密，B6 原生速度场流线与方向，B7 横截面和切向牵引方向。A/B 分别有索引，并有共同 `OPEN_RESULTS.html`。没有新医学外壳、AI 生成图或无分析意义动画。

15. **未完成和停止。** A 低载 G_two/G_three 为 {a['G']['G_two_J_m2']:.12e}/{a['G']['G_three_J_m2']:.12e} J/m²。旧小态高载恢复通过；宽薄高载候选的末十次位置变化及 1 pm 分量差分失败，原冻结门槛保留，没有真正首达阈值区间、高载 A1/A2、材料 Gc。只读 1 nm 差分支持舍入解释，但未重评接受。B 两个解及一次工程筛查完成；表面三角化收敛、接触线极值、到材料点载荷映射、超声 streaming、变形与双向 FSI 不在本轮完成范围。两线按范围停止。

16. **保护情况。** `PROTECTION_AUDIT.json` 对 {protection_audit['checked_files']} 个事先登记的源码、材料卡、原生库和历史结果重新计算 SHA，结果 {protection_audit['status']}。原生环境包版本完全一致，HEAD 未变。新代码只在 `analysis/dual_track_v1`，新数据只在本 run 的独立目录；没有修改 s00、生产载荷或历史接受态，没有 commit/push/reset/clean/stash，也没有停止其他任务。A 的 Y 周期独立能量/HVP 适配仅在隔离诊断目录。

17. **实际 WSL 文档。** `/home/lzy/projects/clot_rupture.md` 状态：`{docstate}`。更新保留 ①–⑨、历史无栓背景流使用记录及各轮负面科研结论，新增含栓解和宽薄 PVCE 的已验证结论；不会把旧碎裂或输运轨迹改称使用了新含栓流场。

18. **备份、diff 与 SHA。** `document_edit/before.md`、`document_edit/after.md`、`document_edit/change.diff`、`document_edit/AUDIT.json` 记录最新原文、并发检查、写入/读回和前后 SHA；任务开始备份另在 `provenance/clot_rupture.before.md`。`SOURCE_MANIFEST.json`、`CODE_ADDITIONS.diff`、`PROTECTION_AUDIT.json`、`provenance/ENVIRONMENT_BEFORE.json` 与 `ENVIRONMENT_AFTER.json`、`REPRODUCE.txt` 提供源码/库身份、代码差异和条件复现说明。`SERVER_ARTIFACT_MANIFEST.json` 清点本轮服务器全部原始数据；`COMPACT_MANIFEST.json` 与 `LOCAL_RETRIEVAL_AUDIT.json` 记录 WSL 精简回收及 SHA 验证。完整 CFD 网格、原生全场、逐步状态和 LAMMPS 原始输入均保留服务器，不降低输出精度。

19. **下一轮唯一优先项。** {nxt} 理由是 B 已得到含栓场及稳定的积分载荷，仍缺少到实际材料点的受力接口验证；A 未识别 Gc 不否定 B 的独立成果。建议只做保力/保力矩离线检查，不自动启动生产加载、材料升级、AI 大扫描或双向耦合。

加密 CFD 代表值：最大速度 {metrics['velocity_max_m_s']*1e3:.9f} mm/s；平均表压 {t['pressure_Pa']['mean']:.9f} Pa；平均/P95 切向剪切 {t['shear_Pa']['mean']:.9f}/{t['shear_Pa']['p95']:.9f} Pa；总力 {t['total_force_N']} N，模长 {t['total_force_norm_N']*1e6:.9f} µN。压力是出口参考表压，力只含暴露流体面，不含壁面粘附反力。
'''
 (COMMON/'REVIEW_DUAL_TRACK_V1.md').write_text(report)
 links=[(A/'OPEN_RESULTS.html','A：宽薄 PVCE 图与数据'),(B/'OPEN_RESULTS.html','B：含栓 CFD 图与数据'),(COMMON/'REVIEW_DUAL_TRACK_V1.md','共同报告：问题 13–19'),(COMMON/'FINAL_STATUS.json','共同状态'),(COMMON/'RESOURCE_AND_COST_SUMMARY.json','实际资源和成本'),(COMMON/'PROTECTION_AUDIT.json','历史与生产 SHA 保护审计'),(COMMON/'REPRODUCE.txt','复现说明'),(COMMON/'CODE_ADDITIONS.diff','本轮代码差异'),(COMMON/'SOURCE_MANIFEST.json','新增源码 SHA'),(COMMON/'document_edit/AUDIT.json','WSL 文档实际更新审计'),(COMMON/'COMPACT_MANIFEST.json','精简回收清单'),(COMMON/'LOCAL_RETRIEVAL_AUDIT.json','本地回收 SHA 核验'),(COMMON/'SERVER_ARTIFACT_MANIFEST.json','完整服务器数据清单')]
 body='<p>run_id：<code>'+RUN+'</code></p><div class="note"><p><b>A：</b>16,000 点；低载三态通过，G ≈ 7.14065×10⁻⁶ J/m²。高载候选未验收并停止；材料 Gc 未识别。</p><p><b>B：</b>原生 Y＋原 V3，两个含栓稳态解通过，15% 工程网格筛查通过。平均剪切约 0.01677 Pa，总力约 0.6638 µN。刚性背景流，不含 streaming/FSI。</p></div><ul>'+''.join('<li><a href="'+href(p,COMMON)+'">'+label+'</a></li>' for p,label in links)+'</ul>'+table(states)+'<p>本轮计算已停止，不自动执行下一轮。</p>'
 (COMMON/'OPEN_RESULTS.html').write_text(page('Y 血栓｜含栓 CFD 与宽薄 PVCE 双线审查',body))
if __name__=='__main__':
 import sys
 if sys.argv[1]=='protection':protection()
 elif sys.argv[1]=='indexes':indexes()
