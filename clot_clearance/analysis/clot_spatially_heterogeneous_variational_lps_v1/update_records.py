"""Stdlib-only WSL update, preserving every byte outside the added block."""
from pathlib import Path
import json,hashlib,difflib,os
R=Path(__file__).resolve().parents[2];TASK='clot_spatially_heterogeneous_variational_lps_v1';RUN='20261017T214917Z';O=R/'verification'/TASK/RUN
def read(p):return json.loads(p.read_text())
def sha(b):return hashlib.sha256(b).hexdigest()
def update(path,block,folder,name,insert=None):
    before=path.read_bytes();old=before.decode();assert block.splitlines()[0] not in old
    folder.mkdir(exist_ok=False);backup=folder/(path.stem+'.before.md');backup.write_bytes(before)
    i=old.index(insert) if insert else 0;new=old[:i]+block+'\n\n'+old[i:]
    assert sha(path.read_bytes())==sha(before)
    tmp=path.with_suffix('.v2.tmp');tmp.write_text(new);os.replace(tmp,path)
    after=path.read_bytes();assert after.decode()==new and new.replace(block+'\n\n','',1)==old
    (folder/(path.stem+'.after.md')).write_bytes(after)
    (folder/(path.stem+'.diff')).write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=path.name+'.before',tofile=path.name+'.after')))
    proof=dict(status='PASS',path=str(path) if name=='WSL_DOCUMENT_UPDATE' else str(path.relative_to(R)),before_sha256=sha(before),after_sha256=sha(after),history_preserved_exactly=True,readback_equal=True)
    (O/(name+'.json')).write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n')
    print(name,proof)
assert read(O/'DELIVERY_VALIDATION.json')['status']=='PASS'
status=read(O/'FINAL_STATUS.json');data=read(O/'COMPARISON.json')['cases'];meta=read(O/'HETEROGENEITY_METADATA.json')
lines=['## Spatially Heterogeneous Variational LPS under Shear',
 f'**三维空间相关软化与同一开环剪切对照（{RUN}）**',
 '研究假设：在不预制碎块、不删除原始非局部键的条件下，空间相关的后峰软化差异是否改变凝聚性多块碎裂。保留1792材料点+336锚点、95666原始键及2736预切；X/m/r0/taper/质量、G=940Pa、K=100000/3Pa、bulk s0=.15及interface .140625/.16875均不改。',
 'Fan 2022异质LPS提供随机场与对称两点参数方法背景；Staessens组织学和三维断层研究支持组织空间异质性，JNIS分析使用的0.2–0.6mm网格不等于材料相关长度。Eyisoylu 2026支持局部组成相关刚度，但本轮G/K不变；Liu 2024强调bulk凝聚与界面黏附区分，不能由刚度或sc直接推断Gc。',
 f"GSTools {meta['gstools_version']}，Gaussian C_z(r)=exp[-(r/300µm)²]、rescale=1、seed=20261017、4096模式。q=tanh(z−c)，材料体积加权均值0；q实际范围{meta['q_min']:.6f}–{meta['q_max']:.6f}。qij=(qi+qj)/2，Δsij=.03exp(.30qij)/Z，scij=.15+Δsij；Z采用全部原始bulk键参考h*r0*V²*(1/mi+1/mj)权重归一化。sc范围{meta['sc_min']:.9f}–{meta['sc_max']:.9f}，加权均值{meta['sc_weighted_mean']:.9f}。", 
 f"lc/dx=6，lc/δ={meta['lc_over_horizon']:.6f}；线性梯度R²={meta['linear_gradient_R2']:.6f}，正负q各1个连通域。多个颜色和边界极值不是多个统计独立组织块；300µm和alpha=.30均是未标定探索值。", 
 '原生新增按原始bond ID冻结的sc字段及restart V4。完整endpoint-star储能与−∂U/∂x公式不改，sc只进入原生κ→d更新；D按固定x的完整储能下降累加，R独立。q=0与旧native逐位回退；异质有限差分、不可逆性、排序及restart通过。原始预切d=1保持，但继承初始化采用κ=sc作为预切失效标记，H预切κ数值随sc变化；不能声称完整初始κ数组逐值相同。未预切κ/d初始不变，预切材料力/储能不变。',
 '同一250µm棋盘涡旋模板、当前受支持表面与流体可达代理、被动滑移拖曳。冻结上一轮checkerboard接受幅值记录为开环函数，取消等功反馈；2ms升、2.4ms降、2.8ms停驱，继续观察到3.6ms。允许两组实际功不同。数千Pa为机械需求代理，真实微泡声流可实现性仍未知。']
for m,q in data.items():
 s=q['native'];o=q['objective'];w=q['work'];f=q['surface']['final']
 lines.append(f"{m}：{s['status']}，{s['accepted_steps']}接受步，净功{w['actual_net_J']*1e9:.7f}nJ（驱动{w['drive_J']*1e9:.7f}，拖曳{w['drag_J']*1e9:.7f}）。峰值实际剪切{w['actual_peak_shear_Pa']:.5f}Pa。20–200点块末态{o['n_fragments_20_200']}、历史最多{s['max_medium_components']}；最大自由块{max(q['final_free_sizes'],default=0)}点；细屑/脱落质量{o['fine_debris_mass_fraction_detached']}；自然bulk失效{s['last_row']['bulk_failed']}。")
 lines.append(f"{m}：最大能量残差{s['max_relative_energy_residual']*100:.6f}%，CPU{s['cpu_s']:.3f}s、wall{s['wall_s']:.3f}s、峰值RSS{s['peak_RSS_bytes']/2**20:.3f}MiB。末态未表示3D湿润表面的自由分量{f['free_components_without_surface']}/{f['true_free_components']}；细屑流体受力未解析，不能等同物理零力。")
lines += [f"科研结论：{status['SCIENTIFIC_CONCLUSION']}。只有一个空间实现；真实图连通性、内部存活键、birth/寿命及卸载观察判定块体，不以损伤区域数或接触邻近代替。相关性不推出因果；不把离线规定路径耗能称GC_EFFECTIVE。", 
 '图1–8及支持图来自冻结输入、明确标记的材料验证/离线路径和真实接受态。三个1080p H.264视频展示同驱动对照、真实失效/存活键、真实分量形成与二次分裂；位移×1，使用实际保存状态，各面板标实际时间，事件采样允许物理时间间隔不同但不插值。目标块不存在时保留TARGET_FRAGMENT_NOT_OBSERVED。',
 '本轮没有RSS/1800s预算，没有追加随机种子、材料调参、BoTorch、微泡CFD或完整V3。下一轮优先建立空间相关尺度与当前horizon/有限试样的可分辨性证据；若观察到完整目标，则先检验独立空间实现的重复性。',
 f'详见[28项科研报告](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{TASK}/{RUN}/REVIEW_CLOT_SPATIALLY_HETEROGENEOUS_VARIATIONAL_LPS_V1.md)、[完整图与三个视频](temp_storage/clot_github_sync_20261015/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)。']
update(Path('/home/lzy/projects/clot_rupture.md'),'\n\n'.join(lines),O/'wsl_document','WSL_DOCUMENT_UPDATE',insert='# ⑥ 外加剪切')
entry=[f'# 当前空间相关变分软化剪切入口（{RUN}）',f'当前任务 {TASK}；analysis/{TASK}/；新隔离native库，原始ID bulk sc字段与restart V4。',
 f"状态={status['SCIENTIFIC_CONCLUSION']}；原生U/H={status['UNIFORM_NATIVE_STATUS']}/{status['HETEROGENEOUS_NATIVE_STATUS']}；目标20–200块={status['TARGET_FRAGMENT_20_200_COUNT']}。", 
 'SPATIALLY_HETEROGENEOUS_EXPLORATORY_UNCALIBRATED。一个固定300µm材料场，G/K/s0/界面及原图不改；同一冻结棋盘涡旋开环输入。完整湿润表面及细屑水动力未解析。预切κ=sc标记的初始化数值差异在INITIAL_HISTORY_AUDIT单列。',
 f'先读verification/{TASK}/{RUN}/FINAL_STATUS.json、REVIEW_CLOT_SPATIALLY_HETEROGENEOUS_VARIATIONAL_LPS_V1.md、DELIVERY_VALIDATION.json；图/三个MP4位于results/{TASK}/{RUN}/OPEN_RESULTS.html。',
 '本轮两个正式工况已完成；新实验必须使用新RUN，保留接受数据和冻结哈希。没有自动追加调参或完整V3。','## 历史入口原文（保留）']
update(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md','\n\n'.join(entry),O/'entrypoints','ENTRYPOINT_UPDATE')
