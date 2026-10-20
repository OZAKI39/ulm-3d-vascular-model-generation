"""WSL-only, byte-preserving research record update with hash/readback evidence."""
from pathlib import Path
import json,hashlib,difflib,os
R=Path(__file__).resolve().parents[2];TASK='clot_biphasic_heterogeneous_shear_fragmentation_v2';RUN='20261017T231900Z';O=R/'verification'/TASK/RUN
def read(p):return json.loads(p.read_text())
def sha(b):return hashlib.sha256(b).hexdigest()
def update(path,block,folder,name,insert=None):
 before=path.read_bytes();old=before.decode();assert block.splitlines()[0] not in old
 folder.mkdir(exist_ok=False);(folder/(path.stem+'.before.md')).write_bytes(before)
 i=old.index(insert) if insert else 0;new=old[:i]+block+'\n\n'+old[i:]
 assert sha(path.read_bytes())==sha(before)
 tmp=path.with_suffix('.biphasic.tmp');tmp.write_text(new);os.replace(tmp,path);after=path.read_bytes()
 assert after.decode()==new and new.replace(block+'\n\n','',1)==old
 (folder/(path.stem+'.after.md')).write_bytes(after);(folder/(path.stem+'.diff')).write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=path.name+'.before',tofile=path.name+'.after')))
 proof=dict(status='PASS',path=str(path) if name=='WSL_DOCUMENT_UPDATE' else str(path.relative_to(R)),before_sha256=sha(before),after_sha256=sha(after),history_preserved_exactly=True,readback_equal=True)
 (O/(name+'.json')).write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n');print(name,proof)
assert read(O/'DELIVERY_VALIDATION.json')['video_full_decode']=='PASS'
status=read(O/'FINAL_STATUS.json');data=read(O/'COMPARISON.json')['cases'];p0=read(O/'PREFLIGHT_STRUCTURE_AUDIT.json');path=read(O/'PRESCRIBED_PATH_ENERGY.json')
lines=['## Biphasic Heterogeneous Shear Fragmentation V2',f'**明确强弱分区、新大试样与同一均匀剪切对照（{RUN}）**',
 '研究问题：在完全相同的流体剪切类机械驱动下，多个强区嵌入连通弱基质能否自然形成至少5个20–200原始材料节点、内部仍有存活材料连接的游离块。区域标签不是碎块定义。',
 f"材料中心包络1.2×1.2×0.6mm，dx=50µm、horizon=150.5µm；8125材料点+1875固定锚点（3层），原始bulk/interface/scaffold键={p0['original_bulk_bonds']}/{p0['original_interface_bonds']}/{p0['anchor_scaffold_bonds']}。共同单边缺口{p0['initial_precrack']}键；奇数z网格把z=0点归下侧，跨z=0且交点x≤−50µm的bulk键预切，无相区相关切割。初始整体仍连通。",
 f"seed=20261018；受约束面连通区域生长得到10个强区，每区270–271点，强相体积分数{p0['strong_volume_fraction']:.8%}，弱基质面邻接分量1。结构中心跨度250–350µm，等体积直径约401µm；相区最小节点间距100µm<horizon，非局部连接可跨弱隔层。原始跨相键{p0['cross_phase_material_bonds']}。P0结构图与数量先行展示，通过后冻结。",
 'WEAK s0/sc=.100/.120、STRONG=.250/.300、MIXED端点算术均值=.175/.210；起裂、完成及区间比均为2.5，无平均Gc归一化。U=.15/.18。G=940Pa、K=100000/3Pa、rho=1000kg/m³及interface=.140625/.16875均沿用。参考端节点仍为V=dx³，总节点体积1.015625mm³，与中心包络0.864mm³不同；未重标质量。',
 '完整endpoint-star变分储能与力核不改，原生新增原始ID双阈值原子写入接口，避免分别设s0/sc时无效中间卡；reference/restart读取上限由4096节点放宽。旧实现对照库仅改两处容量限制，均匀旧/新力、能量、历史逐位一致。代表态梯度、节点排序、不可逆历史和V4 restart通过。U/H预切共同κ=.30、d=1，其余κ/d=0，保证完整初态一致；预切κ是失效标记，不是物理初应变。',
 f"四类完整模型规定开裂路径耗能WEAK/MIXED/STRONG/UNIFORM分别为"+'/'.join(f"{path['summaries'][m]['path_dissipation_J']*1e9:.6f}" for m in ['WEAK','MIXED','STRONG','UNIFORM'])+f"nJ，STRONG/WEAK={path['strong_weak_path_dissipation_ratio']:.6f}；路径梯形积分残差单列，非Gc标定。",
 '加载UNIFORM_SURFACE_SHEAR_X，名义1000Pa沿+X投影到当前切面，不加空间窗口、不使用等功反馈。S(u)=3u²−2u³；0–2ms升、2–2.4ms恒、2.4–2.8ms降、2.8–3.6ms保留被动拖曳。当前受支持表面、裂面可达性及接触不改。外部速度上界25m/s、Re_h=625，仅为机械剪切需求代理，未验证微泡声流功率。']
for m,q in data.items():
 s=q['native'];o=q['objective'];w=q['work']
 lines.append(f"{m}：{s['status']}，{s['time_s']*1e3:.7f}ms，{s['accepted_steps']}接受步、{s['predictive_rejections']}拒绝预测；目标块末态{o['n_fragments_20_200']}、历史最大{s['max_medium_components']}，目标质量占初始{o['mass_fraction_target_initial']:.7%}；总游离质量占初始{o['total_detached_mass_fraction']:.7%}、仍附着{o['attached_body_fraction']:.7%}；最大自由块{max(q['final_free_sizes'],default=0)}点；细屑/脱落质量{o['fine_debris_mass_fraction_detached']}；驱动关闭至终时同成员目标块存活{o['drive_off_target_survivors_to_end']}个。")
 lines.append(f"{m}：驱动/拖曳/净功={w['drive_J']*1e9:.7f}/{w['drag_J']*1e9:.7f}/{w['actual_net_J']*1e9:.7f}nJ；最大相对数值余额{s['max_relative_energy_residual']:.7%}，R不并入D。CPU={s['cpu_s']:.3f}s、wall={s['wall_s']:.3f}s、peakRSS={s['peak_RSS_bytes']/2**20:.3f}MiB、FP64状态{s['saved_state_count']}个。")
if (O/'USER_RESUME_STATE.json').exists():
 resume=read(O/'USER_RESUME_STATE.json');paused=resume['pause_record_to_resume_elapsed_s']
 lines.append(f"H因用户临时关闭本地电脑而SIGSTOP暂停，后经明确指令恢复同一进程，未新增正式尝试。原始wall包括暂停；记录至恢复间隔{paused:.3f}s，扣除此间隔后约{data['biphasic']['native']['wall_s']-paused:.3f}s。记录略晚于实际暂停，因此扣除值仅为近似，CPU不受等待时长影响。")
lines.append('末态目标块加载表示：'+json.dumps(status['TARGET_LOADING_REPRESENTATION'],ensure_ascii=False)+'。存活图上的凝聚性与完整流体作用下的稳定性必须区分。')
lines += [f"科研结论：{status['SCIENTIFIC_CONCLUSION']}。原始存活材料键图与锚层可达性定义自由块，统计1–2、3–13、14–19、20–33、34–113、114–200、>200；20/200包含。通过原始成员集合追踪birth、二次分裂、寿命及卸载存活。空间多处软化和failed-bond网络不能单独认证连续三维裂面。",f"下一项物理优先级：{status['NEXT_PHYSICAL_PRIORITY']}。",'本轮只有一个冻结空间实现，没有增加种子/幅值扫描，没有内存或1800s墙钟预算。细屑和缺少3D受支持表面的分量水动力未解析；不能解释为物理零阻力。材料/界面、接触、整体强度和Gc仍未标定。',f'详见[科研报告](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{TASK}/{RUN}/REVIEW_CLOT_BIPHASIC_HETEROGENEOUS_SHEAR_FRAGMENTATION_V2.md)、[科研图和两个1080p真实视频](temp_storage/clot_github_sync_20261015/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)。视频固定全景、实际坐标位移×1，每面板标真实接受态时间；未观测目标块时显示TARGET_FRAGMENT_NOT_OBSERVED。']
update(Path('/home/lzy/projects/clot_rupture.md'),'\n\n'.join(lines),O/'wsl_document','WSL_DOCUMENT_UPDATE',insert='# ⑥ 外加剪切')
entry=[f'# 当前双相强弱区域均匀剪切入口（{RUN}）',f'当前任务 {TASK}；analysis/{TASK}/；8125材料点+1875固定锚点，10强区、33.329%体积分数，双阈值2.5倍差异。',f"状态={status['SCIENTIFIC_CONCLUSION']}；U/H={status['UNIFORM_NATIVE_STATUS']}/{status['BIPHASIC_NATIVE_STATUS']}；目标20–200块={status['TARGET_FRAGMENT_COUNT_20_200']}。",'EXPLORATORY_UNCALIBRATED_BIPHASIC_MATERIAL；PRESCRIBED_SHEAR_MECHANICAL_SURROGATE。1000Pa均匀+X投影，指定三次smoothstep；没有微泡声流功率验证。',f'先读verification/{TASK}/{RUN}/FINAL_STATUS.json、REVIEW_CLOT_BIPHASIC_HETEROGENEOUS_SHEAR_FRAGMENTATION_V2.md、DELIVERY_VALIDATION.json；图及两个MP4位于results/{TASK}/{RUN}/OPEN_RESULTS.html。','本轮两组正式尝试均已使用，不追加重跑；新物理试验建立新RUN并保留所有原生接受数据与哈希。','## 历史入口原文（保留）']
update(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md','\n\n'.join(entry),O/'entrypoints','ENTRYPOINT_UPDATE')
