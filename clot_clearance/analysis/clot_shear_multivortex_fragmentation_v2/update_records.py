"""Stdlib-only WSL update, preserving every byte outside the added block."""
from pathlib import Path
import json,hashlib,difflib,os
R=Path(__file__).resolve().parents[2];TASK='clot_shear_multivortex_fragmentation_v2';RUN='20261017T175456Z';O=R/'verification'/TASK/RUN
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
status=read(O/'FINAL_STATUS.json')
data={m:read(O/p) for m,p in zip(['counter_rotating','checkerboard'],['SINGLE_VORTEX_NATIVE_RESULT.json','CHECKERBOARD_NATIVE_RESULT.json'])}
lines=['## Current Exposed Surface and Work-Matched Multivortex Shear',
    f'**当前可达裂面与等外功多涡旋剪切（{RUN}，2026-10-17）**',
    '本轮继续声流类剪切的机械需求研究，不再尝试压头：目标是让流体剪切产生多个凝聚性中等碎块，压头迭代不对应这一加载问题。科学标签 SHEAR_MECHANICAL_REQUIREMENT_ONLY，没有把声压、泡半径或质心运动替换成已求解的二阶Streaming。',
    '上一轮182个原顶面片末态只保留87个，约47.802%，未表示的新裂面与碎块表面限制了持续加载解释。本轮先加入有原始节点、存活键及实际开口依据的支撑单元边界；9828个辅助四面体仅用于几何，不承担材料刚度或质量。局部已断开且张开的间隙通过开放通道与外部连通，当前实体遮挡探针排除被覆盖表面。仅凭失效键不生成受力面。',
    '采用 QUASISTATIC_FLUID_ACCESSIBILITY_PROXY；开口1–5 µm内共同平滑门控，按实际当前位置、法向、面积与速度计算τ=(µ/h_eff)P(u_outer−v_surface)，µ=.001Pa·s、h_eff=25µm。每个唯一三角面以正1/3分配至原始ID，游离标签不另加一份力。位于驱动窗口的真实碎块继续受剪切，驱动停止后保留已表示面的拖曳。没有为1–2点细屑制造封闭表面或擅加点阻力。',
    '完整连续表面缺失比例仍不可由点数据唯一确定。闭合、倒置、未连通、缺少三维支撑的面保持未判定；只报告候选面及分量覆盖，不能声称未知水动力为真实零。面积/体积开口和有限遮挡探针不等于渗润动力学或完整流体几何认证。',
    '以历史单向的25.698197nJ实际净外功为共同预算，历史驱动27.443461nJ、拖曳−1.745264nJ。历史A采用旧表面模型，仅提供预算；新B/C才采用同一修补算法。冻结共同幅值反馈，上限6000Pa、变化率8e6Pa/s；只调全场一个幅值，不改空间中心/形状、不针对节点。驱动2.8ms关闭，观察到3.6ms，允许一次计算无法匹配并如实报告。',
    '单对沿用双Gaussian中心、125µm宽度、300µm间距。棋盘唯一新增模式尺度250µm；实际顶面650×350µm，与9个格相交，仅3个涡心在面内，中间1格完整、8格被裁切。它是 CHECKERBOARD_SHEAR_MECHANICAL_PROXY，不是已验证微泡群声流。材料、1792点+336锚点、95666原始键/2736预切均未更改。']
for m,q in data.items():
    w=q['work'];o=q['objective'];s=q['native'];sf=q['surface']
    lines.append(f'{m}：原生状态{s["status"]}，实际到{s["time_s"]*1e3:.8g}ms。净外功{w["actual_net_J"]*1e9:.9g}nJ，驱动{w["drive_J"]*1e9:.9g}nJ、拖曳{w["drag_J"]*1e9:.9g}nJ，目标误差{w["relative_error"]*100:+.6g}%，{w["status"]}。峰值实际剪切{w["actual_peak_shear_Pa"]:.8g}Pa；最大Re_h={w["max_Re_h"]:.7g}。目标20–200块历史最多{s["max_medium_components"]}、最终{o["n_fragments_20_200"]}；最大块{max(q["final_free_sizes"],default=0)}点，细屑/脱离质量={o["fine_debris_mass_fraction_detached"]}，目标质量/初始={o["mass_fraction_target_initial"]}。真实二次分裂{q["secondary_count"]}次，3.4–3.6ms目标块存活{o["cohort_target_survivors_to_end"]}，完整目标={o["full_target"]}。')
    lines.append(f'{m}：累计新受载面{sf["unique_ever_loaded_new_faces"]}（自然新断裂{sf["unique_ever_loaded_natural_new_fracture_faces"]}，预切口重新张开{sf["unique_ever_loaded_reopened_precut_faces"]}）、末态新可达面{sf["final"]["accessible_new_faces"]}；新面净外功{sf["contribution"]["actual_net_work_on_new_accessible_faces_J"]*1e9:.8g}nJ。末态无已表示表面的游离分量{sf["final"]["free_components_without_surface"]}个/{sf["final"]["free_nodes_in_components_without_surface"]}节点。正式墙钟{s["wall_s"]:.3f}s、CPU{s["cpu_s"]:.3f}s、峰值RSS{s["peak_RSS_bytes"]/2**20:.3f}MiB，最大相对能量余额{s["max_relative_energy_residual"]*100:.6g}%。')
lines += ['碎块完全按原始存活bulk/interface及锚可达性定义；contact或几何接近不并块。20–200节点是当前离散工程目标，不是医学认证尺寸。新面工作量直接证明其承担载荷，但没有新模型关闭润湿的额外对照，不能单独归因碎块变化。',
    f'后续优化准入={status["optimization_admission"]}。本轮没有自动BoTorch、第三种物理工况、完整V3或完整泡CFD。高Re薄层代理、分辨率、未表示细屑表面及材料未标定限制仍在。',
    '每接受步保存标量、损伤变化和材料分量，关键事件与约50µs间隔保存FP64状态。实际端点力点乘真实native位移形成外功，R独立于D保留。两个初始化故障已完整归档：一次统计锚点索引导致各一个未接受无损伤微步，一次生成diff哈希防护在动力学前停止；修复后恢复同两组原始协议，未调参。没有4GiB/1800s预算。',
    f'详见[22问报告](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{TASK}/{RUN}/REVIEW_CLOT_SHEAR_MULTIVORTEX_FRAGMENTATION_V2.md)与[科研图表](temp_storage/clot_github_sync_20261015/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)。']
update(Path('/home/lzy/projects/clot_rupture.md'),'\n\n'.join(lines),O/'wsl_document','WSL_DOCUMENT_UPDATE',insert='# ⑥ 外加剪切')
entry=[f'# 当前可达裂面与多涡旋剪切入口（{RUN}）',f'当前任务 {TASK}；analysis/{TASK}/：shear_field_library → current_surface → native → run_cases → analyze。',
    '科学标签 SHEAR_MECHANICAL_REQUIREMENT_ONLY。支撑单元边界、实际开口、外部通道连通和当前遮挡的几何代理已接入；仍有未表示表面，不能称完整真实流体润湿。材料/native内核不改。',
    f'原生状态={status["status"]}；等功={status["work_match"]}；完整碎块目标={status["full_target"]}；优化准入={status["optimization_admission"]}。',
    f'先读 verification/{TASK}/{RUN}/FINAL_STATUS.json、REVIEW_CLOT_SHEAR_MULTIVORTEX_FRAGMENTATION_V2.md、DELIVERY_VALIDATION.json；图表 results/{TASK}/{RUN}/OPEN_RESULTS.html。',
    '历史A仅为预算参考；不追加正式工况或自动优化。所有接受证据与失败初始化记录均保留，新实验使用新RUN。','## 历史入口原文（保留）']
update(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md','\n\n'.join(entry),O/'entrypoints','ENTRYPOINT_UPDATE')
