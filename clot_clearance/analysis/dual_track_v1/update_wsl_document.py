"""Single serial WSL document update with backup, concurrent-change guard and readback."""
from pathlib import Path
import hashlib,json,datetime,difflib,os,re
cfg=json.loads(Path('/home/lzy/projects/.dual_track_current.json').read_text());root=Path(cfg['root']);out=Path(cfg['local'])/'document_edit';doc=Path('/home/lzy/projects/clot_rupture.md');assert not (out/'AUDIT.json').exists(),'Already applied; do not rerun'
a=json.loads((root/cfg['A']/'ANALYSIS_SUMMARY.json').read_text());b=json.loads((root/cfg['B']/'ANALYSIS_SUMMARY.json').read_text());level=b['selected_level'];t=b['traction'][level];m=b['flows'][level]['metrics'];comp=b['comparison'][level]
raw=doc.read_bytes();before=raw.decode();digest=lambda x:hashlib.sha256(x).hexdigest();out.mkdir(parents=True,exist_ok=True);(out/'before.md').write_bytes(raw);started=datetime.datetime.now(datetime.timezone.utc).isoformat();sha_before=digest(raw)
old='⚠️ 这是**没有血栓**时的流场。血栓挡住血流后流场会变，这一点还没算。';new='⚠️ 下面原有参数和历史碎裂、输运轨迹使用的是**没有血栓**时的背景流场。本轮另行完成真正含 V3 血栓的 CFD，结果见本节新增说明；旧轨迹没有重新计算，也没有追溯替换成新含栓解。';assert before.count(old)==1;after=before.replace(old,new)
after=after.replace('两个出口压力都设为 0。','两个出口采用原生零参考牵引条件（Neu、Value=0），不等于出口每个节点的计算压力都恰好为零。',1)
addition_b=fr'''

**本轮真正含栓背景流（CFD 双线任务）**：直接使用原生 Y 管和 V3 原连续高度函数，沿母管 +z 内壁构造流体域，未用点云凸包、球或椭球替代。CFD 三角化血栓体积为 {b['geometry']['CFD_triangulated_clot_volume_m3']*1e9:.9f} mm³；原连续口径约 0.271574 mm³，1927 个材料点的积分体积为 0.240875 mm³，三者分别记录。壁锚点不进入流体域。共享接触线无间隙，内部附着面不作为流体载荷面。

保持原入口抛物线剖面、平均 5 mm/s、流量约 1.581244 mL/min、两个零参考牵引出口及水样参数不变，血管壁和血栓暴露面都无滑移。实际完成基线 35,180 点/177,263 四面体，以及一次局部加密 36,599 点/187,162 四面体；两者均通过原生线性/非线性检查、五个连续稳态区间和冻结质量门槛。加密解质量不平衡为 {m['epsilon_mass']:.3e}，最大流速 **{m['velocity_max_m_s']*1e3:.6f} mm/s**，旧无栓值为 10.074017 mm/s；两出口比例为 {m['outlet_fractions']['OUTLET_A']:.6f}/{m['outlet_fractions']['OUTLET_B']:.6f}。母管 x=−9 至 −3 mm 的截面平均压降从 {comp['parent_pressure_drop'][0]['drop_Pa']:.6f} 增为 {comp['parent_pressure_drop'][1]['drop_Pa']:.6f} Pa。

取血栓指向流体的法向，按 $\boldsymbol\sigma=-p\mathbf I+\mu(\nabla\mathbf u+\nabla\mathbf u^\top)$、$\mathbf t=\boldsymbol\sigma\mathbf n_{{\rm clot}}$ 提取真实暴露面载荷，压力与剪切分开。加密解的面积平均表压为 **{t['pressure_Pa']['mean']:.6f} Pa**，切向剪切平均值/P95 为 **{t['shear_Pa']['mean']:.6f}/{t['shear_Pa']['p95']:.6f} Pa**；总流体力约为 ({t['total_force_N'][0]*1e6:.6f}, {t['total_force_N'][1]*1e6:.6f}, {t['total_force_N'][2]*1e6:.6f}) µN，大小 **{t['total_force_norm_N']*1e6:.6f} µN**，不包含固体壁面粘附反力。一次加密的平均剪切/P95 变化约 6.22%/5.65%，积分力向量变化约 0.067%，通过事先冻结的 15% 工程筛查；共享表面三角化与局部尖角峰值尚未证明渐近收敛。表压不等于绝对血压，当前结果也不是微泡 streaming 或双向 FSI；尚未映射到 LAMMPS 材料点，更未改动下面的工程剪切脉冲。
'''
anchor='> **出处**：';assert anchor in after;after=after.replace(anchor,addition_b+'\n\n'+anchor,1)
addition_a=fr'''

**本轮宽薄周期试样（PVCE 双线任务）**：实际 X×Y×Z=5×0.4×1 mm，16,000 个点、909,392 条完整参考键，其中 131,196 条跨 Y 周期。Y 单元厚度固定，宏观 Y 应变为零，代表理想化周期材料层；X 有限自由，z=±475 µm 的最外层仅控制轴向位移。原生周期影像、ghost、参考键、固定完整 m、膨胀量、力和独立储能均通过静态一致性检查；并非 Liu 完整实验，也不是有限厚度真实血栓。

相同 2 µm 总开度下，a=1.50/1.55/1.60 mm 三态均通过原冻结平衡条件，累计预切 16,416/16,976/17,536 条键，每增量 560 条，$B\Delta a=2.0\times10^{{-8}}$ m²。$G_{{\rm two}}={a['G']['G_two_J_m2']:.8e}$、$G_{{\rm three}}={a['G']['G_three_J_m2']:.8e}$ J/m²；经验数值分辨尺度分别约 ±2.83×10⁻⁸、±5.79×10⁻⁸ J/m²，不是完整离散误差界。比旧小试样约高 15.43%/15.73%，但宽度、裂纹长度、厚度和周期条件同时改变，不能分离归因。**这是试样相关 G_drive，材料 Gc 仍未识别。**

旧 73.5606 µm 小试样保存态仅做一次局部残差修正，最大原生自由力由 1.1291×10⁻¹² 降至 2.4431×10⁻¹⁹ N，原平衡门槛全部通过；历史 V2 的失败记录保留。随后只计算一个宽薄 A0 高载候选（72.554732 µm），自由力通过，但末十次位置变化 1.1061×10⁻¹⁰ m 超过 1×10⁻¹⁰ m，1 pm 分量差分也未通过冻结门槛，故**候选未验收并停止**。其最大 s/s00≈0.91685 位于夹具附近，只能作空间诊断，不能称为已首先起裂。没有第二候选、高载 A1/A2、新动态断裂或 G_high；只读较大差分步长检查支持舍入解释，但没有据此放宽门槛。本线正式静态调用为 5/8，仍无新增证据支持修改 LPS 本构。
'''
anchor='---\n\n# ⑥ 外加剪切';assert after.count(anchor)==1;after=after.replace(anchor,addition_a+'\n\n'+anchor)
after=after.replace('最后严格有效的 A0 仍是 V1 的 2 µm 平衡态。','在 V2 结束时，最后严格有效的 A0 仍是 V1 的 2 µm 平衡态；本轮恢复结果另见后述。',1)
oldnext='下一轮唯一优先项是保持材料、原生力和验收门槛不变，修复并独立验证 A0 高载静态求解器的浮点停机问题，达标后再恢复首达阈值区间搜索；不自动启动。'
links=f"本轮 [宽薄 PVCE 审查]({root/cfg['A']/'REVIEW_LPS_WIDE_PURE_SHEAR_PVCE_FEASIBILITY_V1.md'})、[真正含栓 CFD 审查]({root/cfg['B']/'REVIEW_Y_CLOT_PRESENT_CFD_BASELINE_V1.md'}) 和 [两线共 14 组真实图与数据索引]({Path(cfg['local'])/'OPEN_RESULTS.html'}) 分别记录低载成果、高载停止及含栓牵引结果。下一轮唯一优先项是验证含栓 CFD 暴露面牵引向现有 LAMMPS 材料点的离线保力、保力矩映射，先核对单位、坐标、面身份和面积权重，不直接启动生产加载或双向耦合；不自动启动。"
assert after.count(oldnext)==1;after=after.replace(oldnext,links)
assert re.findall(r'^# .*$',before,re.M)==re.findall(r'^# .*$',after,re.M)
assert all(ord(c)>=32 or c in '\n\t\r' for c in after)
# Explicitly reject stray literal tabs in the new scientific prose.
assert '\t' not in addition_a and '\t' not in addition_b
(out/'proposed.md').write_text(after);(out/'change.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='clot_rupture.md.before',tofile='clot_rupture.md.after')))
assert digest(doc.read_bytes())==sha_before,'Concurrent edit; proposed text preserved, actual document untouched'
tmp=doc.with_name('.clot_rupture.dual_track.tmp');assert not tmp.exists();tmp.write_text(after);os.chmod(tmp,doc.stat().st_mode);assert digest(doc.read_bytes())==sha_before;os.replace(tmp,doc);readback=doc.read_bytes();assert readback==after.encode();(out/'after.md').write_bytes(readback)
audit=dict(status='PASS',actual_WSL_path=str(doc),started_utc=started,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),before_sha256=sha_before,after_sha256=digest(readback),backup_sha256=digest((out/'before.md').read_bytes()),readback_exact=True,concurrent_checks=2,headings_unchanged=True,sections_changed=['① 背景血流','⑤ 失效与 PVCE','碎块统计和本次结果：索引/下一轮建议'],history_preserved=True,old_trajectories_not_recomputed_or_relabelled=True)
(out/'AUDIT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2));print(json.dumps(audit,ensure_ascii=False))
