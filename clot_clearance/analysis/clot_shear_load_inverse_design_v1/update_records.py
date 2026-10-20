"""Stdlib-only local record update; immutable before/after evidence."""
from pathlib import Path
import json,hashlib,difflib
R=Path(__file__).resolve().parents[2];TASK='clot_shear_load_inverse_design_v1';RUN='20261017T165129Z';O=R/'verification'/TASK/RUN
def read(p):return json.loads(p.read_text())
def sha(b):return hashlib.sha256(b).hexdigest()
def update(path,block,folder,name,insert=None):
    before=path.read_bytes();old=before.decode();assert block.splitlines()[0] not in old
    folder.mkdir(exist_ok=True);backup=folder/(path.stem+'.before.md');assert not backup.exists();backup.write_bytes(before)
    i=old.index(insert) if insert else 0;new=old[:i]+block+'\n\n'+old[i:]
    assert path.read_bytes()==before;path.write_text(new);after=path.read_bytes()
    assert after.decode()==new and new.replace(block+'\n\n','',1)==old
    (folder/(path.stem+'.after.md')).write_bytes(after)
    (folder/(path.stem+'.diff')).write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=path.name+'.before',tofile=path.name+'.after')))
    proof=dict(status='PASS',path=str(path) if name=='WSL_DOCUMENT_UPDATE' else str(path.relative_to(R)),before_sha256=sha(before),after_sha256=sha(after),history_preserved_exactly=True,readback_equal=True)
    (O/(name+'.json')).write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n');print(name,proof['after_sha256'])
assert read(O/'DELIVERY_VALIDATION.json')['status']=='PASS'
q=read(O/'ANALYSIS.json');status=read(O/'FINAL_STATUS.json');resources=read(O/'RUN_RESOURCES.json')
lines=[f'**当前表面随动剪切需求逆向研究（{RUN}，2026-10-17）**',
    '为先确定机械剪切需求，本轮不求解完整微泡流场，完成单向与双涡旋两个低维理想化剪切模式。统一标记 SHEAR_MECHANICAL_REQUIREMENT_ONLY；不是声压或R(t)直接换算的真实Streaming。用户FEM网格claude.pdf已读，未采用弱格子、额外预切、大规模扫描或完整V3建议。',
    '双涡旋从相反符号Gaussian表面流函数定义局部切向流；Mobadersany & Sarkar近壁泡研究提供涡流形态动机，其轴对称环涡并非本轮平面Gaussian涡，未移植应力。7维参数冻结名义1000Pa、平面方向0rad、中心(0,0)、宽125µm、间距300µm、驱动2.8ms；共同2ms升、2.4ms降、2.8ms停驱、3.6ms结束。µ=.001Pa·s、h=25µm，外缘上限25m/s，Re_h上限625，只属机械代理。',
    '保留1792材料点+336锚点、95666原始键/2736原预切、均匀G/K/s0/sc/ρ和接触。182个原节点顶面三角形每次力回调更新当前位置、面积、法向和速度。τ=(µ/h)P(u_outer−v_surface)；三个顶点正1/3分配保持离散力、当前力矩和伴随虚功。局部合力不超过正面力预算，初始双涡旋节点与解析方向差最大4.6266°；空间积分尚未网格收敛，旧V3有符号仿射映射未修复/准入。',
    '面片三条原始边存活且几何非退化才保留；失去面片支撑不等于人为删除材料键。脱离不是力开关：保留面片按当前空间和相对流速受力，离开驱动区后仍有拖曳。新裂面、其他原外表面及没有保留面片碎片的水动力未解析，NEW_CRACK_FACE_WETTING_UNRESOLVED；不把全部失效键伪装成润湿裂面。',
    '原生LPS内核不改。速度相关外力以第二半步速度一致解接入旧VV；dt=.25µs/事件最小.03125µs，拒绝预测步在native提交之前缩步，Δd≤.1。实际drift积分外功，边失效根分段，驱动功、被动拖曳耗能、材料D、接触激活J及数值R分列；1%余额只检验数值模型。没有额外数值阻尼。']
for k,v in q.items():
    e=v['final'];lines.append(f'{k}：3.6ms末态本体/界面自然完全失效{int(e["bulk_failed"])}/{int(e["interface_failed"])}，游离{int(e["detached_nodes"])}节点/{int(e["components"])}分量，历史最大块{v["history_max_free_nodes"]}节点；20–200块历史最大/最终{v["history_max_N20_200"]}/{int(e["N20_200"])}，目标带质量/初始={e["mass_fraction20_200"]}，细屑/脱离质量={e["fine_debris_mass_fraction_detached"]}，完整探索目标={v["exploratory_target_met"]}。实际净外功{e["W_J"]*1e9:.8g}nJ，驱动{e["W_drive_J"]*1e9:.8g}nJ、拖曳{e["W_drag_J"]*1e9:.8g}nJ，D={e["D_J"]*1e9:.8g}nJ，最大相对余额{v["max_relative_energy_residual"]:.8g}。原受载参考面积最终保留{e["active_reference_area_fraction"]*100:.6g}%，游离节点中未直接列为保留面片顶点（含受力块内部节点）{int(e["free_nodes_without_resolved_surface"])}个。')
lines += [
    '面片覆盖统计另外按整个分量计算，见FRAGMENT_SURFACE_COVERAGE_AUDIT.json；非表面顶点数包含大块内点，不能当成未耦合碎片数。完整每接受步的原始ID损伤/存活图和锚可达性定义碎片，contact不并块。输出六档数量/质量、参考等效直径、当前三维跨度、内部存活/密度和有限卸载观察寿命。三历史诊断区及阈值/谱分区用于只读残余桥接分析；代理不当作真碎片，未删键。失效网络合并不能证明闭合三维断面。',
    f'科学结论：{status["SCIENTIFIC_CONCLUSION"]} 优化准入={status["OPTIMIZATION_READINESS"]}；本轮两组主要原生、0新优化trial、0完整V3。不同场同材料/包络/名义cap但不同实际外功，不能宣称等能量。旧死载与本轮同时改变场和速度闭合，不能单独归因。',
    '下一轮先审查剪切实际分布、存活三维桥接和表面支持损失，未获目标或明确近分离信号时不盲目启动BoTorch。混合模式本构须独立证据，不因未达标直接修改。声学对接须泡/血栓坐标与动态间隙、泡尺寸相位和径向/非球形响应、一阶复速度梯度、二阶平均流、介质和壁面边界及可追溯表面应力/速度；V3还需局部映射、当前/新裂面和碎片流体覆盖准入。',
    f'两例并行：累计算例墙钟{resources["native_wall_s_sum"]:.3f}s，CPU{resources["native_cpu_s"]:.3f}s；最大单例RSS{resources["max_single_case_peak_RSS_bytes"]/2**20:.3f}MiB。原生结果、能量与所有接受步原始图回放通过。没有4GiB/1800s硬限制。',
    f'详见[20问报告](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{TASK}/{RUN}/REVIEW_CLOT_SHEAR_LOAD_INVERSE_DESIGN_V1.md)及[全部图表](temp_storage/clot_github_sync_20261015/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)。']
update(Path('/home/lzy/projects/clot_rupture.md'),'\n\n'.join(lines),O/'wsl_document','WSL_DOCUMENT_UPDATE',insert='# ⑥ 外加剪切')
entry=[f'# 当前随动剪切开发入口（{RUN}）',f'当前任务 `{TASK}`，入口 `analysis/{TASK}/`：shear_field_library.py → current_surface.py → native.py → run_cases.py → analyze.py。',
    '两组3.6ms小试样原生完成。SHEAR_MECHANICAL_REQUIREMENT_ONLY；182原始顶面随动、正1/3载荷映射和相对速度拖曳。新裂面润湿、未表示碎片表面流体力及旧V3局部映射仍未解决。材料/native内核不改，无压头/弱格子/人为删键。',
    status['SCIENTIFIC_CONCLUSION'],f'准入：{status["OPTIMIZATION_READINESS"]}；没有优化trial或完整V3。',
    f'先读 `verification/{TASK}/{RUN}/FINAL_STATUS.json`、`REVIEW_CLOT_SHEAR_LOAD_INVERSE_DESIGN_V1.md`、`DELIVERY_VALIDATION.json`；图表 `results/{TASK}/{RUN}/OPEN_RESULTS.html`。',
    '所有新计算使用新RUN，已封存结果只读。', '## 历史入口原文（保留）']
update(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md','\n\n'.join(entry),O/'entrypoints','ENTRYPOINT_UPDATE')
