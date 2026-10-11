"""Chinese scientific review and portable static figure index from validated evidence."""
from pathlib import Path
import sys,os,json,html
from design import js,save

def main(base,figs):
 q=js(base/'CRACK_PROPAGATION_AUDIT.json');B=q['M5B'];C=q['M5C'];status=js(base/'FINAL_STATUS.json');checks=js(base/'IMPLEMENTATION_CHECKS.json');detail=js(base/'EVENT_AND_STOP_DETAIL.json');E=js(base/'ENERGY_LEDGER_AUDIT.json');area=js(base/'CRACK_AREA_ATTEMPT.json');kin=js(base/'GRIP_KINEMATIC_AUDIT.json')['kinematics'];obs=js(base/'M5C/OBSERVER_CHECKS.json');resources=js(base/'RUN_RESOURCE.json');start=js(base/'provenance/SERVER_START.json');bp=B['states']['peak'];cp=C['states']['peak_tensile'];ce=C['states']['end'];be=B['matched_common_time'];r=Path(os.path.relpath(figs,base));report='REVIEW_LAMMPS_NOTCH_FRACTURE_GRIP_TRANSFER_V3.md';endd=detail['M5C']['states']['end'];pre=C['states']['pre_first_break']['regional_max_utilization_including_new_breaks'];contact=js(base/'CONTACT_FORCE_AUDIT.json')
 assert checks['status']=='PASS' and status['new_formal_cases']==1
 figures={'Fig1_grips_and_kinematic_bounds':'图1：真实夹具点阵、受控材料键与运动学界限','Fig2_first_break_and_competition':'图2：真实首断、断前风险及竞争断键簇','Fig3_exact_accepted_states':'图3：六个实际接受状态，位移无放大','Fig4a_boundary_force_work':'图4a：边界开度、反力、广义载荷与外功','Fig4b_complete_energy_and_dynamics':'图4b：完整保留时程的能量及动力学','Fig5_topology_and_area_limits':'图5：残留跨面连接、厚度条带及面积限制'}
 figurelinks='\n'.join(f'- [{title}]({(r/(n+".png")).as_posix()})' for n,title in figures.items())
 row=lambda label,b,c:f'|{label}|{b}|{c}|\n'
 table='|指标|M5B（旧数据复用）|M5C（本轮唯一新例）|\n|---|---:|---:|\n'
 table+=row('首断 / ms','5.2300；缺口','5.2555；夹具角点')+row('首个连贯断键簇 / ms','5.6685；触及缺口','5.4460；夹具附近')+row('峰值时间 / ms',f'{B["peak_time_s"]*1e3:.4f}',f'{C["peak_time_s"]*1e3:.4f}')+row('峰值名义牵引 / Pa',f'{B["peak_nominal_stress_Pa"]:.6f}',f'{C["peak_nominal_stress_Pa"]:.6f}')+row('峰值目标带断键',f'{bp["target_band_count"]}/{bp["broken"]}（{bp["target_band_fraction"]*100:.3f}%）',f'{cp["target_band_count"]}/{cp["broken"]}（{cp["target_band_fraction"]*100:.3f}%）')+row('峰值夹具区域断键',f'{bp["region_counts"][1]}/{bp["broken"]}（{bp["grip_fraction"]*100:.3f}%）',f'{cp["region_counts"][1]}/{cp["broken"]}（{cp["grip_fraction"]*100:.3f}%）')+row('同为7.05 ms时总断键',be['broken'],ce['broken'])+row('同为7.05 ms时夹具区域断键',f'{be["region_counts"][1]}（{be["grip_fraction"]*100:.3f}%）',f'{ce["region_counts"][1]}（{ce["grip_fraction"]*100:.3f}%）')+row('本例保留时程末态 / ms','24','7.05；按冻结规则停止')+row('本例末态双端轴向受控键断裂','13788','0')
 text=f'''# 缺口断裂与夹持传力审查 V3

Run ID：`{base.name}`。**实施和运动学检查通过，但本轮缺口优先起裂与可追踪主裂面目标未达到。** 最外端面轴向控制消除了 M5B 的规定性轴向拉断机制；M5C 却首先在夹具角点自然失效，7.05 ms 因单一目标裂纹带解释失效而按预注册规则停止。没有第二套边界试算，没有正式复验。`GC_UNAVAILABLE` 保持。

[静态图与证据入口](OPEN_RESULTS.html) · [分项状态](FINAL_STATUS.json) · [数值实施检查](IMPLEMENTATION_CHECKS.json) · [逐项对照CSV](M5B_M5C_COMPARISON.csv)

## 1. 身份、范围和冻结条件

服务器实际工作树 `/workspace/y_clot/project`，分支 `{start['branch']}`，HEAD `{start['head']}`；WSL 快照分支和起始脏文件清单在 `provenance/WSL_START.json`。没有使用较旧 GitHub 内容覆盖服务器。原 V1、V2 报告、接受态和 SHA 清单保持不变，没有重新计算剪切、M5B、Y血栓或 F7。

复用 `/workspace/y_clot/env.sh`、锁定 LAMMPS 30 Sep 2026 及原生补丁。G=940 Pa、K=33333.333333 Pa、ρ=1000 kg/m³、bulk s00=0.15、interface s00=0.140624859375、α=0、dx=50 µm、δ=150.5 µm、接触距离因子1、数值阻尼500 s⁻¹、dt=0.5 µs。2880点、132610条初始bulk键、2570条预删除键及原穿厚度预缺口完全相同；此试样没有实际interface键，保留参数并不等于执行界面标定。

只在实验室适配器添加 M5C 端面组、被动记录和接受态停止入口，默认旧算例调用不变。`provenance/macro_lab.before.py`、`macro_lab.server.before.py`、`macro_lab.changes.diff` 和 `launch_source/` 保留了变更与实际执行版本。`LAUNCH_SOURCE_IDENTITY.json` 记录正式启动前源码与设计 SHA；正式启动后这些文件均未修改。78项预记录身份中，唯一允许变化的是实验室适配器，其他77项（包括原生库、独立能量、生产模块和旧结果）逐项一致。

## 2. M5B 根因与新端面方案

详见 [只读根因报告](M5B_KINEMATIC_ROOT_CAUSE.md) 和 [运动学审计](GRIP_KINEMATIC_AUDIT.json)。M5B 的864点仅z受规定，2016点xyz完全自由，**没有三维全规定节点**。其 `fix move variable` 在位移和速度均为 NULL 的分量上执行力驱动 Verlet；锁定 `fix_move.cpp` 的初始、最终积分分支及真实输入证明组不重叠，没有把 NULL 当作冻结，也没有叠加 NVE。

对于两端z受规定的键，真实长度仍取决于自由x/y。严谨关系是

`s(t) >= abs((1 + disp(t)/H) * ΔZ) / r0 - 1`。

竖直键因此满足 `s>=disp/H`，H=850 µm。到6.842375006 ms，规定轴向分离达到s00阈值；实际竖直双受控键在6.7135–6.8425 ms间失效，可能因横向变形提前。旧双轴向受控键26308条，其中864条竖直；完整680 µm加载下下界保证失效7720条，旧真实事件中7720条全部断裂。其他键不能仅凭规定轴向位移声称三维伸长已知。

新方案只规定z=±475 µm两个最外端面的轴向位移±disp/2，288点受控、2592点完全自由，全部x/y自由；没有增加夹具质量、连接强度或新接触参数。3204条双轴向受控材料键都在同一端面，规定轴向差恒为0，**无相反端面之间的材料键、无必然轴向失效键**。其保持参考x/y的示意伸长为0；真实x/y仍由力学决定，并非三维长度被锁定。

原M5B均值夹具距离850 µm，新端面距离950 µm。保持实际均值开度 `680 µm ×(1−cos(πt/24ms))/2`；新最大名义应变为0.715789。旧最外层间开度是均值开度的950/850倍，故不把两条名义响应曲线当同一材料本构曲线。只使用一个 `run 0` 做原生组与零态检查；没有先试跑不同候选再择优。

## 3. 真实起裂、夹具污染与连续推进

{table}
M5C 第10511步（5.2555 ms）原生首断两条：body bond ID 130054 / 132610，端点ID分别2659–2660、2879–2880，参考中点 `(275, −275, 450)`、`(275, 275, 450)` µm。它们一端为受控端面、一端为自由内部；这是**力学驱动的端面角点失效**，不是双受控键的必然运动学拉断。第10510步真实断前状态：缺口最大利用率{pre['notch_tip']:.9f}，夹具{pre['grip_transition']:.9f}；首断时夹具已达1.000137。不得把后来缺口断键改称首断。

首个缺口区域断键在{detail['M5C']['first_notch_region_event']['time_s']*1e3:.4f} ms；中央厚度目标带首次断键在5.9015 ms。到终止时中央厚度目标带已有{endd['central_target_band_broken']}条断键，目标带内接近缺口的连通破坏体素簇包含中央断键{endd['tip_connected_targetband_central_broken']}条，说明损伤已进入内部。**但破坏体素云的连接不等于连续自由裂面**：全部12个50 µm厚度条带从缺口起算的跨面完全失连前缀均为0，未达到冻结的中央至少2条带推进100 µm门槛。原{endd['initial_cross_z0_bonds']}条跨z=0材料连接仍存活{endd['surviving_cross_z0_bonds']}条；2880点仍属于一个实际材料连通体。

冻结的V2区域分类基于参考键中点，优先次序仍为缺口→原夹具过渡区→其他自由表面→自由内部，没有把M5C夹具区域重新移动到端面来美化数字。末态5598条断键区域数为2674/1160/930/834；目标带2768条，非跨z=0键4344条（77.599%）；两端完全自由键断裂5114条，双端轴向受控键断裂0条。区域集合、受控端点集合、目标带与跨面集合有不同含义，不可互相替代。

因此，晚期确定性夹持拉断被消除、同一7.05 ms时刻夹具破坏减少，但峰值夹具污染明显高于M5B，缺口优先起裂优势丢失。不能将7.05 ms停止态与旧24 ms末态直接作等时改善比较，也不能外推未计算的M5C后期。

## 4. 提前停止与三维面积尝试

`M5C_DESIGN_FREEZE.json` 在正式启动前保存并哈希锁定。每100步（50 µs）检查当前接受态：非法/非有限状态，单步位移/dx>0.2，严重夹具损伤，目标带解释丢失，分离后低载持续1 ms，最后24 ms上限。当前轨迹只有第14100步首次触发：断键≥1000、目标带占比49.446%<50%、不跨z=0占比77.599%>60%。7.00 ms前一检查点未触发；[停止重放](EVENT_AND_STOP_DETAIL.json)核验全部141个检查点。没有读取未来状态选择停止，也没有裁掉更差阶段；完整0–7.05 ms数据保留。

面积尝试状态为 **`{area['status']}`**，跨面统计本身标为 `CROSSING_CONNECTIVITY_PROXY_ONLY`。最大连通破坏体素云跨厚度550 µm，PCA最小方差方向约(0.99994,0.00061,0.01127)，RMS厚度108.63 µm，远非薄的z=0唯一新增面。单个连通破坏体素云不说明存在单个物理裂面。中央推进、目标带比例、夹具占比、带宽与独立面积误差门槛均有未满足项。连接性面积代理0不表示真实新生表面积为0；真实新增面积为null，不能用断键数乘dx²，也未做平滑外壳、补缝或后删残余键。

## 5. 数值实施、能量、接触与惯性

14,101次原生接受态回调与0–14100步一一对应；所有22个保存的M5C坐标/速度状态均为双精度且拓扑精确匹配全事件。自由分量每一步验证位置及速度Verlet递推，并每100步对照原生最终x/v。最大位置递推误差{obs['max_free_position_recurrence_error_m']:.3g} m、速度递推误差{obs['max_free_velocity_recurrence_error_m_s']:.3g} m/s；原生块末自由速度差{obs['max_native_end_block_free_v_error_m_s']:.3g} m/s。最大实际每步位移/dx={obs['max_step_dx_over_dx']:.8f}，最大点速0.328654 m/s，没有达到冻结的步位移失稳门槛；这不是新增dt/2收敛证明。

沿用V1/V2完整独立LPS储能和原生力换算。第100/1000/5000步，负能量梯度与原生力相对差最高9.58×10⁻¹⁶；有限差分差最高3.00×10⁻⁸。原native bond PE少计偏能，在这三个状态相对独立能量差−14.05%、−47.47%、−47.61%，原样保留作为诊断，**没有将其当完整LPS储能，也没有修改轨迹来迁就native PE**。

外功始终逐实际规定DOF累积：`W=−∫Σ f_material·du_prescribed`，`Q=−Σf_material·∂u_prescribed/∂disp`。M5C可化为 `Q=(Rtop−Rbottom)/2`，实际误差≤2.72×10⁻²⁰ N；M5B可变形三层不能直接用顶部合力乘均值开度替代。账本只记传给材料的外功与自由DOF动能，未将规定运动装置动能混入。

M5C末态（nJ）：外功8.245132，bulk储能4.299974，自由动能0.525915，数值阻尼0.220010，净带符号拓扑项3.198917，保守接触储能0.0000670，残差0.0002495。拓扑正增量总和3.199122 nJ、负增量−0.000205664 nJ，负事件2次；固定构型拓扑储能变化没有截为非负，也未用能量残差反推断裂耗散。最大独立余额/最大外功 **0.011702%**。

真实接触能每步保留，峰值0.00109858 nJ。原生断键步膨胀量存在已知更新滞后，不能在这些步直接把“总力−断后独立bulk力”叫接触力；非零接触断键记录145行显式标为不可分解。另按锁定pair公式重算全部保存状态的保守接触，与原生接触能逐态完全一致：末态有2对接触，保存状态最大节点接触力0.197144 µN；**逐步接触力峰值未保存，仍未知**。见 [接触力审查](CONTACT_FORCE_AUDIT.json)。

首断前一步K/U=0.485316%，首断至峰值最高1.882282%；峰后至停止升到12.230649%。启动阶段0.092 ms还出现188.34%，此时U仅1.1932×10⁻¹⁵ J，不以接近零分母宣称全程准静态。数值阻尼500 s⁻¹不代表血栓生物黏弹耗散。分阶段账本见 [ENERGY_LEDGER_AUDIT](ENERGY_LEDGER_AUDIT.json)，不存在“全程准静态PASS”。

## 6. 科学结论、AI标签与唯一后续事项

本轮没有可靠的新物理裂面面积，也没有独立识别与主裂纹推进对应的不可逆能量分子，**不能给出试样相关表观断裂功/面积候选，更不能给材料固有Gc**。试样外功、带符号拓扑项仍是可审查的数值结果，不能改名为已认证断裂能。`APPARENT_FRACTURE_WORK_STATUS=UNAVAILABLE`，`MATERIAL_GC_STATUS=GC_UNAVAILABLE`。

当前LPS适用于明确边界与离散条件下的小应变响应、断键事件和边界敏感性研究；本轮证明去掉必然运动学拉断仍不足以得到合格缺口裂纹实验。它不能单独证明“缺软化”就是根因，没有充分证据升级软化本构。继承V1网格敏感性：dx50→40 µm曲线差约44.15%、净拓扑能量差约48.13%；本轮未复做网格或时间步研究，也未升级求解器。

可用于未来AI的内容：带边界/点距/步长/可信度标签的原生事件、首断位置、连通分量、试样载荷曲线和数值能量分项；本例应标为“端面加载、夹具先断、未形成可识别单裂面”的负例。不能认证材料Gc、网格无关碎片大小/强度、物理裂面面积、全程准静态、生物患者真值或微泡控制真值。

**下一轮唯一优先项：{status['next_priority']}** 本轮到此停止，未自动启动后续研究。

## 7. 计算预算、交付与WSL文档

实际正式 M5C×1、原生零步检查×1、针对性正式复验×0。唯一边界在运动学预检查后冻结。独立后处理程序的属性访问及native PE校验口径更正均不涉及原生重算，详见 `POSTPROCESS_NOTES.json`。

CPU cgroup仍为15.36核；使用1 MPI rank、OMP/BLAS线程1，GPU未参与，未安装/升级软件，未做压力测试。正式进程墙钟{resources['wall_s']:.3f} s、用户CPU{resources['child_user_s']:.3f} s、最大RSS {resources['max_rss_KiB']/1024:.2f} MiB。独立supervisor进程正常退出，未干扰其他任务。原始数据服务器目录：`{base}/M5C`；图目录：`{figs}`。

实际 WSL 科研文档 `/home/lzy/projects/clot_rupture.md` 只在实施证据通过后做最小增补；写入状态与前后SHA见 [写入回执](provenance/WSL_DOCUMENT_UPDATE.json)，[前文备份](provenance/clot_rupture.before.md)、[逐行diff](provenance/clot_rupture.changes.diff)、[写后副本](provenance/clot_rupture.after.md)。负面科学结果不会写成材料已标定。

大型原始data、逐步CSV、每0.5 ms间隔原生状态留服务器；WSL保留全精度关键状态、真实事件、精简数值源、报告和图。`FILE_MANIFEST.json`给出全部正式交付身份，`provenance/COMPACT_DOWNLOAD.json`说明实际回收范围。没有Git推送、重置、清理或旧结果覆盖。

{figurelinks}

所有图由真实点阵、接受态及原生事件绘制，变形比例1。图1/2/5标出的参考空间用于定位；图3为真实变形坐标，稠密存活键只在注明的中央y切片展示，节点未合并、碎块未放大；没有凭图像推导面积。

## 8. 可复核命令与参考语义

正式命令、源文件SHA、原生输入分别在 `RUN_RESOURCE.json`、`provenance/LAUNCH_SOURCE_IDENTITY.json`、`M5C/in.lammps`。正式结果目录禁止复用重跑；如仅复核分析，可激活环境后执行 `analysis/notch_grip_transfer_v3/analyze.py <V3目录> <V2目录>`、`contact_audit.py <V3目录>`、`detail.py <V3目录>`，只读取原生数据、写入本轮分析文件。原生逐块命令保存在 `M5C/run_commands.txt`，后续`run ... pre no`不会重复setup或初始回调，这一点已对照锁定run.cpp与实际步序。

锁定源码是实际行为依据；公开语义参照 [LAMMPS fix move](https://docs.lammps.org/fix_move.html) 与 [fix nve](https://docs.lammps.org/fix_nve.html)。这里只用其解释NULL分量与Verlet积分，不以在线文档替换锁定版本。原Liu正文和补充材料本轮沿用已有证据；没有拿新文献数值给本次断裂能背书。
'''
 (base/report).write_text(text)
 rootcause=base/'M5B_KINEMATIC_ROOT_CAUSE.md';s=rootcause.read_text();s=s.replace('图1将展示','图1展示');link=f'\n[图1：真实夹具与受控键]({(r/"Fig1_grips_and_kinematic_bounds.png").as_posix()})\n'
 if link not in s:s+=link
 rootcause.write_text(s)
 header='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>夹持传力与缺口断裂 V3</title><style>body{font:17px/1.65 system-ui,sans-serif;background:#eef2f5;color:#173042;margin:0}main{max-width:1280px;margin:auto;padding:28px}section{background:white;border-radius:10px;padding:24px;margin:20px 0}a{color:#146782}img{width:100%;height:auto}code{overflow-wrap:anywhere}.result{border-left:6px solid #b44c55}table{border-collapse:collapse}td,th{padding:8px 18px;border-bottom:1px solid #d3dde3}</style><main>'''
 body=f'<h1>夹持传力与缺口断裂 V3</h1><p>Run {base.name} · 正式M5C 1次 · 0次重试 · 材料和生产参数不变</p><section class="result"><b>运动学缺陷已去除；缺口优先起裂与连续主裂面目标未达到。Gc仍不可识别。</b><p>首断5.2555 ms发生在端面角点；7.05 ms按冻结规则停止。中央有断键，但没有连续跨面完全失连推进，仍为一个2880点连通体。</p><p>'
 links=[(report,'中文审查'),('FINAL_STATUS.json','分项状态'),('M5B_M5C_COMPARISON.csv','对照CSV'),('IMPLEMENTATION_CHECKS.json','实施验证'),('M5B_KINEMATIC_ROOT_CAUSE.md','旧边界根因'),('M5C_DESIGN_FREEZE.json','运行前冻结设计'),('GRIP_KINEMATIC_AUDIT.json','运动学审计'),('CRACK_PROPAGATION_AUDIT.json','裂纹推进审计'),('CRACK_AREA_ATTEMPT.json','面积尝试'),('ENERGY_LEDGER_AUDIT.json','能量账本'),('CONTACT_FORCE_AUDIT.json','接触力'),('EVENT_AND_STOP_DETAIL.json','停止判据重放')]
 body+=' · '.join(f'<a href="{url}">{label}</a>' for url,label in links)+'</p></section>'
 for n,title in figures.items():
  href=(r/(n+'.png')).as_posix();body+=f'<section><h2>{title}</h2><a href="{href}"><img loading="lazy" alt="{title}" src="{href}"></a>'
  if (figs/(n+'.svg')).exists():body+=f'<a href="{(r/(n+".svg")).as_posix()}">矢量SVG</a>'
  body+='</section>'
 body+='<section><h2>数据来源与文档写入</h2><p>科学图使用原生接受态，不生成假想裂面；体素和跨面连接是空间诊断。大型逐步原始文件保存在服务器，本地提供关键精度证据。</p><p><a href="provenance/clot_rupture.before.md">文档备份</a> · <a href="provenance/clot_rupture.changes.diff">文档diff</a> · <a href="provenance/WSL_DOCUMENT_UPDATE.json">写入与SHA回执</a> · <a href="FILE_MANIFEST.json">文件清单与SHA</a></p></section></main></html>'
 (base/'OPEN_RESULTS.html').write_text(header+body);(figs/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><a href="'+os.path.relpath(base/'OPEN_RESULTS.html',figs)+'">打开V3全部审查结果</a>')
 print('Chinese review and static index written')

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
