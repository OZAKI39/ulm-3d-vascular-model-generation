"""Generate the requested Chinese research report from sealed numerical evidence."""
from analyze import *
import html
def fmt(x,scale=1):return 'NOT_OBSERVED / UNAVAILABLE' if x is None else f'{x*scale:.8g}'
def main():
 q=read(O/'ANALYSIS.json');s=read(O/'NATIVE_SHEAR_CASE_RESULTS.json');status=read(O/'FINAL_STATUS.json');mapping=read(O/'SURFACE_MAPPING_AUDIT.json');bridge=read(O/'CRACK_BRIDGE_AUDIT.json');resources=read(O/'RUN_RESOURCES.json');protocol=read(O/'SHEAR_PROTOCOL.json');flow=read(O/'STREAMING_MODEL_SPEC.json');ref=load_npz(O/'provenance/reference.npz');load=ShearLoad(ref,'stokes');a=q['stokes'];b=q['inverse'];fluidcheck=read(O/'FLUID_SHEAR_CHECK.json')
 answers=[]
 def add(question,answer):answers.append(f'### {len(answers)+1}. {question}\n\n{answer}\n')
 add('本轮是否真正恢复了流体剪切加载？','是。新增 `acoustic_streaming_shear.py`、`streaming_interface.py` 和 `shear_loading.py`，实际生成黏性剪切场，经原始节点三角形映射进入已验证 LAMMPS `fix external` 的 N 单位接口，完成两组原生不可逆动力学。没有使用或开发刚性压头，生产本构、MPI/restart 和原始预裂纹均未改变。')
 add('当前声流模型来自哪一种物理方程？',f'本轮计算的是 **PRESCRIBED_STOKES_SHEAR_TEMPLATE / PHYSICS_APPROXIMATION**，不是一阶声学振荡或二阶平均声流。取 u=(sin(ky)v(z,t),0,0)，求解 v_t=ν(v_zz−k²v)，压力增量为零。流体从静止开始，底壁静止无滑移，25 µm 上边界速度为 4 mm/s×q(t)×sin(ky)，横向周期700 µm、沿x无限延伸。该速度严格无散且对流项为零。sinh稳态解加正弦扩散模态给出时变解；独立有限差分最大剪切差/μU/H={fluidcheck["independent_FD_max_error_over_muU_H"]:.4g}。有限血栓边缘、变形反馈不在这个流体解内。')
 add('MEVO/APECSS能提供哪些实际输入？','重新打开并核对了10项实际源文件/数据哈希。F1提供原几何的复声压及梯度；F10N/F10S提供2.5 MHz均匀水控制代理、源幅值与相位约定；APECSS表包含8个有效自适应叶中的体积谐波，半径1.779–1.781 µm、频率2.499–2.501 MHz、声压0.1–40000 Pa的叶内插值支持；F10S提供原坐标系质心轨迹。声源Pa等效量不是PZT电压或校准面压。半径谐波/质心运动均不等于血栓附近流体速度。APECSS黏度2 mPa·s，本轮独立历史Stokes模板1 mPa·s，不能当成同一介质的完整耦合。')
 add('是否完成微泡诱导Streaming计算？','没有。缺少微泡与当前血栓的坐标配准/真实间隙、近壁一阶速度和梯度、完整振荡边界条件及一致二阶边界驱动。没有把声压或R(t)直接换算为三维声流剪切，也没有编造泡位置。新入口遇到微泡模型请求明确返回 STREAMING_LOAD_NOT_YET_DETERMINED 和缺失字段。现有径向模型可生成R(t)，但补齐R(t)本身不足以补齐近壁流体问题。')
 add('实际血栓表面剪切应力大小是多少？',f'在本轮实际小试样参考三角面片上，A的最大剪切为 **{a["peak_shear_Pa"]:.9g} Pa**，B为 **{b["peak_shear_Pa"]:.9g} Pa**。B只标记 INVERSE_SHEAR_REQUIREMENT / MECHANICAL_REQUIREMENT_ONLY。法向压力增量和法向牵引均为零。A不是测得或微泡求得的真实血栓表面应力。')
 add('剪切热点在空间上如何分布？','只选择一种模式：实际顶面上 sin(2πy/700 µm)×eₓ，y两侧方向相反，中心平滑过零。横向350 µm试样跨度决定一个低阶连续模态；它产生近零净力、非零z向力矩，主要检验反向流体剪切造成的扭转，而非沿x整体推动。没有根据目标裂纹或待断键布置载荷。上壁2 ms平滑上升、2.4 ms开始下降、2.8 ms完全停驱，之后黏性流自然衰减，原生观察到3.6 ms。B沿用完全相同形状/时程，仅缩放牵引。历史1000Pa均匀+x剪切在同一参考试样、同一2.275e−7m²顶面产生过882节点块；历史为2ms上升后保持到2.5ms、四节点面片分配。本轮改变了空间方向组织、扩散建立过程、三角积分和卸载观察，不能当成仅改方向的单因素对照，详见HISTORICAL_SHEAR_CONTEXT.json。')
 add('哪些原始材料点真正受载？',f'1792个材料点、336个底部固定锚点的既有三维试样；顶面 **{len(load.loaded_ids)}个原始材料ID** 接收非零剪切，锚点不受外加载。实际面积 **{load.area.sum():.9g} m²**，{len(load.tri)}个原始节点三角形。所有ID、参考坐标、面片、面积、法向、加载向量见 `SHEAR_TRACTION_FIELD.npz`；每个接受态的标量历史和关键FP64态保留实际力，完整力可由冻结空间模式及各接受时刻精确重建。')
 v=mapping['v3'];c=v['local_comparison']
 add('映射是否出现明显局部放大？',f'小试样三角形1/3形函数无外插，精确保持离散总力、力矩和伴随虚功，未引入反向力。其交替三角剖分的节点分摊面积相对规则控制面积为2/3–4/3；按映射面积归一化的节点剪切与节点位置解析值最大误差为峰值的2.8112%。因此局部权重无符号放大不等于点力或破坏结果已网格收敛。真实V3的3611个CFD三角面上，已有仿射权重范围[{v["signed_affine"]["min_weight"]:.4f},{v["signed_affine"]["max_weight"]:.4f}]，最大列绝对和{v["signed_affine"]["max_column_abs_sum"]:.4f}；本模式下节点力模长总和/面力模长总和{v["signed_affine"]["sum_node_force_norm_over_sum_face_norm"]:.4f}。与独立构建的非负局部核相比，共同活跃{c["common_active_nodes"]}点中{c["nodes_direction_difference_over_90deg"]}点合力夹角>90°，最大节点力/非负局部力预算{c["signed_force_over_positive_local_force_budget_quantiles"][-1]:.4f}。这直接来自合力，非由负权重比例推断。非负核保持总力/功率，但牺牲力矩（误差/总绝对矩={v["nonnegative_local_kernel"]["moment_error_over_total_abs_moment"]:.4g}），因此两者都不能凭全局守恒认证V3局部物理精度。没有截断旧负权重、没有运行V3。')
 add('A是否引起真实softening？',f'A末态新软化键 {int(a["final"]["softened"])}，历史最大正伸长κ={max(r["max_kappa"] for r in histories(O/"cases/stokes")):.9g}；首次软化时间={fmt(a["first_softening_s"],1e3)} ms。B首次软化={fmt(b["first_softening_s"],1e3)} ms。这是具体试样/模式结果，不是普适材料阈值。')
 add('是否出现完整断键？','；'.join(f'{name}: 本体{int(x["final"]["bulk_failed"])}、界面{int(x["final"]["interface_failed"])}，首次完整断裂{fmt(x["first_failure_s"],1e3)} ms' for name,x in q.items())+'。初始2736预切键不计自然断裂；每条变化由原始bond ID记录；完整失效键排除力学贡献并在原生连接表中移除。')
 add('是否形成多个空间独立损伤区域？','；'.join(f'{name}: 自然bulk/interface失效网络诞生{x["inclusive_failed_material_networks"]["n_births"]}次，最多同时{x["inclusive_failed_material_networks"]["max_simultaneous_networks"]}个网络，bulk诞生中心最大间距{x["crack_networks"]["max_birth_centroid_separation_m"]*1e6:.6g} µm' for name,x in q.items())+'。失效边共享端点的网络分离/合并是离散诊断，不自动等于连续裂面的闭合或贯通。')
 add('裂纹是否形成贯通分离路径？','以原始存活bulk/interface图及锚可达性判定：'+ '；'.join(f'{name}最终{int(x["final"]["detached_nodes"])}材料点真实游离、{int(x["final"]["components"])}分量' for name,x in q.items())+'。游离块必须没有存活外部材料连接；空间邻近与接触对均不计材料键。不能用完全失效键多或损伤网络合并代替这个判定。')
 add('是否出现20–200节点游离块？','；'.join(f'{name}: 历史最大目标带数量{x["history_max_N20_200"]}，最终{int(x["final"]["N20_200"])}，历史最大游离分量{x["history_max_free_nodes"]}点' for name,x in q.items())+'。≥34且≤200、≥114及内部连接均单独记录。20节点等效参考体积直径约168.4 µm，仅约3.37dx；达到节点数本身不代表充分解析的实体或医学安全。')
 add('是否出现至少5个并通过全部目标？',f'五碎块完整探索目标={status["five_fragment_target_met"]}。目标预先冻结，另需脱离质量≥10%、最大块≤脱离质量50%、1–2点细屑≤20%、内部存活连接及3.4→3.6 ms同一原始成员连通保持。该观察窗开始的剩余牵引为峰值的{protocol["cohort_remaining_load_fraction"]:.5g}，不是严格零流体力；寿命是有限窗右删失观察。')
 add('是否仍以1–2节点碎屑为主？','；'.join(f'{name}: 脱离质量中的细屑比例={fmt(x["final_metrics"]["fine_debris_mass_fraction_detached"],100)}%，最大块占脱离质量={fmt(x["final_metrics"]["largest_free_fragment_fraction_detached"],100)}%' for name,x in q.items())+'。无脱离时分母为零，比例保留null，不把它伪写为0%成功。')
 add('还有多少三维存活连接阻止候选区域脱离？','；'.join(f'{name}最终：'+', '.join(f'{r["region"]}（{r["candidate_nodes"]}节点）{r["surviving_crossing_bonds"]}条存活跨界键、锚可达{r["anchor_reachable_nodes"]}点' for r in bridge['cases'][name]['final']['regions']) for name in q)+'。这些历史诊断体积不是预设碎块或新增切割。完整原始bridge ID、长键、压缩/未达到正伸长阈值和偏离中剖面的锚路径均已保存；图最小割是可达性诊断，不是删除建议或断裂能。')
 add('目标剪切需求与当前声流能力差多少？',f'相同几何、模式、时程下，本轮B/A峰值牵引倍数={status["template_traction_ratio"]:.9g}，相同倍数适用于线性模板的载荷冲量；**不能称为声流不足倍数**。A并非微泡预测，真实声流应力/间隙仍未确定；若B未达到目标，也不能把1000 Pa当成成功阈值。外功来自实际变形，不能按载荷倍数线性推断。机械目标需求是否识别={status["SHEAR_REQUIREMENT_IDENTIFIED"]}。')
 add('主要限制是哪一项？','A能够对这一规定流体模板给出真实材料响应，但不足以判断所有微泡机制。B检验空间反向剪切的机械作用，仍受冻结参考面加载的限制。'+f'B首次保守几何诊断为{s["inverse"]["frozen_geometry_first_diagnostic"]}。该阈值为位移>0.1gap或法向转角>5°，仅提示模型适用性，未作为伪物理停止条件；其后仍是明确的死牵引探索。B在0.37725ms触发诊断，早于0.88925ms首次软化，其末态约99.65%外功成为动能；固定原始ID牵引继续加速已脱离表面节点，不能据此认证真实流体能量输入或裂解效率。末态部分游离点沿x移动数十毫米，仅为不含已求解流体阻力的死载机械结果，并非血管输运预测。A未触发几何阈值也不代表运动壁面反馈或流体阻尼已验证。V3映射局部误差、三维桥接和未校准材料共同限制推断，不能仅由两次计算唯一归因于本构缺少混合模态。')
 add('是否值得重新启用BoTorch做微泡参数反向优化？','目前不宜作为物理闭环优化。旧Optuna/BoTorch软件与独立环境可复用，但缺少声学参数→近壁流体牵引的可计算可信关系，本轮0个优化trial。先完成坐标/间隙/边界输入、微泡一阶振荡与二阶平均流验证，再以真实目标带数量和质量分层约束优化；κ和软化数只作DAMAGE_PRECURSOR_DIAGNOSTIC。')
 add('是否具备真实V3声流断裂条件？','否。虽然实际V3面片、原始ID和native接口可复用，当前尚无其几何上的可信微泡流体解，局部映射与大变形表面更新也未通过物理验证。本轮未完整运行1927点V3、未做双向FSI。下一轮应首先补齐近壁声流与表面传递条件，不能继续以压头或人工删键替代研究主线。本轮更接近目标的是可信剪切耦合链条和对一个低维模式的真实证据；仅两个幅值不能确定最优时空结构，也不能将其变化唯一归因于方向组织。')
 intro=f'''# CLOT_ACOUSTIC_STREAMING_SHEAR_MULTIFRAGMENT_V1

Run: `{RUN}`。软件开发与两组原生3.6 ms计算完成；实际微泡Streaming仍未确定。完整探索目标达成：**{status['five_fragment_target_met']}**。

当前有效入口：`analysis/{TASK}/`。生产native库未改动，SHA `{sha(LIB)}`；1792材料点+336锚点，原始95666键（2736原预切），G=940 Pa、K=100000/3 Pa，本体s0=.15/sc=.18、界面.140625/.16875保持原值，均未为碎块目标调参。

本轮停止迭代历史刚性压头路线：单/双压头属于规定运动的机械法向接触问题，可保留为独立证据，不能解释为声流。下面所有断裂图均来自原生接受态、位移×1。

| 工况 | 面片峰值剪切(Pa) | 最终游离点 | 历史最大块 | 历史最大20–200块数 | 外功(nJ) | 最大能量残差(%) |
|---|---:|---:|---:|---:|---:|---:|
'''
 for name,x in q.items():intro+=f'|{name}|{x["peak_shear_Pa"]:.8g}|{int(x["final"]["detached_nodes"])}|{x["history_max_free_nodes"]}|{x["history_max_N20_200"]}|{x["final"]["W_J"]*1e9:.8g}|{x["max_relative_energy_residual"]*100:.6g}|\n'
 tail=f'''\n## 数值证据和资源

每个接受步检查真实原始ID、κ/d不可逆、Δd≤0.1、原生失效键移除和能量余额。dt=0.25 µs，损伤事件最低0.03125 µs；平滑初始段采用此前验证过的更小步长。W=ΔU+ΔK+D+ΔU_pair−J_pair+R；无压头能、无数值阻尼；R不归入损伤耗能。外功以实际位移路径和时变剪切做Gauss积分。完整重放验证见 DELIVERY_VALIDATION.json。标量history中的net_moment_z_Nm和moment_impulse_z_Nms为参考坐标的名义力矩/积分；关键FP64接受态的实际x(t)×F(t)力矩另存CURRENT_MOMENT_KEY_STATES.csv。实际线冲量及外功均按真实外力和位移计算；没有从稀疏位置样本伪造全时程实际角冲量。

两例合计墙钟{resources['native_wall_s']:.3f} s、CPU {resources['native_cpu_s']:.3f} s、峰值RSS {resources['peak_RSS_bytes']/2**20:.3f} MiB；每例1 rank×1线程，无GPU、无新增墙钟/RSS上限。小量流体解析模态计算和只读映射/拓扑分析另计。

## 出处与证据边界

- [APECSS作者软件说明](https://joss.theoj.org/papers/10.21105/joss.05435)：核对其气泡动力学/声辐射软件范围，没有将其当成已经实现近壁声流。
- [Muller与Bruus时变声流研究摘要](https://arxiv.org/abs/1509.02554)：核对一阶/二阶与声流建立衰减的区分；本轮读取摘要，不声称已阅读全文或复现其边界模型。
- [LAMMPS fix external](https://docs.lammps.org/fix_external.html)：原生外力回调接口。
- Stokes横向模态、sinh稳态解与正弦系数由本轮给定边值问题直接推导，并经独立有限差分验证；不是从缺失论文编造的微泡应力公式。

`ACOUSTIC_AND_BUBBLE_INPUT_AUDIT.json`记录实际文件、HDF5字段和代码行。启动输入审计在可选优化器版本查询时发现其不在native环境，已修正为明确记录独立环境，不重算物理工况。旧封存结果均保留。本轮不自动追加第三例、完整V3或自动优化。
'''
 text=intro+'\n'+'\n'.join(answers)+tail;path=O/'REVIEW_CLOT_ACOUSTIC_STREAMING_SHEAR_MULTIFRAGMENT_V1.md';path.write_text(text)
 figures=read(O/'FIGURE_MANIFEST.json')['figures'];items=''.join(f'<section><h2>{html.escape(f["name"])}</h2><p>{html.escape(f["caption"])}</p><a href="{f["name"]}.svg"><img src="{f["name"]}.png"></a></section>' for f in figures)
 page='<!doctype html><html><meta charset="utf-8"><title>Streaming shear research</title><style>body{font:18px system-ui;max-width:1500px;margin:30px auto;padding:20px;color:#25303c;background:#f4f6f8}section{background:white;padding:20px;margin:25px 0}img{max-width:100%}h1{font-size:30px}p{line-height:1.5}</style><h1>声流／剪切多碎块研究</h1><p>真实微泡Streaming：尚未确定。两组规定剪切原生计算完成。所有断裂图为实际接受态，位移×1。</p>'+f'<p><a href="../../../verification/{TASK}/{RUN}/{path.name}">20问中文报告</a></p>'+items+'</html>'
 (F/'OPEN_RESULTS.html').write_text(page);(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta charset="utf-8"><a href="../../../results/{TASK}/{RUN}/OPEN_RESULTS.html">打开本轮全部科研图表</a>')
 print('REPORT_DONE',len(answers))
if __name__=='__main__':main()
