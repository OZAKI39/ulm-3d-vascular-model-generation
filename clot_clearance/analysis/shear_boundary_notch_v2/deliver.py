"""Build task-scoped scientific review and static index from verified observations."""
from pathlib import Path
import sys,json,csv,os,html,hashlib,time
base=Path(sys.argv[1]);figs=Path(sys.argv[2]);root=Path.cwd();run=base.name
J=lambda p:json.loads(p.read_text())
s=J(base/'SHEAR_BOUNDARY_AUDIT.json');n=J(base/'NOTCH_FAILURE_DIAGNOSTIC.json');checks=J(base/'IMPLEMENTATION_CHECKS.json');rows=list(csv.DictReader((base/'CASE_COMPARISON.csv').open()));row={a['case']:a for a in rows};d=J(base/'DYNAMICS_AND_ENERGY_DETAIL.json');start=J(base/'provenance/SERVER_START.json');end=J(base/'provenance/SERVER_END.json')
assert all(v['finite'] and v['actual_time_index'] and v['native_force_gradient_match'] and v['finite_difference_match'] and v['energy_residual_below5percent'] for v in checks.values())
a=s['baseline'];b=s['periodic'];m=n['baseline'];q=n['improved'];early=n['early_before_forced_collar'];frel=Path(os.path.relpath(figs,base));report='REVIEW_LAMMPS_SHEAR_BOUNDARY_AND_NOTCH_LOCALIZATION_V2.md'
status=dict(task='LAMMPS_SHEAR_BOUNDARY_AND_NOTCH_LOCALIZATION_V2',run_id=run,execution='COMPLETED_WITH_SCIENTIFIC_LIMITATIONS',IMPLEMENTATION_VALIDITY={'status':'VALIDATED_FOR_TWO_DECLARED_CASES','evidence':'PREFLIGHT.json; IMPLEMENTATION_CHECKS.json; source-consistent native force and finite-difference checks'},BOUNDARY_CONDITION_VALIDITY={'M2P':'PERIODIC_GRAPH_AND_SELECTED_DOF_INTEGRATION_VALIDATED','M5B':'DECLARED_BOUNDARY_EXECUTED; EARLY_INITIATION_INFORMATIVE; LATE_COLLAR_STRETCHING_UNSUITABLE_FOR_SINGLE_CRACK_FRACTURE_WORK','forced_collar_threshold_s':n['kinematic_protocol_limit']['threshold_time_s']},SHEAR_MODULUS_MAPPING={'status':'BOUNDARY_SENSITIVITY_EXPLAINED; DISCRETE_RESIDUAL_REMAINS; NOT_FULL_CALIBRATION','M2_ROI_Pa':a['G_ROI_Pa'],'M2P_ROI_Pa':b['G_ROI_Pa'],'M2P_nominal_Pa':b['G_nominal_Pa'],'input_G_Pa':940,'M2P_error_fraction':b['relative_error_to940'],'discrete_prediction_Pa':s['discrete_moment']['predicted_full_neighborhood_small_strain_G_Pa']},NOTCH_INITIATION_LOCALIZATION={'status':'IMPROVED_TO_ACTUAL_NOTCH_FRONT','first_time_s':q['first_time_s'],'peak_target_band_fraction':q['snapshots']['peak']['target_band_fraction'],'gates':q['notch_dominance_checks']},CRACK_PATH_IDENTIFIABILITY={'status':'EARLY_LOCALIZATION_IMPROVED; CONTINUOUS_NEW_SINGLE_SURFACE_UNRESOLVED','early_connected_components':early['component_count'],'end_off_plane_fraction':q['snapshots']['end']['off_projected_plane_fraction']},ENERGY_LEDGER_VALIDITY={'status':'SOURCE_CONSISTENT_INDEPENDENT_LEDGER_WITH_SIGNED_TOPOLOGY','M2P_residual_work_fraction':checks['M2P']['energy_residual'],'M5B_residual_work_fraction':checks['M5B']['energy_residual'],'negative_topology_events':int(row['M5B']['negative_topology_events'])},FRACTURE_WORK_IDENTIFIABILITY={'status':'GC_UNAVAILABLE','Gc_J_m2':None,'apparent_fracture_work_candidate_J_m2':None,'reason':'No unique new physical crack area; early connectivity proxy zero, late distributed and prescribed-collar damage. Signed topology transfer is not all irreversible material dissipation.'},DISCRETIZATION_LIMITATION={'status':'NOT_CONVERGED','evidence':'Reused V1 D1 dx50→40um stress-curve difference44.15%, net topology-energy difference48.13%; fixed-stencil elastic quadrature differs from isotropic continuum.'},PRODUCTION_PARAMETERS_UNCHANGED={'status':'VERIFIED','unchanged_hash_count':end['checked'],'unexpected_changes':end['unexpected_changes'],'laboratory_change':'Only optional observer and two declared boundary cases in macro_lab_v1.py; native solver, independent energy and Y production modules unchanged'},new_nontrivial_LAMMPS_cases=2,zero_step_preflights=2,implementation_bug_retest_cases=0,AI_use={'allowed':'Traceable numerical boundary-sensitivity and initiation-location examples tagged with exact case, discretization and validity limits','not_certified':['Intrinsic clot strength or Gc','Mesh-independent fragment-size labels','Patient-specific response','F7 microbubble cluster/force labels','AI material calibration ground truth']},next_priority='Design and preflight a grip loading protocol that preserves notch initiation without necessarily straining prescribed collar bonds past s00; only after this boundary issue is resolved assess a traceable single crack surface.',next_task_started=False,WSL_document={'path':'/home/lzy/projects/clot_rupture.md','receipt':'provenance/WSL_DOCUMENT_UPDATE.json','state':J(base/'provenance/WSL_DOCUMENT_UPDATE.json')['state'],'after_sha256':J(base/'provenance/WSL_DOCUMENT_UPDATE.json')['after_sha256']},raw_server_root=str(root/base),figures=str(root/figs))
(base/'FINAL_STATUS.json').write_text(json.dumps(status,indent=2,ensure_ascii=False))
raw_files=[p for p in base.rglob('*') if p.is_file()];fig_files=[p for p in figs.rglob('*') if p.is_file()];heavy=[p for p in raw_files if p.name in ['body.data','history.csv','native_energy_density.dat']]
perf=dict(measured_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),cpu_affinity_count=start['affinity'],cgroup_cpu_equivalent=float(start['quota'])/float(start['period']),gpu_start=start['gpu'],gpu_solver_used=False,installed_dependencies=[],execution_model='Two independent processes concurrently; one MPI rank each; OMP/BLAS threads1; existing CPU native library. GPU unused; no pressure test.',cases={k:J(base/k/'EXECUTION.json') for k in ['M2P','M5B']},process_wall=J(base/'CAMPAIGN_EXECUTION.json'),parallel_campaign_seconds=max(x['wall_seconds'] for x in J(base/'CAMPAIGN_EXECUTION.json')),verification_bytes=sum(p.stat().st_size for p in raw_files),figure_bytes=sum(p.stat().st_size for p in fig_files),server_only_heavy_bytes=sum(p.stat().st_size for p in heavy),storage_measurement_scope='Task directories at report generation; final manifest carries exact individual sizes. Uncompressed logical file bytes.',supervisor_end=end['own_supervisor'],resources_start='provenance/SERVER_START.json',resources_end='provenance/SERVER_END.json')
(base/'PERFORMANCE_AND_STORAGE.json').write_text(json.dumps(perf,indent=2))
# Do not edit accepted case input metadata; clarify inherited generic descriptors here.
(base/'INPUT_INTERPRETATION.json').write_text(json.dumps(dict(M2P='No notch or anchor/interface bonds; inherited generic plan descriptions refer to available laboratory builders, not active geometry. Actual removed_notch_bonds=0, anchors=false; topology in geometry.npz is authoritative.',M5B='Pre-cut2570 bonds; no anchor/interface test; all132610 active bonds are bulk. Same geometry and masks as M5; only prescribed axial collar profile changes.',actual_time_and_state='Exact accepted double-precision snapshots and native event steps. Legacy M5 lacks the immediately pre-first-break full field.'),indent=2))
text=f'''# 剪切边界与缺口起裂定位审查 V2

Run ID：`{run}`。本轮完成两个正式算例 M2P、M5B；两个 `run 0` 只用于几何/积分组检查，无第三个正式复验。结果是**剪切边界影响得到解释、缺口起裂得到改善，但单裂面与 Gc 仍不能识别**。生产参数保持不变。

[结果入口（静态图索引）](OPEN_RESULTS.html) · [分项科研状态](FINAL_STATUS.json) · [四算例汇总](CASE_COMPARISON.csv) · [实施检查](IMPLEMENTATION_CHECKS.json)

## 1. 身份、范围与冻结条件

服务器实际目录为 `/workspace/y_clot/project`，是平铺血栓项目；提示中的 `clot_clearance/` 不是服务器额外一级目录。环境为 `/workspace/y_clot/env.sh`。服务器实际分支 `{end['git_branch']}`、HEAD `{end['git_head']}`；WSL GitHub 快照的分支与提交见 `provenance/WSL_START.json`，没有将它覆盖服务器。旧结果来自 `verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z`，没有重跑或改写。提示中的 GitHub README/SOURCE 在 WSL 快照读取，缺失路径及实际入口见 `provenance/PROJECT_ENTRIES.json`。

锁定：G=940 Pa、K=33333.333333 Pa、ρ=1000 kg/m³、bulk s00=0.15、interface s00=0.140624859375、α=0、dx=50 µm、δ=150.5 µm、接触距离因子1、数值阻尼500 s⁻¹、dt=0.5 µs。未改变 native LPS、力密度乘节点体积的 SI 转换、独立能量表达式、密度或阻尼。`provenance/REUSED_INPUT_SHA256.json` 的85项身份中，84项在结束时逐项一致；唯一允许变化的是实验室适配器，旧版与 diff 已保留。正式执行版本见 `provenance/launch_*.py` 与 `LAUNCH_SOURCE_SHA256.json`；后续仅补充只读分析元数据与图示，不改运行版本。

两试样均为600×600×1000 µm、2880节点、三层150 µm夹持区，夹具均值距离H=850 µm、2016个内部节点完全自由。M2/M2P只规定夹具x位移及z位置，y自由；M5/M5B只规定夹具z位移，x/y自由。`fix move` 的自由分量与内部 `nve` 组不重叠。M2P保留原6.8 µm幅度、10 ms加载+10 ms卸载；M5B保留680 µm幅度、24 ms余弦加载、相同缺口和原132610键。

唯一主改动：M2P将x/y外侧改为600 µm固定盒的真实周期边界；M5B将夹持层原刚性位移±disp/2改为 `u_z=disp(t) Z/H`，平均夹具开度不变、夹持层内部允许被规定的轴向应变。没有同时改变尺寸、时程、阈值或接触模型。方案与区域在分析/运行前冻结，见 [ANALYSIS_PLAN](ANALYSIS_PLAN.json)、[CASE_REGISTER](CASE_REGISTER.json)、[先行只读诊断](READ_ONLY_DIAGNOSIS.md)。

## 2. M2为何偏低：反力公式正确，边界场与测量口径影响显著

|同一20%–70%加载拟合区|M2自由侧面|M2P周期侧面|
|---|---:|---:|
|名义剪切牵引 / 实测ROI对称剪切，Pa|{a['G_ROI_Pa']:.6f}|{b['G_ROI_Pa']:.6f}|
|名义剪切牵引 / 夹具开度H，Pa|{a['G_nominal_Pa']:.6f}|{b['G_nominal_Pa']:.6f}|
|ROI斜率相对输入940 Pa|{a['relative_error_to940']*100:.4f}%|{b['relative_error_to940']*100:.4f}%|
|峰值ROI内局部剪切CV|{a['states']['peak']['ROI_local_gamma_CV']*100:.4f}%|{b['states']['peak']['ROI_local_gamma_CV']*100:.4f}%|
|拟合区最大自由DOF K/U|{a['fit_max_KE_U']*100:.4f}%|{b['fit_max_KE_U']*100:.4f}%|
|拟合区末累计数值阻尼能 / 当时U|{a['fit_end_Dnum_over_U']*100:.4f}%|{b['fit_end_Dnum_over_U']*100:.4f}%|

这里ROI量为从真实参考/当前节点坐标最小二乘得到 `γ_xz=F_xz+F_zx`；不是把 `F_xz` 单独当γ。名义应力为原始参考节点单元面面积 `A0=600×600 µm²` 上的夹具合反力。严格地说，759和992仍是**名义牵引对ROI应变**的斜率；M2P均匀程度明显改善，不能据此宣称局部非局部应力具有唯一分配方式。

M2峰值时 F_xz=0.0105838791、F_zx=−0.00235793756、γ=0.00822594151，而名义γ=0.00799995850；存在明显转动及高度方向非均匀场。将应变换成名义口径只会从759升到781 Pa，不能消除19%差异。中央逐节点γ的CV为7.55%；M2P降到0.0177%。图A的剖面保留两个梯度分量，没有挑选另一个更好看的拟合窗口。

M2峰值上夹具外反力2.257080141 µN，下夹具反向，合力相对上反力约2.95×10⁻¹⁴；独立力重建与原CSV相差约1.7×10⁻²¹ N。当前面拟合面积3.599996543×10⁻⁷ m²，参考面积3.6×10⁻⁷ m²，面积变化不足以解释刚度差异。不能用仅点中心凸包的550×550 µm²替代节点单元面以人为抬高模量。中心截面牵引约6.298 Pa、顶部6.270 Pa；中心端点分配的参考键维里代理约7.788 Pa，后者是非唯一局部应力代理，不能直接当均匀材料局部G。

M2P周期图经独立整数偏移枚举完全一致：165600条唯一键，其中30420条跨周期边界；最大邻居数122，无缺失/重复/自键。δ及接触距离均小于半盒；参考长度和运行长度均取一致minimum image；位移只依赖z，因此x/y周期平移与加载兼容。没有将普通自由面旧键表直接用于周期盒。见 [PREFLIGHT](PREFLIGHT.json)。

M2P的ROI模量992.183 Pa比940高5.551%，低于15%工程筛查阈值；名义开度口径1106.782 Pa仍明显不同，说明夹具和测量高度不能忽略。122邻居立方点阵的冻结权重四阶剪切矩为0.0706402020，而连续各向同性值1/15=0.0666666667。由**原独立储能的解析小应变展开**得：

`G_discrete = 15 G × Σ(h rx² rz²/r³) / Σ(h r) = 996.026849 Pa`。

±10⁻⁵能量二阶差分一致；992.183与该预测相差−0.386%。这是一项只读离散积分诊断，**不是新增全部节点仿射运动的模拟，也不是收敛研究**。残余差异与有限夹具、有限小应变、惯性仍相容；本轮没有分别扫描尺寸、horizon、加载速率和阻尼，因此不能给自由表面截断、夹持干扰各自分摊百分比。

M1已有E=2681.349 Pa、ν=0.48314428，按各向同性关系推得G={s['M1_inferred_G_Pa']:.6f} Pa；M2相差{s['M1_inferred_vs_M2_difference']*100:.3f}%，M2P相差{s['M1_inferred_vs_M2P_difference']*100:.3f}%。不同边界试样不构成同一套已标定各向同性常数；M3P受限压缩接近目标也不消除该限制。

[图A1：应力—应变与高度剖面]({frel}/A1_shear_definitions_and_profiles.png) · [图A2：ROI真实场与惯性]({frel}/A2_shear_ROI_and_inertia.png) · [机器可读剪切审查](SHEAR_BOUNDARY_AUDIT.json)

## 3. M5原始夹具起裂的拓扑证据

原预缺口是z=0、x<−100 µm的穿厚度断面：以每条跨面键的实际参考交点判断并移除2570条键，逐条重新构造后与原M5图一致，未发现不应存在的跨预切口存活桥键。缺口前缘到夹持过渡面仅350 µm=2.326δ，尚不能视为远场加载的无限大试样。

4.625 ms首断两键索引130026、132604，原子ID2657–2658、2877–2878；参考中点(275,±275,350) µm，是夹持过渡角点。端点只有46/45邻居，最大完整邻域122；加权体积仅完整值的0.367/0.360。在**首断接受状态**，夹具最大利用率1.000050，缺口0.869405。原M5没有前一接受步完整坐标，故无法追补精确“首断前”场；此处没有把首断状态伪称前一步。

刚性夹持层到自由内部的位移梯度变化，叠加角点缺邻居，形成了比缺口更高的风险；本轮受控边界变化支持此解释，但没有通过独立表面修正试验定量隔离每项原因。最早≥10个六连通dx体素的断键簇在5.1765 ms形成，位于夹具区域；5.488 ms峰值前后缺口随后扩展，不能将全部后期87.15%离面断键都说成夹具损伤。

空间区域固定为：缺口前缘线的径向距离≤δ优先；其余距z=±350 µm夹持过渡面≤δ；其余距外表面≤δ；剩余为自由内部。目标带另定义为 `|mid_z|≤δ/2 且 mid_x≥−100µm−δ/2`。目标带比例与是否跨z=0是不同指标，不能相互替代。

[图B1：三维原生首断与利用率]({frel}/B1_notch_vs_grip_first_break.png) · [图B3：清晰正投影与真实缺口]({frel}/B3_notch_and_grip_projection.png) · [图B2：全事件时空定位]({frel}/B2_native_event_localization.png)

## 4. M5B：早期改善，后期边界本身仍不适合单裂纹工作

|指标|M5|M5B|
|---|---:|---:|
|首断时间 / ms|4.6250|5.2300|
|首断位置|夹持角点|实际缺口前缘(−100,±225,0) µm|
|名义峰值牵引 / Pa|303.227412|306.494994|
|峰值时间 / ms|5.4880|6.1860|
|峰值累计断键|292|704|
|峰值目标带比例|69.863%|78.977%|
|峰值夹持区断键比例|17.123%|0.284%|
|峰值K/U|0.9361%|0.9367%|
|末态断键|8918|31340|
|末态不跨z=0比例|87.1496%|95.4116%|

M5B保存了真实断前一步10459（5.2295 ms）：缺口利用率0.999899、夹持区0.712026；第10460步缺口两键达到1.000093并由原生机制断开。5.6685 ms首个连贯10体素断键簇全部在缺口区，触及预定前缘条带。三个早期定位门槛（首断/连贯簇、峰值目标带≥70%、夹具≤20%）均达到，**完整单裂面和面积门槛仍未达到**。名义峰值只是该边界下的试样响应，不是材料强度标定。

M5B采用的线性夹持层位移还有一个明确限制：夹持层内两端轴向均受规定的竖直键，其伸长至少等于名义轴向应变。加载超过15%后必然触及s00；按固定时程，该时刻为6.842375 ms。这一推论在查看新算例结果前写入 [KINEMATIC_VALIDITY_LIMIT](KINEMATIC_VALIDITY_LIMIT.json)，没有事后改停止时间或阈值。本方案改善了早期过渡，但**不能作为完整24 ms高应变单裂纹边界方案**。

同夹持层双端受轴向规定的26308条初始键中，最终13788条断裂，其中824条在该解析阈值前已断、12964条在其后断；不是所有夹持键都在一个时间整齐断裂，也不能把全部后期失效都归于该下界。末态夹持过渡区域24206/31340断键，区域与“键两端均受规定”集合不同，已分别记录。自由裂纹演化与边界强制夹持层损伤混合，不能使用末态漂亮外观认证裂面。未将这一方案局限冒称程序bug来追加第三次计算。

[图C1：反力、竞争和定位]({frel}/C1_loading_and_localization_comparison.png) · [图C2：同坐标真实首断/峰值/末态]({frel}/C2_actual_fracture_states.png) · [完整断裂诊断](NOTCH_FAILURE_DIAGNOSTIC.json)

## 5. 新增面积、独立能量与Gc

提前定义的1 ms接受态中，选取严格早于6.842375 ms的最后一个，即6 ms。此时368条断键，84.239%在目标带、夹持区0、26.630%不跨z=0；整个点阵仍为一个2880节点连通分量。沿原缺口前缘逐50 µm厚度条带检查跨面存活键，尚无从前缘连续推进的完全失连投影单元，代理面积为0。这是连接性判据的分辨率限制，不能解释成物理上“裂纹面积恰好为零”。

末态仅两个边缘条带出现一个50 µm步长的完全失连，投影连接损失代理为5×10⁻⁹ m²；其余仍有跨面桥接且大量离面/夹持区损伤。名义50 µm格点面积分辨率为2.5×10⁻⁹ m²，单格已相当于此代理的一半；空间分辨率和表面唯一性均不足。没有用断键数乘dx²，没有将投影代理当真实新自由面，也没有用平滑外壳补造裂面。见 [EARLY_AREA_ATTEMPT](analysis/M5B/EARLY_AREA_ATTEMPT.json)、[AREA_ATTEMPT](analysis/M5B/AREA_ATTEMPT.json)、[图D2]({frel}/D2_crossing_connectivity_limits.png)。

复用V1锁定源码一致的独立势能：每一步实际规定DOF外功、当前bulk储能、保守接触储能、自由DOF动能、数值阻尼、固定构型断前/断后拓扑储能差及独立余额。原native PE仍只作诊断，不当完整LPS储能；负拓扑增量原样保留，残差没有塞入断裂项。

M5B：净带符号拓扑项17.655581 nJ，正增量和17.945217 nJ、负增量和−0.289636 nJ，共714次负增量；最大独立余额/最大外功为0.060807%（M5为0.038378%）。M2P为1.6490×10⁻¹¹的无量纲比例。每个新算例都通过3个实际状态的native力对独立能量梯度与位移有限差分检查。小余额表明账本相容，**不把带符号拓扑项认证为不可逆材料耗散**。

M5B末外功为{d['M5B']['energy_end']['W_ext_J']*1e9:.6f} nJ，末数值阻尼为{d['M5B']['energy_end']['D_num_J']*1e9:.6f} nJ；数值阻尼不是血栓物理黏弹耗散。主要早期峰值K/U约0.94%，但断裂开始后在7.9655 ms达到22.286%（K=0.69847 nJ、U=3.13407 nJ），不能把整段断裂视为准静态。启动时0.093 ms还出现K/U=206.94%，此时不能用接近零的储能作全程准静态证据；全时程及峰值/末态分账见 [DYNAMICS_AND_ENERGY_DETAIL](DYNAMICS_AND_ENERGY_DETAIL.json) 和 [图D1]({frel}/D_independent_energy_and_limits.png)。在可变形夹持层，功共轭广义载荷 `Q=−Σ f_material·(∂u_prescribed/∂disp)` 不等于顶部合反力，图C1同时保留两者；外功逐DOF计算，没有用顶部力—均值开度面积替代。

因此 **GC_UNAVAILABLE**；也未生成 `SPECIMEN_DEPENDENT_APPARENT_FRACTURE_WORK_CANDIDATE`。即使未来面积和区间能识别，该候选也只能描述特定试样、dx、δ、边界、加载速率，不能自动等同固有Gc。

## 6. 模型适用性与AI范围

当前LPS适合在明示离散和边界的前提下研究小应变弹性响应及原生断键位置的数值趋势；本轮发现边界和立方积分方向性足以显著改变表观G，而早期非目标起裂可通过边界改善。它尚未获得本体强度、纯界面强度或材料Gc标定。现有证据不足以把非目标破坏单独归因于“材料太脆”，也不足以认定必须增加软化；本轮不升级本构。

复用V1限制：M5 dt/2整条曲线差约4.78%；dx50→40 µm曲线差约44.15%、净拓扑能量差约48.13%，没有离散收敛。本轮M2P与M5B没有各自再做dt/2验证。F7微泡内部构型仍未时间收敛，本轮未重算。

允许的AI用途仅为带算例身份与局限标签的数值边界敏感性、起裂定位和失败模式示例；不可认证真实材料强度/Gc、网格无关碎块尺寸、生物实验真值、微泡群集载荷或患者预测标签。科研状态分为9项记录在FINAL_STATUS，没有一个总PASS替代科学结论。

## 7. 资源、实现差异与复现

实际cgroup为{perf['cgroup_cpu_equivalent']:.2f}核，CPU affinity {perf['cpu_affinity_count']}；开始可用RAM约52.79 GB、磁盘296.22 GB，RTX4090约24 GB显存空闲。两个独立单rank CPU算例并发，OMP/BLAS各1线程，没有GPU内核、压力测试或依赖升级。M2P 40000步、wall {perf['cases']['M2P']['wall_seconds']:.3f}s；M5B 48000步、wall {perf['cases']['M5B']['wall_seconds']:.3f}s；包含进程启动的并行时长约{perf['parallel_campaign_seconds']:.3f}s。此数字是求解阶段墙钟，不是把全部审查时间伪报成求解时间。任务supervisor已EXITED。

报告生成时verification约{perf['verification_bytes']/1e6:.2f} MB，图约{perf['figure_bytes']/1e6:.2f} MB；其中{perf['server_only_heavy_bytes']/1e6:.2f} MB完整逐步历史/原始data/native能量仅留服务器。最终逐文件字节及SHA见FILE_MANIFEST，WSL精简下载见 `provenance/COMPACT_DOWNLOAD.json`。没有制作4K视频。

只增补原 `macro_lab_v1.py` 的可选观察器和两种边界入口，默认V1/生产调用不引入这两种新边界。真实输入、geometry、日志、原生事件和接受态都保存在M2P/M5B目录；变更前源码、逐行diff及正式启动源码均保留。源文件身份和历史内容再次核验见SERVER_END；生产Y模块、native库、独立energy库与旧数据一致。

只读复现（在已有服务器环境内）：

```bash
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
python analysis/shear_boundary_notch_v2/finish.py verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z {base} {figs}
python analysis/shear_boundary_notch_v2/supplement.py {base} verification/lammps_material_macro_calibration_and_timestep_audit_v1/20261010T215825Z {figs}
```

正式运行命令已逐项记录在 `CAMPAIGN_EXECUTION.json` 和 `run_campaign.py/sh`；本轮完成后不再执行这些求解命令。复现求解应使用新的独立输出目录，不能覆盖本轮接受态。

图源：剪切场在 `baseline/M2/*_fields.npz`、`analysis/M2P/*_fields.npz`；原生事件空间映射在两个NOTCH分析目录的 `EVENT_LOCALIZATION.npz`；B图利用率来自真实first-break风险文件；C图使用实际首断/峰值/末态；D2面积图源为 `analysis/M5B/AREA_GRID_SOURCE.npz`。轻量曲线每20个接受步加所有拓扑事件、精确峰值和末态保存，没有插值出破坏状态。完整CSV在服务器可逐步复算。

## 8. WSL更新、资料入口与唯一下一步

实际WSL文档 `/home/lzy/projects/clot_rupture.md` 的本轮修改只覆盖④弹性边界证据、⑤缺口与Gc状态及末尾下一步/新结果链接。备份 `provenance/clot_rupture.before.md`、逐行diff `provenance/clot_rupture.changes.diff`、写入前并发SHA核查及写入后SHA在 `provenance/WSL_DOCUMENT_UPDATE.json`。以该回执核验实际写入结果，不把服务器同名文件视为WSL更新。既有AGENTS已经有后续更新惯例，本轮不重复追加。

服务器原始数据绝对路径：`{root/base}`。本地关键结果通过显式文件列表及SHA回收，完整原始轨迹仍在服务器；没有Git推送或删除历史成果。

**下一轮唯一优先项**：设计并先行审查一种不会必然把受规定夹持层键拉过s00、同时保留缺口起裂优势的夹持传力边界，再检验是否形成可追踪的新增单裂面。本轮不继续试错，也不自动开始软化、材料库、含栓CFD或AI扫描。

参考实现以锁定源码为准；公开说明可见 [LAMMPS bpm/peri](https://docs.lammps.org/bond_bpm_peri.html)、[fix move](https://docs.lammps.org/fix_move.html)、[boundary](https://docs.lammps.org/boundary.html)。这些链接用于边界命令语义，不替代本项目native修改身份。
'''
(base/report).write_text(text)
# Preserve readable static scientific gallery; links are relative and work on WSL/server.
fig_names={'A1_shear_definitions_and_profiles':'A1 剪切曲线、测量口径与剖面','A2_shear_ROI_and_inertia':'A2 ROI与实际剪切场','B1_notch_vs_grip_first_break':'B1 三维原生首断与利用率','B3_notch_and_grip_projection':'B3 正投影：缺口与夹具','B2_native_event_localization':'B2 原生断键时空分布','C1_loading_and_localization_comparison':'C1 反力与定位比例','C2_actual_fracture_states':'C2 真实首断、峰值与末态','D_independent_energy_and_limits':'D1 独立能量账本','D2_crossing_connectivity_limits':'D2 跨面桥接与面积限制'}
header='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>剪切边界与缺口定位 V2</title><style>body{font:17px/1.65 system-ui,sans-serif;background:#eef2f5;color:#173042;margin:0}main{max-width:1180px;margin:auto;padding:30px}section{background:white;border-radius:10px;padding:22px;margin:20px 0}a{color:#146782}img{width:100%;height:auto}code{overflow-wrap:anywhere}table{border-collapse:collapse}td,th{padding:8px 20px;border-bottom:1px solid #d3dde3}</style><main>'''
body=f'<h1>剪切边界与缺口起裂定位 V2</h1><p>Run {run} · 2个正式算例 · 生产参数不变</p><section><b>G_ROI：759.18 → 992.18 Pa；缺口起裂改善；Gc仍不可识别。</b><p>M5B后期夹持层被规定拉伸，不能作为完整单裂纹断裂能实验。图中全部为真实接受态或原生事件。</p><p><a href="{report}">中文完整审查</a> · <a href="FINAL_STATUS.json">分项状态</a> · <a href="CASE_COMPARISON.csv">汇总CSV</a> · <a href="SHEAR_BOUNDARY_AUDIT.json">剪切审查</a> · <a href="NOTCH_FAILURE_DIAGNOSTIC.json">断裂诊断</a> · <a href="PERFORMANCE_AND_STORAGE.json">性能与存储</a></p></section>'
for name,title in fig_names.items():
 href=(frel/(name+'.png')).as_posix();sv=(frel/(name+'.svg')).as_posix();body+=f'<section><h2>{title}</h2><a href="{href}"><img loading="lazy" alt="{title}" src="{href}"></a>'+(f'<a href="{sv}">矢量SVG</a>' if (figs/(name+'.svg')).exists() else '')+'</section>'
body+='<section><h2>来源与限制</h2><p>前后场没有插值生成；原M5缺少首断前一步坐标。完整逐步历史在服务器，轻量状态与关键图已回收到WSL。</p><p><a href="provenance/clot_rupture.changes.diff">WSL文档diff</a> · <a href="provenance/clot_rupture.before.md">修改前备份</a> · <a href="provenance/WSL_DOCUMENT_UPDATE.json">写入回执</a> · <a href="FILE_MANIFEST.json">文件SHA</a></p></section></main></html>'
(base/'OPEN_RESULTS.html').write_text(header+body)
(figs/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><a href="'+os.path.relpath(base/'OPEN_RESULTS.html',figs)+'">打开剪切与缺口审查全部结果</a>')
print(json.dumps(dict(report=str(base/report),figures=len(fig_names),performance=perf),ensure_ascii=False,indent=2))
