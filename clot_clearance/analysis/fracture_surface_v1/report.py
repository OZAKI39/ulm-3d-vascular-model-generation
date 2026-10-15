"""Assemble Chinese evidence report and granular status; no physical computation."""
from common import *
from datetime import datetime,timezone

def load(name):return json.loads((O/(name+'.json')).read_text())
t=load('THEORY_FORMULA_AUDIT'); d=load('DISCRETE_SURFACE_ENERGY_RESULTS'); q=d['cases'][0]; c=load('CROSS_PLANE_REFERENCE_TOPOLOGY_AUDIT')['cases'][0]; e=load('INDEPENDENT_ENERGY_FORM_CROSSCHECK')['cases'][0]; sens=load('GEOMETRIC_SENSITIVITY'); perf=load('INTEGRATION_RESOURCES'); prep=load('PREPARATION_RESOURCES'); plots=load('VISUALIZATION_RESOURCES'); doc=load('WSL_DOCUMENT_UPDATE'); hist=load('HISTORICAL_PROTECTION_AFTER'); checks=load('IMPLEMENTATION_CHECKS'); source=load('SOURCE_LPS_CONSISTENCY_AUDIT'); freeze=load('DESIGN_FREEZE')
assert checks['status']==hist['status']=='PASS' and doc['status']=='PASS_ACTUALLY_UPDATED_AND_READ_BACK'
def record(status,scope,evidence,**kwargs):return dict(status=status,scope=scope,evidence=evidence,**kwargs)
statuses={
 'PAPER_3D_FORMULA_VERIFIED':record('PASS','WSL完整19页PDF及三维公式核查；不是服务器端全文读取','PAPER_EQUATION_AUDIT.md'),
 'INFLUENCE_FUNCTION_COMPATIBILITY':record('PASS','ω与m一致缩放时常数δ消去；不是忽略taper','SOURCE_LPS_CONSISTENCY_AUDIT.json'),
 'FIXED_REFERENCE_M_COMPATIBILITY':record('PASS_INTACT_REFERENCE','完整邻域原m按真实ID对应；不外推为删键后Eq2/Eq9恒等','PRISTINE_REFERENCE_TOPOLOGY_AUDIT.json'),
 'PARTIAL_VOLUME_WEIGHT_COMPATIBILITY':record('PASS_NATIVE_RULE','原h进入全部相关积分；有限点阵不等于连续球','TAPER_SHELL_AUDIT.json'),
 'PRISTINE_TOPOLOGY_IDENTITY':record('PASS','完整909392键，原manifest SHA匹配，原生数据逐ID读回','NATIVE_INPUT_IDENTITY_AUDIT.json'),
 'PERIODIC_REFERENCE_BONDS':record('PASS','131196跨缝键，唯一最小镜像，完整整数邻居枚举一致','PRISTINE_REFERENCE_TOPOLOGY_AUDIT.json'),
 'CROSS_PLANE_CLASSIFICATION':record('PASS','参考z<平面为下侧，等号归上侧；实际无点落平面','CROSS_PLANE_REFERENCE_TOPOLOGY_AUDIT.json'),
 'TWO_SIDED_AREA_NORMALIZATION':record('PASS','两侧各720有贡献节点，A=0.6mm²，无额外2','CROSS_PLANE_REFERENCE_TOPOLOGY_AUDIT.json'),
 'LAMBDA_BOUNDS':record('PASS','未裁剪任何lambda；完整范围在[0,1]','LAYER_PROFILES.json'),
 'DISCRETE_ENERGY_IDENTITY':record('PASS','lambda式与所有邻居的体积+偏能式一致','INDEPENDENT_ENERGY_FORM_CROSSCHECK.json',relative_error=e['relative_global_error']),
 'G0_THEORY':record('ACCEPTED_IDEAL_CONTINUUM_VALUE','论文理想三维连续标称值','THEORY_FORMULA_AUDIT.json',value_J_m2=t['G0_theory_J_m2']),
 'G0_DISCRETE':record('ACCEPTED_NOMINAL_VALUE','实际50µm完整点阵上的理想点状态积分；不是已标定Gc','DISCRETE_SURFACE_ENERGY_RESULTS.json',value_J_m2=q['G0_discrete_J_m2']),
 'DISCRETIZATION_SENSITIVITY':record('PASS_ONE_INTERIOR_TRANSLATION','同网格相位的相邻间隙；没有网格/相位收敛认证','GEOMETRIC_SENSITIVITY.json'),
 'GC_EFFECTIVE_NOT_MEASURED':record('CONFIRMED_NOT_MEASURED','没有自然裂纹推进与分离功/真实新生面积','MODEL_VS_EXPERIMENT_SCOPE.md'),
 'NATIVE_LAMMPS_UNCHANGED':record('PASS','最终双修复源码/库SHA未变；本轮零LAMMPS实例','HISTORICAL_PROTECTION_AFTER.json'),
 'PRODUCTION_LIBRARY_UNCHANGED':record('PASS','生产共享库和env.sh未替换','HISTORICAL_PROTECTION_AFTER.json'),
 'HISTORICAL_RESULTS_PROTECTED':record('PASS',f'{hist["files_checked"]}个历史、源码、库文件SHA不变','HISTORICAL_PROTECTION_AFTER.json'),
 'WSL_DOCUMENT_UPDATE':record('PASS_ACTUALLY_UPDATED_AND_READ_BACK','真实文档已原子写入、读回；旧文本逐字保留','WSL_DOCUMENT_UPDATE.json',after_sha256=doc['after_sha256'])}
assert len(statuses)==18
dump(O/'FINAL_STATUS.json',{'task':C['task'],'run_id':C['run'],'completed_utc':datetime.now(timezone.utc).isoformat(),'outcome':'NOMINAL_THEORY_AND_DISCRETE_IDENTIFICATION_COMPLETE','statuses':statuses,'G0_theory_J_m2':t['G0_theory_J_m2'],'G0_discrete_J_m2':q['G0_discrete_J_m2'],'ratio':q['G0_over_theory'],'GC_EFFECTIVE':None,'GC_EXPERIMENTALLY_CALIBRATED':False,'experiment_comparison':'NOMINAL_ENERGY_SCALE_BELOW_EXPERIMENTAL_REFERENCES','new_native_runs':0,'formal_integrals':2,'next_unique_recommendation':'Controlled Mode-I separation with actual crack area and complete energy accounting; not executed','stopped_after_this_task':True,'delivery_verification':'See WSL_DELIVERY_VERIFICATION.json, independently generated after transfer'})
scope=f'''# 标称能量与实验定义的边界

本轮 G0_THEORY={t['G0_theory_J_m2']:.16g} J/m²，G0_DISCRETE={q['G0_discrete_J_m2']:.16g} J/m²。它们是理想临界跨面 point-state 的标称能量，未测得实际等效扩展阻力 GC_EFFECTIVE；历史 PVCE G_DRIVE 是给定载荷下试样驱动力，也不是材料 Gc。

| 参考 | 量值（J/m²） | 材料、定义与局限 |
|---|---:|---|
| [Liu 2021](https://doi.org/10.1016/j.eml.2021.101444)，PDF第1、4页 | 5.90±1.18 | 人体全血凝块，实验总断裂能；不能与理想弹性非局部分配逐项等同 |
| [Liu 2024](https://doi.org/10.1016/j.jmps.2024.105858)，PDF第6页 | 18.27±1.45 | 牛血栓、最低加载率0.01 s⁻¹下的纯剪切总撕裂能；明显依赖速率 |
| 同一Liu 2024，PDF第10–11页Fig.8 | 约6–8 | Yeoh/Mullins与双线性黏聚模型约束反演；结合强度及最低速率总能量，不是另一组直接Gc实测值 |

本轮已读取这两篇用户提供的合法本地PDF相关正文页，路径与SHA见LITERATURE_SOURCE_IDENTITY.json。6–8 J/m²依赖模型参数关联和强度条件，不能把它作为与18.27 J/m²相同定义的误差带。细节可在第11页关于S_bulk超过约5.2 kPa时的推断段落核对。

本轮0.0697 J/m²明显低于上述参考尺度，约为前两项均值的{100*q['G0_discrete_J_m2']/5.90:.4f}%和{100*q['G0_discrete_J_m2']/18.27:.4f}%。这里的比值只表达数值尺度，**不是材料Gc实验误差或已经认证的韧性偏差**。结论仅为 `NOMINAL_ENERGY_SCALE_BELOW_EXPERIMENTAL_REFERENCES`。材料种属、制样、加载率、有限变形、Mullins/黏性耗散及内禀功定义均不同。

当前LPS只保留由Liu材料卡映射的小应变K/G、固定伸长门槛，没有据本轮结果调参，也没有足够证据判定必须新增渐进软化。state-based膨胀量耦合邻居；逐条相互作用不能当成独立保守弹簧。本轮没有追加固定坐标断键诊断，历史带符号拓扑项继续保留其原定义。

下一轮唯一建议：冻结一套受控Mode-I分离试验，测量是否形成可追踪的自然裂纹面积，记录载荷—开度及外功、初末弹性储能、动能、接触储能、数值阻尼和带符号拓扑项；据明确的面积和能量分解识别实际等效阻力及其加载依赖。要验证的是实际传播/耗能与本轮理想标称值的关系，不能预设二者相等。本轮未启动该试验。
'''
(O/'MODEL_VS_EXPERIMENT_SCOPE.md').write_text(scope)
# A stable snapshot whose own manifest/report files do not recursively count themselves.
excluded={'OUTPUT_SIZE_AT_REVIEW.json','ARTIFACT_MANIFEST.json','WSL_DELIVERY_VERIFICATION.json','SERVER_FINAL_PACKAGE_CHECK.json','FINAL_READBACK.json','REVIEW_LPS_3D_FRACTURE_SURFACE_ENERGY_IDENTIFICATION_V1.md'}
payload=[p for base in [O,F,R/'analysis/fracture_surface_v1'] for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in excluded]
size={'captured_utc':datetime.now(timezone.utc).isoformat(),'scope':'Current new-task source/data/figures/evidence snapshot, excluding this size record, final report and self-referential transfer/manifest records. Exact final deliverable inventory is ARTIFACT_MANIFEST.json.','file_count':len(payload),'bytes':sum(p.stat().st_size for p in payload),'figure_bytes':sum(p.stat().st_size for p in F.glob('Fig*')),'full_precision_npz_bytes':(O/'LAMBDA_AND_ENERGY_DENSITY.npz').stat().st_size,'excluded_names':sorted(excluded)}
dump(O/'OUTPUT_SIZE_AT_REVIEW.json',size)
layer=[a for a in load('LAYER_PROFILES') if a['plane_index']==0 and 0<a['z_m']<DELTA]
layertext='；'.join(f"距面{a['distance_m']*1e6:.0f} µm：λ={a['lambda_mean']:.12f}" for a in layer)
qanda=[
('三维状态型公式是否从2018原文核查？','是。WSL完整19页PDF已阅读，重点核对Eq.(2)–(4)、(9)–(12)、(17)–(22)、(26)–(27)，并直接查看第3、5、6页公式和Fig.2。使用三维OSB-PD，不用二维或键型公式；详见PAPER_EQUATION_AUDIT.md。'),
('理论G0准确值？',f"{t['G0_theory_J_m2']:.17g} J/m²（float64第一式）；第二式{t['formula2_J_m2']:.17g}，差{t['difference_J_m2']:.3e}。解析β/δ=1/5、β′/δ=26/405用Fraction独立积分核对。"),
('体积和剪切贡献？',f"连续K项{t['K_contribution_J_m2']:.14g}、G项{t['G_contribution_J_m2']:.14g} J/m²，K占{100*t['K_fraction']:.6f}%。离散K项{q['K_contribution_J_m2']:.14g}、G项{q['G_contribution_J_m2']:.14g} J/m²，K占{100*q['K_fraction']:.6f}%。K−5G/3是混合系数，未误标纯体积项。"),
('影响函数与论文是否相容？','在本轮完整参考邻域上相容。源码ω=1/r，论文ω=δ/r；同时缩放m及theta/偏能的全部积分，归一化抵消常数。复核同一保存弹性态，点状态能量与旧独立C++能量相对差3.874×10⁻¹³，旧原生力与独立力差1.864×10⁻¹⁴。本轮没有新的原生力调用。'),
('partial-volume taper如何进入？','实际中心键先满足r0≤δ，再在|r0−δ|≤Δx/2时用h=(δ+Δx/2−r0)/Δx，否则h=1。h同时进入m_full、m_cross、theta和全部偏伸长储能。√8Δx壳12邻居、h≈0.681572875；3Δx壳30邻居、h=0.51。未补造δ外键，也未做切换taper的额外能量扫描。'),
('是否真正完整未切参考表？','是。实际原reference_topology.npz的SHA为ef145b6d6db6c17584dce6cbea1fa1b9d95847d71a69265a73df6334d1695c03，与20261011T154854Z原manifest一致。909392个唯一键及独立整数邻域枚举完全一致。NPZ内A0/A1/A2的alive掩码不用于本轮；它们分别删16416/16976/17536键，不能代替完整参考。'),
('m是否对应真实atom ID？','是。原生body.data逐行读回16000个实际atom ID及全精度坐标、909392个原bond ID/端点，均与原NPZ逐位一致。ID=((ix×8+iy)×20+iz)+1，原m按此顺序保存。完整邻域重算m与原数组相符；最终双修复将m/ready作为每原子属性随排序、迁移、restart保存，本轮未修改此机制。'),
('周期Y相互作用是否正确？','是。B=0.4 mm>2δ，采用唯一最小镜像；131196条跨周期缝键已核对，Y平移一个点距后全键集合与m不变。无额外镜像键重复加入。材料体积外界面不同于LAMMPS非周期方向的带余量计算盒。'),
('真实测量面积？','X中央30列，半开窗口[-0.75,0.75) mm；完整Y周期8层、B=0.4 mm。A=(30×50 µm)×0.4 mm=0.6 mm²=6×10⁻⁷ m²。不是历史PVCE的BΔa=2×10⁻⁸ m²。未事后调整ROI。'),
('两侧是否重复计算？','没有。4800个ROI点中1440点非零、上下各720。每条唯一键为两端点的各自邻域提供贡献，最终各点W_i V_i计一次；下侧与上侧各0.03484811867116883 J/m²，求和即主结果，没有再乘2。'),
('lambda范围与分布？',f"全试样λ范围{c['lambda_all_range'][0]:.6g}–{c['lambda_all_range'][1]:.12f}（含X边缘点）；非零ROI点范围{c['lambda_ROI_active_range'][0]:.12f}–{c['lambda_ROI_active_range'][1]:.12f}。{layertext}，上下镜像各一层，每层240点。其他ROI层λ=0。未进行裁剪。"),
('跨面键数量？','全试样跨z=0唯一键55408条，其中7324条同时跨Y缝；与ROI至少一个端点相连的跨面键17392条（其中2296条跨Y缝），两端都在ROI内16208条；ROI端点的有向跨面贡献总数33600。这些口径分别记录，不用任何一个键数替代面积。'),
('当前G0_discrete？',f"{q['G0_discrete_J_m2']:.17g} J/m²，ROI理想能量{q['U_ideal_ROI_J']:.17g} J。来自实际保存的完整参考邻域和原生求积规则；只验收理想标称能量。"),
('两种独立表达式是否一致？',f"是。第一式从m_cross/m_full计算lambda；第二式从全部原始相互作用e构造theta与偏伸长，保留非跨面e=0但偏伸长不为零的项。主结果总能量相对差{e['relative_global_error']:.6e}，逐点差按max(W)归一化为{e['pointwise_error_normalized_max_W']:.6e}。非跨面偏能贡献0.00192221583766 J/m²，不能遗漏。"),
('与0.07170070833相差多少？',f"按未四舍五入理论值，离散−理论={q['difference_from_theory_J_m2']:.14g} J/m²（约−2.004471 mJ/m²）；比值{q['G0_over_theory']:.14g}，相对−2.79560835%。没有人为5%验收标准，也不据接近程度认证物理模型。"),
('哪些差异来源已确认？','实际离散β=29.1389395961 µm、β′=9.39752186164 µm；连续分别30.1和9.66172839506 µm。由这些求积矩可精确分解总差：K项−0.00178339410058、G项−0.000221076890414 J/m²。原生122邻居的有限方向/壳采样、最外层taper和单侧3层z求积是已存在的离散规则；本轮未隔离各因素的独立因果份额。所有有贡献ROI点均为完整122邻域，主面最小外界面距离375 µm、对照325 µm，均大于150.5 µm；计数、面积及归一化错误已排除。'),
('一次位置/相位检查结果？',f"z=50 µm为另一真实间隙，G0={d['cases'][1]['G0_discrete_J_m2']:.17g} J/m²，与主面差{sens['G_difference_J_m2']:.6e}；按Z平移一层对齐后λ最大差{sens['lambda_mapped_max_abs_difference']:.3e}。两者使用不同实际键ID集合，但局部拓扑平移等价、同网格相位；因此只支持内部平移稳定，不能宣称独立相位或网格收敛。"),
('发现新的明确源码错误了吗？','没有。最终双修复源码在本轮完整弹性邻域上与独立势能逐项对应。已知native PE偏能计量不足及断键时lagged theta限制仍在；未将理论/离散差异当成新bug。Eq2/Eq9的本轮恒等依赖完整m与完整求积家庭，不外推到删除键而m仍固定的情况。'),
('改动生产C++或参数了吗？','没有。K、G、s00、delta、dx、alpha不变，生产库、最终双修复库、原能量.so、源码、材料文件和env.sh均经SHA保护。新增的是独立Python只读积分、审查和绘图代码；本地PDF解析在独立临时venv，未改服务器环境。'),
('能称已得到材料固有Gc吗？','不能。理想跨面伸长状态未必有全局可实现位移，本轮没有真实自然传播、增量裂面识别与完整分离功账本。G0_THEORY、G0_DISCRETE、历史G_DRIVE与GC_EFFECTIVE必须分开；最后一个仍未测量。'),
('与Liu实验可比和不可比之处？','只能比较J/m²量级。本轮0.0697低于人体全血总断裂能5.90±1.18、牛血栓低速率总撕裂能18.27±1.45，以及特定模型反演内禀功6–8。物种、制样、加载率、有限变形及耗散定义不同，不计算“已认证Gc误差”。详见MODEL_VS_EXPERIMENT_SCOPE.md及本地原文身份。'),
('是否必须升级渐进软化？','本轮证据不足。数值尺度偏低只触发后续实际分离阻力识别；没有改变s00或加入新的耗散/软化机制。'),
('下一轮Mode-I应验证什么？','唯一建议是受控Mode-I分离：明确并追踪自然新生裂面，测量载荷—开度、外功、储能、动能、接触、数值阻尼及带符号拓扑项，检验实际等效阻力与理想标称值的关系。不能先假定二者一致。本轮到此停止，未启动下一轮试验。'),
('真实耗时、资源和输出规模？',f"服务器cgroup为15.36 CPU核配额、42.07 GiB内存上限；本任务单进程，OMP/BLAS各1线程，不用4090。输入读取{perf['input_read_s']:.6f} s；两项真实积分{perf['integration_s']:.6f} s；含读写的积分进程{perf['process_wall_s']:.6f} s、CPU {perf['process_CPU_s']:.6f} s、峰值{perf['peak_RSS_MiB']:.3f} MiB。准备/身份核对{prep['wall_s']:.6f} s、峰值{prep['peak_RSS_MiB']:.3f} MiB，另有原文件SHA检查及绘图时间，不能把0.068 s说成全任务耗时。审查时输出快照{size['file_count']}文件、{size['bytes']} B；最终清单见ARTIFACT_MANIFEST.json。未观察不可复现问题，但未额外重复正式积分。"),
('生成哪些科研图？','四组均有PNG/SVG及OPEN_RESULTS.html：①真实点阵、z=0、ROI和真实周期参考邻域；②lambda、跨面m、K/G/总密度与层剖面；③从零起轴的理论/主值/对照及贡献差；④源码对应、真实taper壳层、周期/双侧计数与已知限制。绘图只使用保存结果，图1全16000点，键显示有明示的确定性子集和原ID清单；无连续医学壳或AI生成数据。'),
('历史数据完整吗？',f"是。对{hist['files_checked']}个保护文件、{hist['bytes_checked']} B重新核对SHA，缺失0、变化0；PVCE V1/V2、宽薄历史A0/A1/A2和最新CFD接口均保留。没有git reset/clean、stash/push、rsync删除或生产部署。"),
('WSL clot_rupture.md实际更新了吗？','是。通过科学检查后实际更新/home/lzy/projects/clot_rupture.md，在⑤段内增加本轮公式、数值、定义限制和Mode-I建议，文末追加本轮链接；原①–⑨标题及全部旧正文逐字保留。双次并发检查、原子写入和实际读回均通过。'),
('备份、diff与SHA在哪里？',f"本run的wsl_document/含clot_rupture.before.md、clot_rupture.after.md、clot_rupture.diff及DOCUMENT_UPDATE_AUDIT.json；入口WSL_DOCUMENT_UPDATE.json。before SHA `{doc['before_sha256']}`；after/readback SHA `{doc['after_sha256']}`。服务器和WSL各保留这套交付证据。")]
report=f'''# LPS_3D_FRACTURE_SURFACE_ENERGY_IDENTIFICATION_V1

run_id：`{C['run']}`。本轮已完成服务器只读理论/拓扑/储能核对、主裂面及一次对照积分、四组图和实际WSL文档更新。

**当前 s00=0.15、50 µm 真实点阵的标称裂面能量为 {q['G0_discrete_J_m2']:.14f} J/m²。这个数值不足以代表真实自然裂纹扩展阻力；GC_EFFECTIVE 未测量，材料Gc未完成实验认证。**

| 结果口径 | K项（J/m²） | G项（J/m²） | 合计（J/m²） |
|---|---:|---:|---:|
| G0_THEORY | {t['K_contribution_J_m2']:.14f} | {t['G_contribution_J_m2']:.14f} | {t['G0_theory_J_m2']:.14f} |
| G0_DISCRETE，z=0 | {q['K_contribution_J_m2']:.14f} | {q['G_contribution_J_m2']:.14f} | {q['G0_discrete_J_m2']:.14f} |
| 同相位内部位置对照，z=50 µm | {d['cases'][1]['K_contribution_J_m2']:.14f} | {d['cases'][1]['G_contribution_J_m2']:.14f} | {d['cases'][1]['G0_discrete_J_m2']:.14f} |

离散/理论={q['G0_over_theory']:.12f}；差约−2.7956%。两种独立离散表达式相对差3.16×10⁻¹⁶。离散定义的数值一致性通过，不等于有限点阵已收敛或实验韧性已匹配。

原始server工作树分支`migration-20261009`、HEAD `cc5ee77149f55e58d715bb5a4660c6467231a826`，与GitHub选择性快照HEAD `32e95bf69b76d1856e5f872f5e5d099522a14959`分开记录。预存改动见PREFLIGHT.json，未清理。最终双修复库SHA `{EXPECTED}`；生产库SHA `{hist['production_library_sha256']}`。加载最终库仅为ctypes身份读取，没有创建LAMMPS实例。

## 28项科研审查

'''+ '\n\n'.join(f'### {i}. {question}\n\n{answer}' for i,(question,answer) in enumerate(qanda,1))+'\n\n## 分项验收状态\n\n| 项目 | 状态 | 认证边界 |\n|---|---|---|\n'+ '\n'.join(f'| {k} | {v["status"]} | {v["scope"]} |' for k,v in statuses.items())+f'''

## 准备过程与复现

首次准备脚本在迁移后的SOURCE.json旧路径不存在处停止。已保留PREFLIGHT_INITIAL_PATH_PROBE.json、PREPARATION_NOTE.json及日志，定位到真实isolated_checkout历史来源后完成准备；这发生在设计冻结和正式积分之前，不是求解器失败，也没有物理重跑。实际最后库和源码SHA来自服务器当前双修复目录，SOURCE.json只证明上游/迁移补丁来源，不冒充最终双修复manifest。

设计于{freeze['frozen_utc']}冻结，SHA `{sha(O/'DESIGN_FREEZE.json')}`；正式积分于{perf['completed_utc']}完成。严格为一主一对照，无额外材料、网格、相位、taper模型扫描。旧保存非仿射态仅使用已有独立.so评估一次以复核能量对应，不是新增平衡或动力运行。图稿只因可读性调整，不变更计算结果。

复现命令、不得覆盖旧run的约定见REPRODUCE.txt；真实输入/源码见SOURCE_AND_LIBRARY_IDENTITY.json、HISTORICAL_INPUT_MANIFEST_EXTRACT.json、NATIVE_INPUT_IDENTITY_AUDIT.json及source_inputs/。节点输出保留float64，包含原atom ID、X、V、m、两套lambda/theta/W、ROI和原键序的跨面掩码；交付不复制11716992 B原生库或完整论文，保留精确路径与SHA。

四组图入口：[OPEN_RESULTS.html](OPEN_RESULTS.html)。最终文件SHA清单：[ARTIFACT_MANIFEST.json](ARTIFACT_MANIFEST.json)。回收与本地逐文件核对：[WSL_DELIVERY_VERIFICATION.json](WSL_DELIVERY_VERIFICATION.json)。后续仅建议受控Mode-I分离，本轮停止。
'''
(O/'REVIEW_LPS_3D_FRACTURE_SURFACE_ENERGY_IDENTIFICATION_V1.md').write_text(report)
print(json.dumps({'report_questions':len(qanda),'statuses':len(statuses),'output_snapshot':size},ensure_ascii=False))
