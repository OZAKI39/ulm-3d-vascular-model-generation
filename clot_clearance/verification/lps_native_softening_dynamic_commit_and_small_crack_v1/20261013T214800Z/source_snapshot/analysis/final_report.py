from common import *
import shutil,sys,time,subprocess,re
m=json.loads((O/'DYNAMIC_341_NODE_AUDIT.json').read_text());e=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());h=json.loads((O/'DT_HALF_WINDOW_COMPARISON.json').read_text());r=json.loads((O/'RUN_RESOURCES.json').read_text());i=json.loads((O/'FINAL_LIBRARY_IDENTITY.json').read_text());doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text());protect=json.loads((O/'HISTORICAL_PROTECTION_AFTER.json').read_text());assert doc['status']=='PASS' and protect['status']=='PASS'
assert sha(O/'clot_rupture.after.md')==doc['after_SHA256'];assert sha(LIB)==i['library_sha256']
items={
'DEPRECATED_CODE_AUDITED':'PASS_28_CANDIDATES',
'PROVEN_DEAD_CODE_REMOVED':'NO_PROVEN_SAFE_DELETION_0_FILES',
'CLEANUP_REVERSIBLE':'PASS_NO_DELETION_RESTORE_UNNEEDED',
'ACTIVE_CODE_ENTRYPOINT_CREATED':'PASS',
'HISTORICAL_EVIDENCE_PRESERVED':'PASS_31001_FILES',
'PRODUCTION_LIBRARY_UNCHANGED':'PASS',
'FROZEN_SOFTENING_MODEL_PRESERVED':'PASS',
'DYNAMIC_INTEGRATOR_CONNECTED':'PASS_ISOLATED_NATIVE_VERLET',
'ACTUAL_FREE_NODE_MOTION':'PASS_245_FREE_NODES',
'PHYSICAL_TIMESTEP_IDENTITY':'PASS_OBSERVED_244_MAIN_STEPS',
'ONE_COMMIT_PER_ACCEPTED_STEP':'PASS_OBSERVED_WINDOW',
'RUN0_AND_TRIAL_READ_ONLY':'PASS_AFTER_PYTHON_ASSERTION_REPAIR',
'DAMAGE_HISTORY_IRREVERSIBLE':'PASS_OBSERVED_WINDOW_UNLOAD_NOT_REACHED',
'FIRST_SOFTENING_EVENT':'OBSERVED_STEP170_85US_BOND3393',
'FIRST_COMPLETE_FAILURE_EVENT':'OBSERVED_STEP185_92_5US_7_BONDS',
'ORIGINAL_BOND_IDENTITY':'PASS_SINGLE_RANK_DYNAMIC_AND_RESTART',
'DYNAMIC_FORCE_CONSISTENCY':'PASS_ALL_245_SAVED_STATES',
'DYNAMIC_ENERGY_LEDGER':'PASS_WITHIN_STOPPED_WINDOW_RESIDUAL_0_144631_PERCENT',
'MOMENTUM_CONSISTENCY':'PASS_OBSERVED_WINDOW',
'TIMESTEP_SAFETY':'STOPPED_MAX_DD_0_102554_GREATER_THAN_0_10',
'POSTBREAK_CONTACT_CONSISTENCY':'ZERO_FAILED_BOND_TENSION_VERIFIED_ACTIVE_CLOSURE_UNOBSERVED',
'SMALL_NONAFFINE_CRACK_VALIDITY':'NOT_RUN_P3_ADMISSION_NOT_PASSED',
'CRACK_AREA_IDENTIFIABILITY':'UNAVAILABLE',
'MATERIAL_GC_NOT_CERTIFIED':'NOT_CERTIFIED_GC_UNAVAILABLE',
'MACROSCOPIC_STRENGTH_NOT_CERTIFIED':'NOT_CERTIFIED',
'FULL_V3_NOT_RUN':'NOT_RUN',
'PRODUCTION_SCALABILITY_NOT_CERTIFIED':'NOT_CERTIFIED',
'HISTORICAL_RESULTS_PROTECTED':'PASS',
'WSL_DOCUMENT_UPDATE':'PASS_ACTUAL_FILE_UPDATED'
}
final=dict(task=TASK,run_id=RUN,overall='IMPLEMENTATION_AND_EVIDENCE_DELIVERED_P3_STOPPED_AS_REGISTERED_P4_NOT_RUN',phases=dict(P0='PASS_ZERO_SAFE_DELETIONS',P1='FROZEN_BEFORE_CPP',P2='ISOLATED_IMPLEMENTATION_VERIFIED',P3='STOPPED_TIMESTEP_DAMAGE_INCREMENT_GATE',P4='NOT_RUN',P5='EVIDENCE_FIGURES_REPORT_DOCUMENT_DELIVERED'),items=items,main=m,short_window=h,material_Gc='UNAVAILABLE',macroscopic_strength='NOT_CERTIFIED',new_crack_area='UNAVAILABLE',production_deployed=False,further_computation_started=False,next_only_priority='同一341点冻结工况的损伤时间步准入与全时程能量验证；材料、夹具与加载不变，先预注册完整较小步长协议，通过后再审查P4。',answer_1='在已保存0–122µs窗口内，真实原生VV产生不可逆自然软化与完全失效，力、历史、动量一致，能量误差可解释且已量化；但主时程触发冻结Δd门槛而停止，闭合接触和完整时程未认证，不能判为整体动力学验收通过。',answer_2='没有足够证据进入更大规模非均匀裂纹或宏观强度/断裂能标定。P4未运行，Gc不可用。')
dump(O/'FINAL_STATUS.json',final)
answers=[
('P0之前实际存在多少候选废弃代码？','受限入口审计共28个候选源码文件；“候选”不是已经证实废弃。包括当前两套analysis模块、两份冻结原生源码及CMake编译探针。详见CODE_CLEANUP_INVENTORY.json。'),
('哪些代码被确认已经无效？','没有文件同时满足“确定不被调用且无复现/历史价值”的安全删除条件。正式源和复现代码继续保留，编译探针用途不确定。'),
('实际删除了哪些文件？','删除0个；DELETED_CODE_MANIFEST.json路径列表为空。未通过改名、移动或删除历史数据制造清理成果。'),
('删除依据是什么？','白名单为空；实际根据限定Python imports、C++ includes、冻结manifest、构建依赖、Git状态及科研来源判断。无充分依据的文件不删。'),
('替代代码在哪里？','没有需要替代的已删文件。本轮新活跃入口analysis/dynamic_softening_v1/，独立原生源码在本轮source/src/BPM；冻结离线参考仍为analysis/softening_feasibility_v1/model.py。'),
('是否仍有不确定的历史代码？','有，清单中的UNCERTAIN记录保留，尤其CMake编译探针。未遍历历史实验全文作过度清理。'),
('删除是否可恢复？','没有删除，无需恢复。原始候选SHA、Git before/after及所有受保护历史SHA均记录；旧科学源快照完整保留。'),
('是否建立了精简的活跃入口？','已建立根目录ACTIVE_DEVELOPMENT_ENTRYPOINTS.md。P0版27行单独保留；最终入口补充停止状态与唯一后续优先项，总行数远低于150。'),
('是否减少了下一轮需要读取的代码范围？','入口明确当前四个C++文件、冻结model.py、动力学驱动和关键审计。后续可按需读取，不必展开旧PROJECT_ACTIVE_CONTEXT或全部verification。未虚构量化的上下文节省比例。'),
('是否保护了历史科研证据？','实际逐文件复核31001个历史文件和34项受保护资产，冲突0；已有tracked Git diff逐字相同。P0先于任何新C++实现。'),
('实际采用什么Velocity-Verlet调用路径？','锁定LAMMPS Verlet::run → modify.initial_integrate（free组fix nve两次半步中的第一半步与漂移、grip组fix move）→新fix POST_INTEGRATE→邻居重建→pair/bond→reverse_comm→POST_FORCE力密度乘V→final_integrate第二半步。245个自由节点三分量真实积分。'),
('损伤在哪个物理时间层提交？','在x_(n+1)已漂移、第二次速度半步尚未执行时。先用d_n计算同坐标U_old，再提交κ/d，计算U_new与带符号D增量；新力使用d_(n+1)。这是明确的损伤算子分裂。'),
('每个接受时间步是否只提交一次？','主时程真实step1–244与commit计数、token逐步相同；从170步restart续算计数不重置。短窗口另有68个真实步，原生整数step211–278，对应物理105–122µs。无force调用即接受的混淆。'),
('run0和重复force evaluation是否保持只读？','是，实际重复run0、改变试探坐标求力和恢复均不改变κ/d/D/事件/提交token；有动态fix时显式外部commit被拒绝。第170步的Python断言曾错误比较不同坐标的U，失败记录保留；修正后同接受态续算，未更改C++或物理协议。'),
('是否发生真实自由节点运动？',f'发生。245个free点位置和速度由原生力更新，停止态相对参考最大自由位移{m["free_motion_max_m"]*1e6:.9f}µm。独立核对每一步两次半更新方程，位置最大误差{m["position_equation_max_m"]:.5e}m，速度误差{m["velocity_equation_max_m_s"]:.5e}m/s。'),
('实际dt是多少？','主dt=0.5µs；初始完整切线给出dt*omega_max=0.07190393、线性临界dt=13.90744591µs。实际损伤门槛触发后，仅105–122µs一个窗口使用dt=0.25µs。没有dt/4或全程dt/2重跑。'),
('首次软化时间与原始bond ID是什么？','第170步、85µs，原始键3393。全部历史重建与原生κ/d逐值相同。'),
('首次完全失效时间与bond ID是什么？','第185步、92.5µs，同步首失效7条：2180、2194、3368、3393、3904、9197、10264。所有新增失效时的实际拉伸均达到sc=0.18。停止时520个唯一原始键各失效一次。'),
('κ/d是否有任何回退？','已保存主时程及短窗口无回退，固定坐标D增量非负。卸载段没有达到，不能把此范围描述成动态全循环不可逆性验证。'),
('完全失效后是否仍有不合理非接触传力？','已失效边的两个端点标量t均为0，原始边非接触拉力为0；固定停止态梯度最大相对误差1.3468e-9。整个试样仍连通，不宣称本轮形成了完全分离的大块。'),
('接触何时启用？','d首次到1的POST_INTEGRATE内删除相应owned键和1–2屏蔽，随后同物理步邻居重建、pair计算可启用1.0dx接触。当前保存的自然张开时程没有实际闭合接触，活动接触步数0，因此动态闭合覆盖状态为UNOBSERVED。'),
('动力学速度与步间位移是否合理？','最大速度0.4443383054m/s，步位移/dx=0.0044333927，相对端点步位移/dx=0.0078199333，均低于冻结数值门槛。但最大Δd=0.1025542483超过0.10，故停止；这不是已认证生物加载速率或全程步长。'),
('是否存在数值异常冲量？',f'已观察步的实际总动量与含夹持惯性反力的外部冲量余额最大{m["momentum_residual_max_Ns"]:.6e}N·s，未发现额外异常冲量。损伤后的正常力重分配和有限时间分裂误差另计。'),
('是否执行P4？','NOT_RUN。主P3损伤增量门槛未过，严格停止，没有创建或运行第二种试样。'),
('为什么选择该试样？','本轮唯一主试样复用历史341点/12047键，能直接对比冻结材料和原始身份。由于P3停止，没有对P4几何作虚构选择。'),
('实际预裂纹与加载边界是什么？','P4缺口与加载均未创建。P3没有新增预切键；原12047条键初始完整，左右x端49/47点规定对称位移，其余245点自由，夹持点y/z固定。幅值每侧60µm，升载300µs、保持50µs、卸载300µs、末保持50µs，冻结700µs；实际只到122µs。'),
('损伤是否从预期区域发生？','P4 NOT_RUN，没有预期裂尖可验。P3实际损伤分布见Fig3和FP64数据；该341点阶段只做动力学机制诊断，不追认它为缺口定位实验。'),
('是否出现真实裂纹扩展？','观察到真实非均匀损伤和原始键自然失效；没有证据将其认证为从预定缺口前缘推进的裂纹。P4没有运行。'),
('是否形成可辨识完整裂面？','没有建立可辨识新裂面，停止态为一个341点连通体。部分d、520条失效键或节点颜色不能代替裂面识别。'),
('是否发生非目标损伤？','P4未执行，目标/非目标带指标UNAVAILABLE。P3只保存实际分散损伤，不能直接用P4的目标裂尖语义分类。'),
('是否仍有跨面存活键？','有。作为明确标注的只读拓扑诊断，x=20µm平面仍有2041条存活跨面键；这不是预注册裂纹面积测量面。完整存活键11527条。'),
('是否能够报告面积？若不能，原因是什么？','新增裂面面积UNAVAILABLE，不能以0替代。P4没有准入、无已识别连续完全分离面；不用断键数乘dx²。'),
('实际边界外功是多少？',f'到122µs为{e["terms"]["external_work_J"]:.16e}J，即5.013269603nJ。使用真实边界位移与反力的梯形功，加规定边界动能增量；是数值路径功估计，不是终点F·u。'),
('完整变分储能如何变化？',f'从0到{e["terms"]["stored_J"]:.16e}J。来自完整端点邻域自由能；逐态独立重建最大误差1.0960e-23J，未使用旧瞬断不完整PE替代。'),
('实际动能如何变化？',f'从0到{e["terms"]["kinetic_J"]:.16e}J，全部材料点含规定边界；独立sum(0.5*rho*V*v²)核对native KE，最大差6.2039e-24J。free点动能也独立存表。'),
('不可逆损伤耗散是多少？',f'累计{e["terms"]["damage_J"]:.16e}J，即0.461696165nJ；每步由同一新坐标的U_old−U_new实际计算，独立增量最大差6.2039e-24J。未按残差反推D。'),
('数值阻尼和接触贡献是多少？','没有阻尼，阻尼功恰为0；实际活动接触为0，Uc和接触拓扑注能为0。提取Ftotal−Fbond的舍入残量所积接触功约8.30e-28J，不把它称为实际闭合接触。带符号材料拓扑变化−0.461696165nJ与D只计一次。'),
('动态能量账本是否闭合？',f'在停止窗口按W=ΔU+ΔK+D+ΔUc−Jc+R核对，R={e["terms"]["residual_J"]:.16e}J；最大相对余额{100*e["relative_residual_max"]:.9f}%，低于冻结2%门槛。它是有有限误差的窗口通过，不能称整个700µs协议已通过。'),
('新模型是否出现未经解释的能量注入？','未发现未解释来源。余额主要来自损伤后力切换功+7.226909pJ，另有Verlet求积+0.046242pJ和固定d空间求积−0.011886pJ；其负和重建R，差2.82e-23J。这仍是有限时间误差，不宣称精确守恒。'),
('当前是否可以报告材料Gc？','不可以，UNAVAILABLE。旧G0_THEORY=0.07170070833333J/m²和G0_DISCRETE=0.06969623734234J/m²保留为旧瞬断模型的标称量；不等于新软化材料Gc。历史173.156623/230.776817nJ也不直接换算成本轮Gc。'),
('新模型是否已证明宏观强度不变？','没有。冻结未损伤刚度与s0不等于宏观峰值强度已认证，341点短动力学窗口不能替代材料实验和尺度验证。'),
('本轮是否有必要修改sc或进一步升级本构？','没有支持这种修改的证据。main sc=0.18保持；sc=0.24仅一个固定状态回归，未另跑动态。历史新旧二值力L2差8.814180703270e-6N原样保留；本次停止坐标同二值掩码差1.239512740840495e-5N属于两种自由能的本构差异。当前优先问题是时间分辨率。'),
('实际修改了哪些C++源码？','本轮独立副本中的bond_bpm_peri_variational.h/.cpp，加动态提交入口、同坐标缓存及完整全失效快速返回；新增fix_lps_variational_dynamic.h/.cpp负责hook、物理力换算和运行约束。pair、Verlet、NVE及生产源均不改。src/.gitignore仅补包文件登记。'),
('新隔离共享库SHA是多少？',i['library_sha256']+'。本轮全部正式步实际加载此SHA；FINAL_LIBRARY_IDENTITY.json还记录四个源码、最终补丁及旧生产库SHA。'),
('编译和运行总成本是多少？',f'P0 {r["P0_wall_s"]:.3f}s；首次编译与一次生命周期析构防护修正的增量编译合计{r["build_wall_s"]:.3f}s，12-job上限。主244步{r["main_wall_s"]:.3f}s，短窗口68步{r["half_window_wall_s"]:.3f}s，正式动力合计{r["formal_dynamic_wall_s"]:.3f}s、CPU {r["formal_dynamic_CPU_s"]:.3f}s，峰值{r["peak_dynamic_RSS_bytes"]/1024**2:.3f}MiB。主平均{r["main_wall_per_step_s"]*1e3:.3f}ms/步包含观测与保存，不是裸内核性能。零步表达式失败、静态回归和图形费用在资源表分开；没有虚报未计时编排耗时。'),
('是否发生任何历史文件意外变化？','31001个历史文件及34项关键资产SHA无冲突，已有tracked diff未变。只有本轮新路径、活跃入口和明确授权的真实WSL文档发生预期修改。'),
('是否保留生产库？','保留。生产库SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a；旧双修复库00b246...和上轮b80bfc...亦完全保持，env.sh未换库。'),
('是否生成必要科研图？','实际生成5组PNG+SVG及OPEN_RESULTS.html：代码结构、物理时间/提交/损伤、P3实际停止点云、完整能量、时间步与分项准入。P4未运行，因此第三组明确为P3停止证据，未伪造缺口推进或医学连续外壳。没有动画。'),
('WSL文档是否真正更新？',f'实际修改/home/lzy/projects/clot_rupture.md，读回一致；SHA由{doc["before_SHA256"]}变为{doc["after_SHA256"]}。保留原①–⑨结构、中文文风、材料参数和全部旧结果。'),
('备份、diff、删除清单和SHA在哪里？','本轮verification目录内：clot_rupture.before.md / after.md / diff、WSL_DOCUMENT_UPDATE.json、DELETED_CODE_MANIFEST.json、NATIVE_SOURCE_PATCH.diff、FINAL_LIBRARY_IDENTITY.json、HISTORICAL_PROTECTION_AFTER.json、PAYLOAD_MANIFEST.json。source_snapshot保留驱动/四个C++与读取的核心依赖，binary保存实际库。'),
('下一轮唯一最合理的开发任务是什么？','同一341点冻结工况的损伤时间步准入与全时程能量验证。先预注册完整较小dt协议，保持材料和夹具；不能把一个17µs窗口当作全程通过。通过后再审查P4，本轮不自动续算或开展更大裂纹/宏观标定。')]
assert len(answers)==51
intro=f'''# LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1

运行目录：`{O}`。本轮严格先P0、再P1/P2、一个341点真实动力学工况，最后一次由实际风险触发的短dt/2窗口；P4未准入。

**核心结论：在0–122µs实际窗口内，原生Velocity-Verlet已产生不可逆自然软化和完全失效，自由节点真实运动，力、历史及动量一致，能量余额被独立重建。但主时程第244步的最大Δd=0.102554超过冻结0.10门槛，已停止；不能判为完整动力学验收通过。没有足够证据进入更大裂纹、宏观强度或材料Gc标定。**

## 运行前冻结与实际边界

G=940Pa，K=100000/3Pa，ρ=1000kg/m³，dx=50µm，horizon=150.5µm，s0=.15、主sc=.18，α=0。真实物理质量ρV=1.25e-10kg；原生BPM力密度由新fix在第二半步前转换为N。245自由点、左右49/47夹持点。完整加载预注册700µs，实际122µs后严格停止，没有临时增幅。图和表只用真实已保存数据。

短窗口105–122µs从同一接受restart出发，0.5/0.25µs对应34/68步。最大Δd由{h['coarse_max_dd']:.10f}降为{h['fine_max_dd']:.10f}；能量余额−6.051331549pJ降为−3.002602740pJ，比值{h['residual_ratio_fine_coarse']:.9f}；同一终点失效ID集合完全相同。最大坐标分量差1.226177nm、速度分量差0.000123467m/s。仅认证局部比较，不外推全局时间收敛。

## 完整能量账本

$W=\\Delta U+\\Delta K+D+\\Delta U_c-J_c+R$。初态U=K=Uc=D=0，数值阻尼不存在。

| 项目 | 停止值 / J | 定义与计量 |
|---|---:|---|
| 实际边界功W | {e['terms']['external_work_J']:.16e} | 边界实际位移×梯形材料反力，加精确边界ΔK；数值求积估计 |
| 完整变分U | {e['terms']['stored_J']:.16e} | 原生完整端点邻域自由能，独立重建 |
| 实际K | {e['terms']['kinetic_J']:.16e} | 全节点物理质量、接受整步速度 |
| 不可逆D | {e['terms']['damage_J']:.16e} | 同新坐标的U_old−U_new累计 |
| 数值阻尼 | 0 | 未设置阻尼 |
| 接触Uc / 接触拓扑注能Jc | 0 / 0 | 本窗口没有活动接触 |
| 带符号材料拓扑变化 | {-e['terms']['damage_J']:.16e} | −D，仅展示，不重复累计 |
| 残差R | {e['terms']['residual_J']:.16e} | 真实未闭合部分，不改写为耗散 |

相对冻结累计能量尺度最大残差{100*e['relative_residual_max']:.9f}%。为独立识别时间层误差，令d_n下的末态力为F_old(x_(n+1))：

- `Verlet defect = ΔK_free − Σ_free ((F_n+F_(n+1))/2)·Δx`。
- `Spatial defect = U(x_(n+1),d_n)−U(x_n,d_n) + Σ_all ((F_n+F_old(x_(n+1)))/2)·Δx`。
- `Switch work = 0.5 Σ_all (F_(n+1)−F_old(x_(n+1)))·Δx`。
- 接触缺陷按相同求积定义，显式扣Jc；其本轮只剩舍入量。

三个非接触项累计分别为+0.046242415809pJ、−0.011885524449pJ、+7.226909035287pJ；负和重建R，差2.8192e-23J。主要来源为损伤后的力切换，dt/2窗口约减半；这支持有限分裂误差的解释，不表示持续时间/速率/尺度独立的物理认证。

## 分项状态

| 项目 | 本轮状态 |
|---|---|
'''
report=intro+'\n'.join('| '+k+' | '+v+' |' for k,v in items.items())+'\n\n## 51项逐项回答\n\n'
for num,(question,answer) in enumerate(answers,1):report+=f'### {num}. {question}\n\n{answer}\n\n'
report+='''## 复现、中断与停止证据

固定状态回归只使用intact、sc=.18部分、sc=.24部分、全失效四个状态，未形成第二动态族。C++最初编译后检查到LAMMPS析构时先释放force、后释放fix，增量修正析构防护，正式科学步均使用最终同一SHA。初次变量min/max语法失败有0接受步；改为数学等价ternary。第170步Python只读断言错误并非材料历史错误，INSTRUMENTATION_RECOVERY逐值证明同一接受态恢复；没有重放前170步。第244步Δd越界为真实停止条件。短窗口没有继续到更晚物理时间，没有第二窗口。

科学局限：本轮1MPI，无新增跨rank动力学迁移认证；旧2MPI静态/运动学证据保留。当前小图复制实现不证明大规模伸缩性。邻居列表本工况每步重建，不能把它写成任意邻居刷新策略都已覆盖。活动闭合接触、卸载段及P4均未覆盖。新旧二值本构差异保留，不把旧0.069696J/m²标称量用于新材料认证。

关键证据：DYNAMIC_DAMAGE_COMMIT_DESIGN_FREEZE.json、DYNAMIC_INTEGRATOR_IMPLEMENTATION_AUDIT.md、DYNAMIC_341_NODE_AUDIT.json、DYNAMIC_FORCE_ENERGY_AUDIT.json、STOP_STATE_TOPOLOGY_GRADIENT_AUDIT.json、DYNAMIC_ENERGY_LEDGER.json、DT_HALF_WINDOW_COMPARISON.json、RUN_RESOURCES.json、REPRODUCE.txt及PAYLOAD_MANIFEST.json。五组图见OPEN_RESULTS.html。

本轮开发及科学交付结束，等待下一轮科研审查；不部署生产、不运行完整V3、CFD、微泡或AI。
'''
(O/'REVIEW_LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1.md').write_text(report)
(O/'DYNAMIC_INTEGRATOR_IMPLEMENTATION_AUDIT.md').write_text('''# 原生动力学实现与实际覆盖

只修改独立的四个BPM研究C++文件；pair_bpm_peri、verlet、fix_nve、fix_move原源码SHA不变。新内部style为LPS/VARIATIONAL/DYNAMIC all rho；要求SI、正密度、全组力换算、一个free组nve、仅允许move/property与该协调fix，禁止未审查的force modifiers。实际驱动明确不重叠的49/47/245边界与自由组。复制图上限4096节点不代表该规模已认证。

## 时间层

1. 原生free组NVE首次半kick和drift；边界move以连续物理time读取冻结解析曲线。
2. POST_INTEGRATE只在whichflag=1且非setup触发，以真实ntimestep为token；首次之后必须连续+1。old/new同坐标能差记D，只在新完全失效时更新owned键和special exclusions。
3. 此时早于Neighbor::decide。设置force_reneighbor及本步next_reneighbor；当前实测输入每步完整重建，pair在bond前计算但已看到新排斥拓扑。
4. 新损伤下pair/bond力密度汇总并reverse_comm，POST_FORCE只对owned乘同一V。setup也只转换本次新计算的初始力，不提交损伤。不改变原真实mass。
5. 原生NVE最后半kick。Python只读取接受结果/求积/停止，不做自由节点位置或速度推进。

## 状态与缓存

κ/d、原始键ID、固定m/r/taper、token、计数、D、接受坐标继续使用上一轮版本化bond restart。诊断求力次数不序列化，所以进程重启后它会重置；科学commit/token不重置。force/run0/output均不提交。静态显式C API在动态fix存在时返回−7。固定同坐标新状态缓存避免紧随commit的重复局部消元；d未变无需第二次消元。仅完整全图d=1可零能零力快速返回，部分损伤星邻域没有跳过被释放的内部自由度。

## 守恒与误差

native exported force为N，U与D为J；native epair乘V，native KE使用物理质量后已为J。实际记录全体K和K_free。夹持外功为−F_material的梯形路径功加ΔK_grips；外部冲量含ΔP_grips。储能与每步D由冻结离线模型在所有245保存态独立核对。已存终态3分量固定d能量差分最大相对误差1.3468e-9。

本分裂仅保证明确且一次提交；d改变时不宣称整体二阶准确。残差分解详见DYNAMIC_ENERGY_LEDGER.json和energy_residual_decomposition.npz。主工况Δd越界时停止，有限窗口dt/2改善不构成全程通过。

## 已实际执行与未覆盖

实际244个主物理步、一次同接受态restart续算、一次68步短dt/2、固定状态与只读试探。完整失效原键t=0；自然轨迹未闭合接触，未运行P4、2MPI动态迁移或生产规模。历史静态接触和2MPI证据不替代本轮未覆盖项。析构函数防护确保LAMMPS force先释放时不解引用悬挂material指针；正式测试前已修正。
''')
checks=dict(P0_before_CPP=True,compile=True,frozen_state_regression=json.loads((O/'P2_FIXED_STATE_REGRESSION.json').read_text())['status'],read_only_repair_recorded=True,source_mathematics_unchanged=True,physical_force_and_mass_verified=True,all_main_steps_independently_reconstructed=True,original_native_endpoint_topology_verified=True,history_run0_and_restart_verified=True,force_gradient_relative=1.3468153377642423e-9,dt_safety_main='STOPPED',P4='NOT_RUN',active_closure='NOT_OBSERVED',MPI_dynamic='NOT_RUN',required_51_answers=len(answers),all_five_figures_present=all(len(list(F.glob('*.'+ext)))==5 for ext in ['png','svg']),documentation=doc['status'],production_protection=protect['status'])
dump(O/'IMPLEMENTATION_CHECKS.json',checks)
# Preserve P0 index, then update this round's own entrypoint with final scope.
if not (O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS_P0.md').exists():shutil.copy2(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS_P0.md')
entry=f'''# 当前clot原生软化开发入口

本轮：{TASK} / {RUN}。

- **结论：P3第244步(122µs)因最大Δd=0.102554>0.10停止；P4 NOT_RUN。**
- 已证明：原生VV真实自由积分、一次提交、自然软化/全失效、当前窗口力/历史/动量一致，能量余额0.144631%。
- 唯一dt/2窗口105–122µs局部改善；没有完整时间步认证或后续自动续算。
- 先读本轮 `verification/{TASK}/{RUN}/FINAL_STATUS.json` 与51项REVIEW报告。
- 服务器环境：`source /workspace/y_clot/env.sh`；项目 `/workspace/y_clot/project`。
- 当前C++：`analysis/dynamic_softening_v1/bond_bpm_peri_variational.h/.cpp`。
- 动态协调：同目录 `fix_lps_variational_dynamic.h/.cpp`（内部style LPS/VARIATIONAL/DYNAMIC）。
- 实际隔离库：`{LIB}`。
- 该库SHA：`{sha(LIB)}`；仅研究使用，未写入env.sh或生产入口。
- 冻结数学：`analysis/softening_feasibility_v1/model.py`，SHA `{sha(R/'analysis/softening_feasibility_v1/model.py')}`。
- 动态驱动：`analysis/dynamic_softening_v1/native.py`、`dynamic_run.py`；读事件/能量用 `postprocess.py`。
- 本轮源快照/补丁/复现步骤：verification本轮的source_snapshot、NATIVE_SOURCE_PATCH.diff、REPRODUCE.txt。
- 核心FP64：main_accepted_trajectory.npz；dt_half_window_trajectory.npz；energy_residual_decomposition.npz。
- 只读证据：DYNAMIC_FORCE_ENERGY_AUDIT、STOP_STATE_TOPOLOGY_GRADIENT_AUDIT、RUN0_AND_TRIAL_STATE_SAFETY。
- 图：results本轮OPEN_RESULTS.html（5组PNG/SVG，P3实际停止点云，无P4假图）。
- 参考341点/12047原键；左右49/47受控，245自由，dx50µm，horizon150.5µm。
- G940Pa，K100000/3Pa，rho1000，s0=.15，主sc=.18；完整m/r/taper不变。
- native BPM f为N/m³，新协调fix在POST_FORCE乘V，原真实质量rhoV保持；导出force/U为N/J。
- 只在真实POST_INTEGRATE提交，run0和试探不提交；accepted token在restart延续。
- 上轮受控运动学基线：verification/lps_variational_softening_native_isolated_v1/20261013T201800Z。
- 上轮库SHA b80bfc0a99f4160d62567c0e2d9ddc611004ad40f79c8e49370293c681560773，保留。
- 冻结新旧二值本构力L2差8.814180703270e-6N保留，不靠经验修正消除。
- 本轮未覆盖：卸载段、活动闭合接触、P4、2MPI动态迁移；宏观强度/Gc/生产伸缩性未认证。
- 生产库 tools/lammps-native-y-wall-shared-v1-build/liblammps.so SHA d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a，不变。
- P0审计28候选，安全删除0；UNCERTAIN及历史科研证据均保留。31001历史文件保护通过。
- 实际WSL文档：/home/lzy/projects/clot_rupture.md；before/after/diff/SHA见本轮verification。
- 下一轮唯一优先：同一341点冻结工况的损伤时间步准入与全时程能量验证；先预注册，不自动运行。
- 不必展开旧PROJECT_ACTIVE_CONTEXT及全部历史目录；按需读取以上少量入口即可。
'''
assert len(entry.splitlines())<=150;(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md').write_text(entry);(O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md').write_text(entry)
(O/'REPRODUCE.txt').write_text(f'''EVIDENCE-FIRST REPRODUCTION; no automatic new scientific run.
1. Read FINAL_STATUS.json, design freeze, this round's source patch and report.
2. Actual environment: source /workspace/y_clot/env.sh
   Set LD_LIBRARY_PATH={O}/build:$LD_LIBRARY_PATH
   Set PYTHONPATH={O}/source/python:$PYTHONPATH
   Original production paths and env.sh remain unchanged.
3. Actual source: {O}/source; exact four BPM C++ files in source_snapshot and NATIVE_SOURCE_PATCH.diff.
   Upstream base is prior accepted source {PREV}/source. Do not overwrite it.
   Rebuild only in a NEW isolated directory with the CMake command saved in BUILD_RESOURCES.json.
   The final patch includes src/.gitignore bookkeeping; C++ used in all dynamic runs is unchanged.
   Do not rerun p0_cleanup.py/p1_freeze.py or overwrite the sealed run directory.
4. Actual scientific execution sequence:
   P0 p0_cleanup.py -> P1 p1_freeze.py -> build.py -> p2_checks.py.
   dynamic_run.py used final library SHA {sha(LIB)}.
   Zero-step input min/max syntax failure archived under setup_attempt_min_syntax.
   main trajectory step0..170 in dynamic341_trajectory.npz.
   Trial assertion repaired (history invariant, U changes under changed x); continue_main.py resumes saved step170 without physical replay.
   continuation trajectory step170..244; stop at registered max-damage increment.
   half_window.py runs only saved step210 at dt/2 for68 steps, ending same physical122us.
   postprocess.py merges the same-case segments and independently reconstructs every main state.
   stop_state_checks.py reloads frozen native endpoint topologies and does fixed-history energy finite differences; no new dynamic steps.
5. Raw FP64 main_accepted_trajectory.npz is the combined primary evidence. Physical time is the time_s column; dt/2 integer step IDs differ.
   d/kappa/force/x/v remain full precision. Complete failure event IDs in FIRST_SOFTENING_AND_FAILURE_EVENTS.json and segment EVENTS files.
6. For a future authorized replay, choose a NEW run_id and adjust only output paths in a private copy of common.py;
   preserve scientific design/material/input and point to immutable reference files. Existing scripts deliberately refuse overwriting completed main runs.
   Do not silently bypass P0 active-index conflict guards; preserve existing active index and make a separate reproduction sandbox.
7. Payload verification (read-only): compare each relative path to SHA256 in PAYLOAD_MANIFEST.json.
   Server retains the full source/build tree; local package omits compiled object caches and duplicate upstream source, retaining actual binary and touched sources/core hooks.
   Frozen offline model.py and original geometry are copied into source_snapshot/provenance for dependency clarity.
8. No P4, no dt/4, no full V3, no production deployment. A fresh longer scientific protocol needs the next research instruction.
''')
dump(O/'DELIVERY_READINESS.json',dict(status='READY_FOR_MANIFEST',report_questions=51,figures=5,P4='NOT_RUN',doc_readback_verified=True,active_entry_lines=len(entry.splitlines()),final_library_sha256=sha(LIB)))
print('REPORT_51_AND_FINAL_STATUS_WRITTEN',flush=True)
