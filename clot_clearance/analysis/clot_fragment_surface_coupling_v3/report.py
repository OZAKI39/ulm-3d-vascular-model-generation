"""Evidence-driven Chinese review and portable static results entry."""
from common import *
from analyze import read_history
import csv,html

def main():
    result=json.loads((O/'NATIVE_V3_RESULT.json').read_text());q=result['objective'];status=json.loads((O/'FINAL_STATUS.json').read_text());lineage=json.loads((O/'FRAGMENT_LINEAGE.json').read_text())
    old=json.loads((V2_O/'NATIVE_BIPHASIC_RESULT.json').read_text());u=json.loads((V2_O/'NATIVE_UNIFORM_RESULT.json').read_text());native=result['native'];last=result['last_row'];cost=json.loads((O/'COMPUTE_COST.json').read_text())
    audit=json.loads((O/'SURFACE_GEOMETRY_AND_MAPPING_AUDIT.json').read_text());cover=json.loads((O/'SURFACE_COVERAGE_DEFINITIONS.json').read_text())['targets'];p0=json.loads((O/'P0_TARGET_SURFACE_DIAGNOSIS.json').read_text())['targets']
    def pc(x):return '未定义' if x is None else f'{100*x:.3f}%'
    targetcount=int(q['n_fragments_20_200']);outcome='末态达到了至少 5 个目标块的节点数条件' if targetcount>=5 else '末态没有达到至少 5 个目标块的节点数条件'
    motion=json.loads((O/'SURFACE_MOTION_COMPARISON.json').read_text());same_time=motion['stages'][-1]
    motiontext=f"相同物理时刻 {same_time['time_s']*1e3:.6f} ms，新旧 H 的全部 8125 个原始材料节点位置 RMS 差为 {same_time['all_material_RMS_position_difference_m']*1e6:.6g} µm，速度 RMS 差为 {same_time['all_material_RMS_velocity_difference_m_s']:.6g} m/s。另按旧末态三个目标块及新末态目标块的原始成员集合，分别保存新旧轨迹中的质心、速度、动能、内部相对动能、实际节点载荷和真实分量归属。这些固定成员组在另一个版本或另一个时刻可能已分裂，不能当作同一个完整碎块。详见 SURFACE_MOTION_COMPARISON.json。"
    text=f'''# CLOT_FRAGMENT_SURFACE_COUPLING_AND_PHYSICAL_HETEROGENEITY_V3

唯一正式 H 算例状态：**{native['status']}**，最后接受时间 **{native['time_s']*1e3:.6f} ms**，目标区间 20–200 节点的末态游离块数 **{targetcount}**，时程最大数 **{q['max_target_count']}**。{outcome}。细屑、形貌、撤载后存活和表面完整性继续分别报告。

科学标签：`EXPLORATORY_UNCALIBRATED_BIPHASIC_MATERIAL`、`PRESCRIBED_SHEAR_MECHANICAL_SURROGATE`、`CURRENT_SURFACE_APPROXIMATE_COUPLING`、`MICROSTRUCTURE_SCALE_UNCALIBRATED`。

## 1. 五个零加载目标碎块究竟为什么没有面

直接读取 V2 FP64 末态，没有重新执行旧动力学。U 1164、U 1804、H 3809、H 3846 的全部存活原始材料图仍各自连通，但没有任何完整存活、正体积的辅助四面体，属于原因 A。U 1848 有 5 个这样的四面体，但 16 个候选面全未通过开放判定，属于原因 B。U 1884 的 12 个候选中只有 2 个加载；H 3769 的 22 个候选中只有 1 个加载，后者另有一个已开放面未通过外部流体连通检查。不能把七块都归咎于“没有四面体”。

未发现这七块的零加载是节点 ID／所有权映射错误造成的；旧面公式对原始节点力的回放通过。P0 记录也没有支持“当前遮挡”或退化法向是这五个零面的主因。完整原始节点、参考／当前坐标、内部存活键、失效外部键和实际候选面均在 `p0_targets/`；`P0_TARGET_SURFACE_DIAGNOSIS.json` 保存逐类证据。

## 2. 新表面求积有哪些物理与离散支持

保留旧 L1 受支持三角面路径。L2 使用原始 dx³ 对偶控制体积以及存活原始非局部邻域的仿射几何，补充真正游离分量的参考边界面。近邻键失效不等于裂面开放：所有原始键穿过对偶面的参考交叉表、当前分量身份、开口、当前单元障碍与真实 L1 四面体遮挡都参与筛选。旧 L1 面附近的补充单元保守排除，以减少重复面积。没有凸包外壳、人工连边、键恢复或按目标尺寸切分材料。

这是局部控制体积的**部分近似表面求积**，不是经过验证的连续流体包络。每片采用中心一点求积；可达性是外向直线与局部探针的代理，不等同于真实曲折流路，亦未逐点认证整个面片的连续润湿。秩亏／低维碎屑、桥接、几何反转、外推过大、遮挡和重叠歧义保留为未解析。真实遗漏的连续表面积仍未知。

同一 V2 七块的离线回放如下。有效面积包含开口权重；两种覆盖率的分母均为该碎块节点数，均不是连续表面积覆盖率。

|碎块|节点|L1 面|补充 L2 面|求积所有者／顶点覆盖|受力权重支撑覆盖|有效面积 mm²|
|---|---:|---:|---:|---:|---:|---:|
'''
    for c in cover:text+=f'|{c["case"]} {c["component_id"]}|{c["nodes"]}|{c["L1_active_facets"]}|{c["L2_active_dual_faces"]}|{pc(c["quadrature_owner_or_vertex_fraction"])}|{pc(c["mapped_force_support_fraction"])}|{c["effective_loaded_area_m2"]*1e6:.6g}|\n'
    maxforce=max(c['mapping']['force_N'] for c in audit['checks']);maxmoment=max(c['mapping']['moment_Nm'] for c in audit['checks']);maxpower=max(c['mapping']['power_W'] for c in audit['checks'])
    text+=f'''
## 3. 力、力矩、功率误差及限制

对同一面使用同一组仿射权重插值中心与速度，并伴随分配节点力。七块所在真实状态上的最大合力误差 {maxforce:.6g} N、力矩误差 {maxmoment:.6g} N·m、功率配对误差 {maxpower:.6g} W。恒定 1000 Pa 剪切、刚体平移阻力及随动旋转检查通过。

L1 权重非负。L2 外侧面中心常在节点支撑凸包之外，强制非负标量权重会破坏该中心的力矩／虚功配对，故采用有仿射意义的受限带符号权重：最低 −0.25、绝对权重和不超过 2；超过限制不激活。这控制了外推，**没有认证局部牵引分布精度或网格收敛**。

从真实 V2 末态进行最小原生微步验证，含新 L2 非零阻力、只接受态提交、非提交试态和一段重启续步。最后版本的 2 个连续微步与 1 个重启续步通过；重启后的坐标、速度、力、κ/d 和物理能量标量与连续路径相同。开发时因增加真实 L1 障碍校核和记录字段，对这组微验证进行了必要复核，没有扩大物理试验矩阵。正式运行的每步门限仍为 Δd≤0.1、相对能量余额≤1%，dt 基准 0.25 µs、最小 0.03125 µs。

## 4. 新 H 的运动、二次碎裂与撤载状态

|指标|旧 V2 H|新 V3 H|
|---|---:|---:|
|末态 20–200 节点游离块|{old['objective']['n_fragments_20_200']}|{targetcount}|
|最大游离块节点|{max(old['final_free_sizes'])}|{max(result['final_free_sizes'],default=0)}|
|目标质量／初始质量|{pc(old['objective']['mass_fraction_target_initial'])}|{pc(q['mass_fraction_target_initial'])}|
|1–2 节点细屑／游离质量|{pc(old['objective']['fine_debris_mass_fraction_detached'])}|{pc(q['fine_debris_mass_fraction_detached'])}|
|游离质量／初始质量|{pc(old['objective']['total_detached_mass_fraction'])}|{pc(q['total_detached_mass_fraction'])}|

新 H 保存 {len(native['secondary_events'])} 次原始成员分裂事件。其中至少两个子分量均有 ≥20 节点的描述性子集有 {len(lineage['cohesive_secondary_splits'])} 次；这包含 >200 节点子块，不改变 20–200 的目标定义。事件前后内存活键、实际表面力／功率和成员闭合记录保存在 `FRAGMENT_LINEAGE.json`。

2.8 ms 撤去外部驱动，仍按真实运动保留被动阻力。撤载时的目标块中，末态原成员集合完整不变者为 {q['exact_target_sets_surviving_drive_off']} 个；这些目标块的末态目标区间后代共 {q['target_descendants_of_drive_off_targets']} 个。掉点会改变精确成员 ID，不能据此把仍连通的较小后代判为完全消失。逐块最大后代保留率、内存活键、运动、剪切和功率均列入谱系文件。

同一 H 的新旧差异与表面／流体载荷表示变化一致，但单次实现不足以证明特定二次断裂由阻力单独造成。对比旧 Uniform 时，材料和表面版本均不同，不是双相材料的独立因果对照。低维碎屑和部分凝聚体的流体表面仍未解析，不能把零数值载荷写成零物理阻力，也不能仅凭撤载后图连通宣称真实水动力稳定。

{motiontext}

补充表面尚未起作用的前缀核对保存在 PRE_SUPPLEMENT_PREFIX_COMPARISON.json：初态逐位相同，首次软化和首次自然失效时刻分别仍为 0.5505 ms 和 0.591 ms；对应保存态的坐标最大差低于 1×10⁻¹⁸ m、速度最大差低于 6×10⁻¹⁴ m/s。新旧 L1 插值求和次序存在舍入差异，因此没有声称整个前缀逐位相同。

## 5. 至少五个目标块是否出现

末态 {targetcount} 个，时程最多 {q['max_target_count']} 个；目标区间包含 20 和 200 两端。撤载观察窗最低目标数量为 {q['minimum_target_count_in_drive_off']}。{outcome}。全部计数来自存活原始 bulk/interface 材料图与锚层可达性，空间接近、接触、相区标签和视频颜色均不增加材料边。

`FRAGMENT_MORPHOLOGY.csv` 分别给出参考／当前包围盒、回转半径、主轴、长宽比、内部键密度和细长／片状特征。三维点云秩或节点数不是紧实连续固体的认证。

## 6. 材料相关长度是否有独立证据

现有资料支持真实组织的异质性，但不足以标定当前强相区域的相关长度。Liu 的 0.1–1 mm 是断裂过程区／非线性尺度；2022 年 microCT/SEM 研究描述数百微米级 RBC 聚集区及数十微米多孔边界，不能直接当作本模型强纤维蛋白／血小板区的已标定直径。4–7 µm 单 RBC 尺度也不是相区相关长度。

因此正式 H 原样沿用 10 强区、270–271 节点／区、强相体积分数 33.32923%。体积等效直径约 401 µm，与中心跨度 250–350 µm 分别报告。原始跨相键 91576，跨不同强区直接非局部键 3237（初始存活 3004），说明强区并非相互独立的预设碎块。弱基质最小参考单元间隙约 50 µm，小于 horizon；它是离散几何估计，不是实测组织隔层厚度。

新增参数化生成器以米制相关长度／各向异性、强相概率和固定 seed 为输入，没有目标块数或目标节点数参数。ell=100/200/400 µm 是明确未标定的离线假设，分别检查相同连续随机场在 dx=25/50/100 µm 的取样差异；不把 ell 当成团直径，没有按碎裂结果选择 seed，也未运行任何新材料尺度动力学。

资料与测量对象详见 `PHYSICAL_DOMAIN_SCALE_ASSESSMENT.md` 和 `DOMAIN_SCALE_PROVENANCE.json`。外部原始来源：[2022 microCT/SEM](https://www.nature.com/articles/s41598-022-06623-8)；表面牵引离散的一致性背景：[Yu–You–Trask 2021](https://arxiv.org/abs/2101.01434)、[Yu–Huang–Zhou 2023](https://doi.org/10.3390/ma16062252)。本实现没有声称继承这些论文的全部模型或收敛证明。

## 7. 下一步优先事项

先依据本轮真实后代谱系与未解析面清单，改善仍被桥接／几何不确定性挡住的表面表示，并验证局部载荷分配、分辨率与遮挡的敏感性；不能把受力权重扩散到更多节点等同于表面已完整。再获取与试样对应的三维组织分区及相关尺度，独立校准强弱相损伤／断裂响应。当前数据不足以直接指定“缩小强区到某个目标节点数”。

声流对接仍需要可信流体速度／牵引来源及真实近壁几何，而非把名义 1000 Pa 直接解释为微泡能力。本轮没有新增材料／幅值／随机种子扫描，没有新 Uniform 全时程，没有完整 Navier–Stokes 双向 FSI。

## 完整能量与成本

Wdrive + Wdrag = ΔU + ΔK + D + ΔUc − Jc + R。末态驱动功 {last['W_drive_J']*1e9:.6g} nJ，阻力功 {last['W_drag_J']*1e9:.6g} nJ，最大归一化余额 {native['max_relative_energy_residual']*100:.6g}%。R 单列，没有并入 D。L2 与可达新裂面的功是可能重叠的诊断子通道，未重复加入总功。表面拓扑变化时仍以真实节点位移积分，避免虚拟求积中心跳变造成伪功。

正式原生接受步 {native['accepted_steps']}，预测拒绝 {native['predictive_rejections']}，CPU {cost['cpu_s']/3600:.4f} h，墙钟 {cost['wall_s']/3600:.4f} h，峰值 RSS {cost['peak_RSS_bytes']/2**30:.4f} GiB。只有一个新正式 H，无墙钟或 4 GiB 上限。

`FRAGMENT_FORCE_AND_POWER.csv` 是每接受态的 **20–200 节点目标块**明细；所有尺寸分量的完整逐态受力、功率、面积、累计功、未解析标记仍保存在 `cases/biphasic/CONNECTED_COMPONENT_HISTORY.jsonl.gz`，原始成员集合在 `COMPONENT_REGISTRY.json`。累计成员功从 t=0 沿原始节点继承，精确成员集合出生后的增量另列，避免分裂时重置／重复分配历史功。

源码／输入冻结：`SOURCE_FREEZE.json`。原生库继续使用 V2，SHA256 `{sha(LIB)}`。没有修改历史 V2 文件。没有删除用途不明的历史材料或代码。原生状态、日志、重启、全部诊断与图／视频都保留；本轮 Git 提交只选取代码、必要依赖、紧凑证据和可复现记录，不打包无关大文件。
'''
    report=O/'REVIEW_CLOT_FRAGMENT_SURFACE_COUPLING_V3.md';report.write_text(text)
    reproduction=f'''# Reproduce / replay\n\nTask {TASK}, run {RUN}. One formal H attempt only.\n\nActivate `/workspace/y_clot/env.sh`; set OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1. Runtime root defaults to `/workspace/y_clot/project`; `Y_CLOT_PROJECT_ROOT` overrides it.\n\nRead-only replay, with full archived data present:\n\n```bash\npython analysis/{TASK}/analyze.py\npython analysis/{TASK}/coverage_audit.py\npython analysis/{TASK}/figures.py\npython analysis/{TASK}/render_videos.py --smoke\n# Review the five actual-frame smoke previews before final video rendering.\npython analysis/{TASK}/render_videos.py\npython analysis/{TASK}/report.py\n```\n\nThese commands read accepted dynamics and do not restart the physical experiment. The launch script deliberately refuses a second formal attempt in the frozen directory. A justified new experiment requires a new independent RUN and a new frozen request.\n\nV2 dependency root: `{V2_O.relative_to(R)}`. Required original geometry / material / restart / states and native library hashes are in SOURCE_FREEZE.json and V2_DEPENDENCY_FREEZE.json. The raw library, dense trajectories, large native restart/data and paper PDFs are archival dependencies, not incidental Git payload. The native library SHA is `{sha(LIB)}`. V2 BUILD_AUDIT.json records compiler/object provenance; upstream LAMMPS commit is `8de817dd79bfe4525d5d39246a212d833e6dee07`. This round changes no LAMMPS material/force kernel. The small ray helper is compiled from ray_cells.cpp with g++ -O3 -std=c++17 -fPIC -shared and its source SHA is recorded.\n\nThe binary restart is accompanied by surface_restart.npz and fluid_work_restart.npz; reference ID/time/step checks are mandatory. A failed-unaccepted native restart, if present, is explicitly not a checkpoint of last_accepted.npz. Continuation equivalence of the new surface sidecar is tested in NATIVE_SURFACE_HISTORY_RESTART_AUDIT.json.\n\nCoordinates are metres in FP64 archives; figures/videos display micrometres without displacement magnification. Frame-time CSV files point to actual accepted native keys. No synthetic trajectory interpolation.\n\nThe physical-scale candidates are offline random fields with fixed physical Fourier coefficients, never formal dynamics. The active material remains the frozen V2 BIPHASIC_FIELD.npz.\n'''
    (O/'REPRODUCE.md').write_text(reproduction)
    rel=lambda p:os.path.relpath(p,F)
    links=[('中文科研审查',report),('复现说明',O/'REPRODUCE.md'),('最终状态',O/'FINAL_STATUS.json'),('表面方法',O/'SURFACE_METHOD_SPEC.md'),('物理尺度证据',O/'PHYSICAL_DOMAIN_SCALE_ASSESSMENT.md'),('原生验证',O/'DELIVERY_VALIDATION.json'),('逐态尺寸统计',O/'FRAGMENT_SIZE_HISTORY.csv'),('目标块受力／功率',O/'FRAGMENT_FORCE_AND_POWER.csv'),('完整能量账本',O/'ENERGY_LEDGER.csv'),('碎块谱系',O/'FRAGMENT_LINEAGE.json'),('同成员运动比较',O/'SURFACE_MOTION_COMPARISON.json'),('原始数据范围',O/'DATA_INVENTORY.json')]
    pieces=['<!doctype html><html lang="zh"><meta charset="utf-8"><title>Clot surface coupling V3</title><style>body{max-width:1360px;margin:32px auto;padding:0 24px;font:18px/1.65 system-ui;color:#243d4c;background:#f8fafb}img,video{width:100%;height:auto;background:white}section{background:white;margin:28px 0;padding:20px;border:1px solid #dce5e9}a{color:#126b89}table{border-collapse:collapse}td,th{padding:8px 18px;border:1px solid #cdd9df}.tag{color:#9b4b52}</style>',f'<h1>碎块表面耦合 V3</h1><p>{html.escape(outcome)}；末态 {targetcount} 块，时程最大 {q["max_target_count"]} 块。原生状态：{native["status"]}。</p><p class="tag">探索性未标定材料 · 规定剪切机械代理 · 部分近似表面耦合</p>','<p>'+'　'.join(f'<a href="{html.escape(rel(p))}">{label}</a>' for label,p in links)+'</p>']
    for p in sorted((F/'figures').glob('FIG_*.png')):pieces.append(f'<section><h2>{html.escape(p.stem)}</h2><a href="figures/{p.with_suffix(".pdf").name}">矢量 PDF</a><img loading="lazy" src="figures/{p.name}" alt="{p.stem}"></section>')
    for p in sorted((F/'videos').glob('VIDEO_*.mp4')):pieces.append(f'<section><h2>{html.escape(p.stem)}</h2><video controls preload="metadata" src="videos/{p.name}"></video><p>1920×1080 · 24 fps · 原生接受态真实坐标；具体时间见 <a href="videos/{p.stem}_FRAME_TIMES.csv">帧来源表</a>。</p></section>')
    pieces.append('</html>');(F/'OPEN_RESULTS.html').write_text('\n'.join(pieces))
    print('REVIEW_AND_HTML_WRITTEN',flush=True)

if __name__=='__main__':main()
