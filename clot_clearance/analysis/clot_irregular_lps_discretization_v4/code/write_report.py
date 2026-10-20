"""Evidence-driven final report and static index, only after the trajectory ends."""
from common import *
import html,shutil,subprocess

def main():
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE'];outcome=json.loads((O/'SCIENTIFIC_OUTCOME.json').read_text());geo=json.loads((O/'IRREGULAR_REFERENCE_GEOMETRY.json').read_text());elastic=json.loads((O/'ELASTIC_PATCH_CHECK.json').read_text());perf=json.loads((O/'PERFORMANCE_OPTIMIZATION.json').read_text());base=json.loads((O/'PERFORMANCE_BASELINE_V3_COMPACT.json').read_text());native=json.loads((O/'NATIVE_STATE_SURFACE_CHECK.json').read_text());crack=json.loads((O/'CRACK_METRICS.json').read_text());survival=json.loads((O/'FRAGMENT_SURVIVAL.json').read_text());complete=result['status']=='COMPLETED';p=perf['comparisons']['damaged'];ncheck=native['history']
 timing='\n'.join(f'| {name} | {v["wall_s"]:.6f} | {v["cpu_s"]:.6f} | {v["calls"]} | {v["mean_wall_s"]:.6g} | {100*v["wall_fraction"]:.2f}% |' for name,v in sorted(base['windows']['damaged']['python_exclusive'].items(),key=lambda x:-x[1]['wall_s']))
 elastic_rows='\n'.join(f'| {r["mode"]} | {r["modulus_energy_Pa"]:.6f} | {s["modulus_energy_Pa"]:.6f} | {100*elastic["relative_modulus_changes"][r["mode"]]:+.4f}% |' for r,s in zip(elastic['patches']['regular']['rows'],elastic['patches']['irregular']['rows']))
 physics_status='完整单次时程已完成。' if complete else '正式轨迹触发冻结数值门槛后停止；未追加物理算例、改阈值或把失败态标为接受态。'
 scientific='只能依据这一组参考点云，评价整个离散方案的变化；当前结果不支持统计稳健性、真实微观结构或声流机制认证。'
 failure='无。' if not result['exception'] else '```text\n'+result['exception']+'\n```'
 pilot_path=O/'cases/H_IRREGULAR_LPS_V4/PILOT_RESULT.json'
 if pilot_path.exists():
  pilot=json.loads(pilot_path.read_text());pilot_note=f"2.0ms试运行状态{pilot['status']}，wall={pilot['wall_s']:.3f}s、CPU={pilot['cpu_s']:.3f}s、接受步{pilot['accepted_steps']}。它仅是数值准入窗口。按这一窗口平均成本线性外推3.6ms约{pilot['wall_s']*1.8/3600:.3f}h；这只是当时的成本尺度，后续损伤/表面复杂度不同，实际总成本以上述终态计时为准。后段读取同一pilot checkpoint延续，没有重算前2ms。"
 else:pilot_note='未完成2.0ms数值试运行窗口，因此没有准入后段正式续算，也不从这一短时窗口推断最终碎块目标失败。'
 text=f'''# REVIEW_CLOT_IRREGULAR_LPS_V4

任务：V4 不规则 LPS 点云、性能优化与碎块稳定性研究。Run `{RUN}`。

状态 **{result['status']}**；合法接受到 **{result['accepted_time_s']*1e3:.8f} ms / 3.6 ms**。{physics_status}

本轮对象为8125个材料点、1875个固定锚点的标准试样。模型标记 **EXPLORATORY_UNCALIBRATED**；加载标记 **MECHANICAL_SHEAR_SURROGATE**。

## 1. 冻结来源和研究边界

冻结本地V3 HEAD `43ee28e7719e26189757ea472ee52752b63fa981`；服务器Git记录属于旧迁移元数据，实际源码逐文件与最新V3清单核对。详见 `SOURCE_FREEZE.json`。正式V4实现、材料、参考和库由 `V4_IMPLEMENTATION_FREEZE.json` 锁定。原V3正式3.6ms数据只读取，没有重跑旧完整基线。

用户明确选择保持V3积分体积和质量：材料点中心包络1.2×1.2×0.6mm = 0.864mm³；积分域1.25×1.25×0.65mm = **1.015625mm³**；材料质量 **1.015625mg**。锚点体积另列，不计入血块质量。两个体积口径不能混写。

G=940Pa，K=100000/3Pa，ρ=1000kg/m³；δ=150.5µm，名义dx=50µm。弱/混合/强相阈值分别(0.100,0.120)、(0.175,0.210)、(0.250,0.300)，原10强区及全部原材料ID相区不变。相同预制缺口物理切断规则与1000Pa名义当前表面切向驱动：2.0ms升载结束、2.4ms降载开始、2.8ms驱动归零、3.6ms观察结束；被动相对速度阻力保留。

## 2. Q1：V3哪些模块最耗CPU？

完整态与真实2.650125ms损伤态各取8个原生接受步。损伤窗口原总wall={p['before_wall_s']:.6f}s、CPU={p['before_cpu_s']:.6f}s。Python独占分项如下；内部C++计时是对应调用的细分，不能再次加到总时间。

| 独占分项 | wall s | CPU s | 调用数 | 平均wall s | 占总wall |
|---|---:|---:|---:|---:|---:|
{timing}

最大热点为碎块谱系/诊断，随后为磁盘checkpoint和首次事件输出。原生力/积分在这一代表性窗口不是最大瓶颈。完整原生细分、调用次数、均值保留于 `PERFORMANCE_BASELINE_V3.json`。这是单次配对计时，没有统计误差区间。

## 3. Q2：保持物理语义的性能提升

损伤窗口 **{p['before_wall_s']:.6f}→{p['after_wall_s']:.6f}s，{p['speedup']:.4f}倍**。CPU={p['before_cpu_s']:.6f}→{p['after_cpu_s']:.6f}s；峰值RSS={p['before_peak_RSS_bytes']/2**20:.2f}→{p['after_peak_RSS_bytes']/2**20:.2f}MiB。8个接受步和1次预测拒绝保持。

把“旧分量对全部新分量求集合交集”改为原成员ID直接查当前标签；同时测试严格按相同参考查询和当前坐标/实体四面体集合缓存。初末物理数组、全量标量CSV及父子事件完全一致。完整态窗口 **{perf['comparisons']['intact']['before_wall_s']:.6f}→{perf['comparisons']['intact']['after_wall_s']:.6f}s**，没有测得加速。

不把短窗1.55倍当作原11.7小时的已实测整体提升。V4正式入口每次损伤拓扑改变均更新真实成员与分裂事件，完整形貌记录降为25µs及协议端点；受力表面始终跟随当前原生坐标/拓扑，不按每N步延迟更新。预测拒绝发生在原生步之前，不会产生需撤销的损伤提交；静态试算和表面缓存回退另经验证。

## 4. 不规则参考与积分

种子20261019；初始提案为三维均匀球内随机位移，模长上限12.5µm。固定外边界、预缺口边界和锚点；通过很小的几何投影保留原强相积分体积。实际分布受这些约束条件化，不能称为无条件独立均匀随机点云。

- 实际移动 **{geo['displaced_nodes']}** 个材料点；最大位移{geo['max_displacement_m']*1e6:.6f}µm；最小间距{geo['minimum_point_spacing_m']*1e6:.6f}µm。
- 强相节点2708个，强相体积分数{100*geo['strong_volume_fraction']:.8f}%；投影最大修正{geo['projection_max_correction_m']*1e6:.6f}µm。
- 每个映射控制体用6个辅助四面体积分，全部正体积；最小有向体积比{geo['minimum_tet_signed_volume_ratio']:.6f}。
- 材料体积范围{geo['nodal_volume_min_max_m3'][0]:.9e}～{geo['nodal_volume_min_max_m3'][1]:.9e}m³；没有把移动后的点继续当作等体积点，也没有事后统一缩放体积凑总量。
- KD-tree按新参考位置重建原始材料图：**{geo['original_bonds']}→{geo['new_reference_bonds']}**条；新增{geo['added_bonds']}、移除{geo['removed_bonds']}；预缺口切断{geo['initial_cut_bonds']}条；初始图连通、无非预期孤点、无初始未连接接触对。
- 所有键通过端点原始ID对匹配；新数字bond ID不冒充旧bond ID。断后不按当前位置补键。

完整endpoint-star储能沿实际节点体积改写：`Ai=Vi*(K-5G/3)`，`Lib=3*Vj*hb/mi`，`Dib=15*G*Vi*Vj*hb/(rb*mi)`，`kib=9*K*Vi*Vj*hb/(rb*mi)`；`mi=sum(hb*rb*Vj)`包含完整原始邻域。维持完整局部秩一变分消元，材料力是该能量的负梯度。原生内部力和势能统一为N/J；nve质量为逐ID的ρVi，不再乘全局体积。

接触保留原线性排斥、50µm起始距离、存活材料键排除和对称性；非等体积时把原平均体积约定延伸为`((Vi+Vj)/2)^2`物理力/势系数，等体积极限还原原V3。它不是接触精度或真实纤维直径的认证。短参考键断后可能激活50µm接触，新增接触势Jc独立记账。

L1支撑四面体仍由原格点身份拓扑定义，但其坐标、面积、法向、体积与键身份均使用新参考。L2改为真实映射控制体的面中心、面积向量、原始全键跨面列表和六四面体射线遮挡，映射权重使用完整三维参考偏移；当前位置仍为局部仿射近似。力、力矩、功率配对守恒。L2依然是保守的部分受力面，不等于连续裂面或完整水动力润湿。

## 5. Q3：弹性与原生状态是否合理？

729节点同成员小试样，全部坐标按ε=10⁻⁴体积应变或γ=10⁻⁴简单剪切规定。不是自由边界宏观试验，也没有重拟合K/G。

| 模式 | V3能量响应Pa | V4能量响应Pa | V4相对变化 |
|---|---:|---:|---:|
{elastic_rows}

不规则剪切方向间spread={100*elastic['shear_directional_spread']['irregular']:.5f}%，通过预先冻结的15%比较门槛。规则小试样的表观剪切响应约1317Pa，本身不等于材料卡940Pa；有限邻域、边界截断和全场仿射约束会影响这一诊断。这里证明的是相对离散变化受控，不能把它写成宏观剪切模量已标定。

原生与独立力/能量、能量有限差分均通过；完整不规则参考态内力/储能/接触势为零。软化验证中真实提交{ncheck['softened_bonds']}条软化键，最大d={ncheck['max_damage']:.8f}；固定损伤的梯度相对误差{ncheck['softened_force_gradient_relative_error']:.3e}。代表性历史回退及restart续步物理状态误差为零或机器舍入量级；非均匀接触势及梯度验证通过。仅映射测试使用人工分离夹具，其数据完全排除正式碎裂结果。详细见 `LPS_DISCRETIZATION_CHECK.json`、`NATIVE_STATE_SURFACE_CHECK.json`。

## 6. 正式接受轨迹与账本

唯一工况 `H_IRREGULAR_LPS_V4`。接受步{result['accepted_steps']}，预测拒绝{result['rejected_predictions']}，接受细分步{result['substeps']}。累计wall={result['wall_s']:.6f}s、CPU={result['cpu_s']:.6f}s、峰值RSS={result['peak_RSS_bytes']/2**30:.4f}GiB。未设置4GiB或1800s预算。

{pilot_note}

原协议BASE=0.25µs、最小事件步0.03125µs、单步Δd≤0.1。外功沿真实接受端点力与真实位移累计，Wdrive/Wdrag分列；完整账本`W=ΔU+ΔK+D+ΔUc−Jc+R`，没有把R加到D或隐藏成阻尼。最大已接受归一化残差 **{100*result['max_relative_R']:.7f}%**。

标量账本逐接受步保存；损伤变化逐事件以FP64记录；完整节点坐标、速度、力和损伤快照按约25µs、首次事件及协议端点保存，并非逐步全量坐标视频流。合法checkpoint带原生restart、接受态、表面历史、节点功和谱系/总账sidecar；仅当2ms试运行通过后才从同一checkpoint继续3.6ms，实际覆盖时间以上述原生终态为准。失败原生态与最后接受态分开，`LAST_LEGAL_CHECKPOINT.json`指向最近合法可恢复磁盘状态；若它早于最后接受NPZ，两者明确区分。

异常记录：

{failure}

## 7. Q4：邻域、断裂路径和方向变化

初始拓扑变化如第4节。实际断裂指标见 `CRACK_METRICS.json`；仅对本体新断键统计，排除原预缺口。用断键方向二阶张量、四阶坐标轴对齐量和断键中点空间分布评价，不能只看渲染是否“不方正”。

匹配接受时间V3/V4={crack['times_s']['V3']*1e3:.7f}/{crack['times_s']['V4']*1e3:.7f}ms；本体新断键={crack['metrics']['V3']['bulk_failed']}/{crack['metrics']['V4']['bulk_failed']}。四阶轴对齐量V3/V4={crack['metrics']['V3']['axis_alignment_fourth_moment']}/{crack['metrics']['V4']['axis_alignment_fourth_moment']}，该量只是原始键方向代理。相对各自初始本体键方向背景扣除后的V3/V4差值为{outcome['final_axis_alignment_excess_over_reference']['V3']}/{outcome['final_axis_alignment_excess_over_reference']['V4']}；必须同时看背景与断裂演化，不能把初始邻域方向差异直接解释为裂纹偏置改善。

输出的网络路径比是最大断键网络上选定端点的最短图路径/直线距离，端点规则及坐标均保存；它不是识别出的连续裂纹曲折度。由于参考位置、邻域和积分权重共同变化，不能把差异归因于单纯的裂纹外观。

## 8. Q5：20–200节点游离块是否增加？

按 **SURVIVING_ORIGINAL_MATERIAL_GRAPH** 统计；匹配接受时间的20–200节点游离分量V3/V4为 **{outcome['final_target_count']['V3']}/{outcome['final_target_count']['V4']}**；V4时程最多 **{result['max_N20_200']}** 个。匹配观察窗口的历史峰值V3/V4为{outcome['matched_window_max_target_count']['V3']}/{outcome['matched_window_max_target_count']['V4']}。对应物理质量V3/V4={outcome['target_mass_kg']['V3']*1e6:.8f}/{outcome['target_mass_kg']['V4']*1e6:.8f}mg。

`FRAGMENT_HISTORY.csv`、`FRAGMENT_SPECTRA_HISTORY.csv`、`FRAGMENT_SIZE_SPECTRA.json`与`FRAGMENT_MORPHOLOGY.csv`同时给出原始节点数、参考体积、质量和以50µm立方体体积换算的等效点数。20–200仍以原始节点数定义，并不等于固定物理体积。形貌依据质量加权协方差和主轴比例作描述性COHESIVE_BLOCK/THIN_SHEET/FILAMENT/SMALL_DEBRIS分类；不改变真实图成员，也不把空间相邻点并成块。

## 9. Q6：卸载后是否更稳定？

V4是否到达2.8ms撤载：**{survival['V4']['drive_off_reached']}**。撤载时目标成员集合数、末端完整存活数分别为V3 **{survival['V3']['unload_cohort_count']}/{survival['V3']['unload_cohort_survived']}**，V4 **{survival['V4']['unload_cohort_count']}/{survival['V4']['unload_cohort_survived']}**。未到达撤载时用null表示不可评价，不填0冒充已观察失败。

全部成员集合出生/真实分裂时间和父子关系见 `FRAGMENT_LINEAGE.json`、`FRAGMENT_SURVIVAL.json`。右删失表示观察结束时尚未分裂，不代表永久稳定；一次轨迹也不能给出统计意义上的稳定性提升。

## 10. Q7与Q8：物理结论及后续价值

{scientific}

当前不规则实现已通过相对弹性、完整变分梯度、非等体积原生质量/接触、状态与表面映射检查，具备作为探索性离散方案继续研究的基础。是否优于规则方案须综合本轮实际碎块谱、撤载存活和成本；若本轮未走完整协议，先处理记录的数值停止原因，不从短时没有目标块推断失败，更不靠改强度、增加种子或载荷扫描自动补结果。

视频只使用真实原生接受态、位移×1、同一世界坐标范围与固定Z向上相机。用球体表现节点积分体积，不是连续碎块外壳或真实球形纤维；完整存活原始图决定分量，局部目标块显示全部内部存活键。两组按真实保存时间匹配，逐面板标明实际时间，不插值伪造中间断裂态。

## 11. 交付与复现

十张主图、两个1080p MP4见 `OPEN_RESULTS.html`。源码在 `code/`；参考/阈值/原始ID、体积和旧新键映射、全部接受快照及损伤记录在本轮verification归档。源码、必要证据与轻量日志进入独立Git提交；原生大库、对象、密集状态和重复事件载荷保留服务器。`PERFORMANCE_BASELINE_V3.json`在交付目录为保留完整计时的紧凑版，记录原始完整文件SHA和服务器地址，没有删改原始档案。

详细执行顺序、环境和复制边界见 `REPRODUCE.md`。本轮没有自动推送GitHub，没有新CFD、声场、微泡动力学或AI训练。
'''
 (O/'REVIEW_CLOT_IRREGULAR_LPS_V4.md').write_text(text)
 reproduce=f'''# V4 reproduction

Use the isolated server environment `/workspace/y_clot/env.sh`, single MPI rank and one BLAS/OpenMP thread. Native library SHA: `{sha(LIB)}`. The original V3 baseline archive and V2 native build objects are required: this compact Git snapshot is not a standalone binary distribution.

1. Preserve existing accepted archives. Run `python code/restore_export.py` in a compact Git checkout to restore and verify the losslessly compressed native reference, body data and correspondence tables. Reproduction must use a fresh project copy or a new RUN configured in `code/common.py`; output guards deliberately refuse replacing a frozen reference or trajectory.
2. Verify `SOURCE_FREEZE.json`, `V4_IMPLEMENTATION_FREEZE.json`, environment package records and original V3 sources.
3. Run `run_p0.sh`, `run_p1.sh`, then `performance_summary.py` to reproduce the paired eight-step profiling, before generating the irregular reference. Do not repeat the historical full V3 dynamics.
4. Run `prepare_reference.py`, `build_native.py`, and `build_surface.py`. The native builder verifies and reuses the frozen V2 object hashes, recompiles changed classes/contact/dispatch, and writes isolated V4 outputs. Full commands and source/object hashes are in `BUILD_AUDIT.json`.
5. Prepend the new output `build` directory to `LD_LIBRARY_PATH` **before** importing LAMMPS. Run `check_elastic.py`, `check_native.py`, then `admit.py`. The tests use separate small validation coupons/fixtures; they are not additional formal cases.
6. `run_formal.sh` runs `run_case.py pilot` to 2.0ms and `run_case.py continue` from its legal checkpoint. There is one physical case. If a numerical gate fails, stop: do not weaken gates, refit materials or add expensive cases.
7. After terminal status, run `terminal_audit.py`, `analyze_result.py`, `render_videos.py`, `write_report.py`, `validate_delivery.py`. `reference_figures.py` and `compact_performance.py` prepare the reference-only outputs and small performance payload.

Current service: `supervisorctl status clot_v4_formal_20261019T221000Z`. Log: `{O}/formal.supervisor.log`. Postprocessing is a separate managed service, never another physical trajectory.

The initial generation uses exactly seed20261019, bounded isotropic ball proposals, fixed exterior/notch/anchor positions, a documented geometric strong-volume projection, actual mapped-cell volumes, and KD-tree original neighborhoods. Never substitute nominal equal volumes in a V4 replay.

Full original FP64 states, native restarts, build objects and dense damage streams remain at `{O}` on the server. The Git delivery includes code, manifests, small scientific tables, plots and selected evidence; large source payloads are identified explicitly in export manifests. Videos and dense data must be copied from the result archive if excluded by the local Git size policy.
'''
 (O/'REPRODUCE.md').write_text(reproduce)
 names=['SOURCE_FREEZE.json','PERFORMANCE_BOTTLENECKS.md','PERFORMANCE_OPTIMIZATION.json','IRREGULAR_REFERENCE_GEOMETRY.json','LPS_DISCRETIZATION_CHECK.json','NATIVE_V4_RESULT.json','FRAGMENT_HISTORY.csv','FRAGMENT_SPECTRA_HISTORY.csv','FRAGMENT_MORPHOLOGY.csv','FRAGMENT_LINEAGE.json','ENERGY_LEDGER.csv','COMPUTE_COST.json','REVIEW_CLOT_IRREGULAR_LPS_V4.md','REPRODUCE.md']
 for name in names:shutil.copy2(O/name,D/name)
 shutil.copy2(O/'PERFORMANCE_BASELINE_V3_COMPACT.json',D/'PERFORMANCE_BASELINE_V3.json')
 for name in ['REFERENCE_BOND_TOPOLOGY.csv','NODAL_VOLUME_AUDIT.csv']:
  dest=D/name
  if not dest.exists():dest.symlink_to(os.path.relpath(O/name,D))
 for name,dest in [('figures',F/'figures'),('videos',F/'videos'),('states',O/'cases/H_IRREGULAR_LPS_V4/keys'),('logs',O)]:
  path=D/name
  if not path.exists():path.symlink_to(os.path.relpath(dest,D))
 cards='\n'.join(f'<figure><a href="figures/{p.name}"><img src="figures/{p.name}" loading="lazy"></a><figcaption>{html.escape(p.stem)}</figcaption></figure>' for p in sorted((F/'figures').glob('FIG_*.png')) if p.name not in ['FIG_01_NATIVE_TIME_BREAKDOWN.png','FIG_02_BEFORE_AFTER_TIMING.png','FIG_03_REGULAR_VS_IRREGULAR_LATTICE.png'])
 vids='\n'.join(f'<figure><video controls preload="metadata" src="videos/{p.name}"></video><figcaption>{html.escape(p.stem)}</figcaption></figure>' for p in sorted((F/'videos').glob('*.mp4')))
 page=f'''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>V4 irregular LPS</title><style>body{{font:17px system-ui;max-width:1500px;margin:auto;padding:24px;background:#fafafa;color:#26313b}}h1{{color:#991e3e}}figure{{margin:20px 0;background:white;padding:16px}}img,video{{width:100%;height:auto}}figcaption{{padding:10px;color:#56636c}}.status{{padding:20px;background:#f3e7eb}}a{{color:#8d2040}}</style><h1>V4 不规则LPS：实际计算与离散研究</h1><div class="status">{html.escape(result['status'])} · 接受到 {result['accepted_time_s']*1e3:.8f} ms / 3.6 ms<br>EXPLORATORY_UNCALIBRATED · MECHANICAL_SHEAR_SURROGATE</div><p><a href="{os.path.relpath(D/'REVIEW_CLOT_IRREGULAR_LPS_V4.md',F)}">开发与科研报告</a> · <a href="{os.path.relpath(D/'REPRODUCE.md',F)}">复现说明</a></p>{cards}{vids}</html>'''
 (F/'OPEN_RESULTS.html').write_text(page)
 (D/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url='+os.path.relpath(F/'OPEN_RESULTS.html',D)+'"><a href="'+os.path.relpath(F/'OPEN_RESULTS.html',D)+'">打开V4结果</a>')
 dump(O/'FINAL_STATUS.json',dict(status=outcome['status'],native_status=result['status'],accepted_time_s=result['accepted_time_s'],full_protocol_completed=complete,formal_cases=1,model_label=HET_LABEL,load_label=LABEL,report='REVIEW_CLOT_IRREGULAR_LPS_V4.md',delivery_validation='PENDING'));shutil.copy2(O/'FINAL_STATUS.json',D/'FINAL_STATUS.json');print('V4_REPORT_WRITTEN',flush=True)
if __name__=='__main__':main()
