from pathlib import Path
import json,html,hashlib
C=json.loads(Path(__file__).with_name('config.json').read_text());L=Path(C['local']);O=L/'verification'/C['task']/C['run'];F=L/'results'/C['task']/C['run'];load=lambda n:json.loads((O/n).read_text());q=load('SCIENTIFIC_SUMMARY.json');p=load('PATCH_PROVENANCE.json');d=load('DYNAMICS_AUDIT.json');perf=load('PERFORMANCE_AND_OUTPUT_SIZE.json');doc=load('DOCUMENT_UPDATE_AUDIT.json');status=load('FINAL_STATUS.json');mi=load('MPI_SORT_AND_MIGRATION_EVENTS.json');a=load('NATIVE_FORCE_AND_TOPOLOGY_AUDIT.json');nf=q['new_final'];of=q['old_final'];ns=a['key_state_checks'];rn='REVIEW_LPS_FINAL_NATIVE_BUGFIX_FULL_V3_REGRESSION_V1.md'
sha_rows='\n'.join(f"| {k} | `{v['final_sha256']}` |" for k,v in p['final_source_sha256'].items() if isinstance(v,dict))
text=f'''# LPS_FINAL_NATIVE_BUGFIX_FULL_V3_REGRESSION_V1 科研审查

Run ID：{C['run']}。最终双修复诊断库完成唯一一次原C1的4660步、2.33 ms回归，两MPI退出码均为0。原始输入、材料、几何、载荷、接触、排序及背景流阻力未改。**FULL_V3_NATIVE_REGRESSION_PASS仅指本次原生诊断回归；生产环境没有切换，材料断裂物理尚未得到实验认证。**

## 结果对照

| 指标 | 历史接受C1 | 最终双修复 |
|---|---:|---:|
| 首次bulk断键 | 241步 / 0.1205 ms | 241步 / 0.1205 ms |
| 首次interface断键 | 621步 / 0.3105 ms | 621步 / 0.3105 ms |
| 首次真实游离及阻力开启 | 593步 / 0.2965 ms | 611步 / 0.3055 ms |
| 累计bulk / interface断键 | 11771 / 718 | 10808 / 643 |
| 末态附着材料点 / 质量比例 | 1710 / 88.739% | 1728 / 89.673% |
| 末态游离节点 / 组件 | 217 / 208 | 199 / 190 |
| 单点 / 双点 / 三点 / 四点游离组件 | 201 / 6 / 0 / 1 | 184 / 4 / 1 / 1 |
| ≥3点真实游离组件 | 1 | 2 |
| 1–2点碎屑占游离质量 | 98.157% | 96.482% |
| 末态最大真实游离块 | 4点，98.4745 µm | 4点，98.4745 µm |
| 全过程最大游离组件 | 13点 | 20点，168.389 µm，1313步 |
| 全程≥200 / ≥300 µm游离组件 | 0 / 0 | 0 / 0 |

末态最大游离块质量5.0×10⁻¹⁰ kg；新游离总质量2.4875×10⁻⁸ kg，其中1–2点未解析碎屑2.4×10⁻⁸ kg。全过程满足附着材料点＋游离材料点＝1927，组件成员和＝游离节点数。壁锚点1561个不计入血栓质量。旧`largest_free_fraction`和AI标签均保留；本报告只使用ATTACHMENT_AWARE_V2解释“游离”。

## 37项问题的逐项回答

### 原生修复（1–10）

1. **实际库**：是最终双修复库。两个实际MPI Python进程的`/proc/self/maps`均指向上轮`build_owner_fix/liblammps.so.0`，SHA一致；不是生产或m-only中间库。Python3.13.11、LAMMPS20260930、私有OpenMPI4.1.6；Python包装模块来自既有锁定源码，实际共享库由本次进程专属LD_LIBRARY_PATH选定。`ldd`和实际libmpi路径/SHA已记录。
2. **完整身份**：见下表及`SOURCE_AND_LIBRARY_IDENTITY.json`、`PATCH_PROVENANCE.json`。本轮没有新增物理C++补丁、没有重编译主库，只编译小型只读观察器。
3. **步数**：完整4660步，终态物理时间0.00233 s；234个原生接受帧和完整精度终态一致。
4. **MPI退出**：两个rank的包装进程实际退出码均0，MPI launcher退出0；不是只凭输出文件判成功。
5. **m身份**：包括初始化在内4661个状态，对全部owned和ghost逐ID验证固定完整参考邻域。最大相对误差1.12031267×10⁻¹⁵，最大绝对误差1.57772181×10⁻³⁰ m⁴，门槛2×10⁻¹²。分母是该ID真实正参考m，没有隐蔽下限。断后m不重定义。
6. **实际排序**：初始化非恒等排列为rank0的1938点、rank1的1544点；运行中真实排序在1018、2004、3200、4023步，两rank都有非恒等排列。用实际nextsort变化及原生permute记录核验，不能用计划1000/2000/3000/4000步替代。
7. **ownership迁移**：已触发，第一次逐ID证实为455步，全程34次原子跨rank转移，涉及{len(mi['ownership_unique_atom_ids'])}个不同ID。迁移、同rank顺序变化、排序及ghost列表重建分别记录。
8. **跨rank断键**：1039条，首条315步。事件CSV保存真实端点ID、owner及在两个rank中的owned/ghost下标。
9. **访问及拓扑安全**：冻结源码包含两处0≤index<nlocal检查；每次拓扑变化核对全部owned端点的ID、类型、唯一有向键及逐ID度数，global端点数＝2×nbonds。11451条原生断键事件与拓扑删除逐事件一致，没有重复、漏端点或rank计数冲突。未启用ASan，不能把“本次未发现非法访问”扩大为形式化内存安全证明。
10. **新原生错误**：本轮没有发现；主力学源码没有再改。首次静态检查的报告代码误把标量能量返回值按数组读取，修正后重复同一微型静态检查；另修正一次离线JSON的NumPy整数序列化。两者均保留记录，均非LAMMPS缺陷或正式V3失败。

### 真实断裂（11–20）

11. **首次bulk**：241步＝0.1205 ms。
12. **首次interface**：621步＝0.3105 ms。
13. **首次游离**：611步＝0.3055 ms，并在同一步开启阻力；没有强制沿用旧593步。
14. **末态真实游离节点**：199个，10.3269%初始材料质量。
15. **末态游离组件**：190个。
16. **最大真实游离组件**：末态4点；全程最大20点，首次达到20点为1313步（0.6565 ms）。
17. **等效球径**：末态98.4745 µm；全程已证实最大168.3890 µm。
18. **≥200 µm**：全程没有；≥300 µm也没有，首次阈值时间为“不发生”，不是0 ms。
19. **碎屑比例**：192个游离点属于1–2点组件，占游离质量96.4824%；占初始总血栓质量9.96367%。不能称其为实验认证的微栓。
20. **速度与动能**：每步原生final_integrate之后读取的全程vmax为3.9133943818 m/s，759步（0.3795 ms）、ID364。最大步间位移1.95708138 µm，760步、ID364，约0.03914Δx。全体材料点动能最大4.37534617×10⁻⁸ J（911步）；游离材料点动能最大2.06164087×10⁻⁸ J（954步）。旧3.9199684 m/s仅为2104个已保存力时刻的完整步速度重建最大值；旧234个轨迹帧最大3.9130832 m/s。三种采样口径明确分开。

### 历史差异（21–25）

21. **首次可核实力学分叉**：第455步（0.2275 ms）。历史保存到454步的可对齐记录相同；455步坐标、半步速度、拓扑仍逐位相同，pre-flow节点力最大差11.6499638 µN，最差ID386。第455步接受完整速度的最大差已为0.0232999 m/s（旧完整速度由保存的半步速度与力重建，新值为直接原生读取）；456步首次保存位置差23.2999 nm、半步速度差0.0465999 m/s；462步首次原生断键事件分叉。历史缺失时刻不插值、不补造事件。
22. **归因**：455步真实迁移且未排序。依据旧裸数组源码及新运行当时实际ID顺序，重建旧局部m槽位，推得ID386和3480的m错配；该重建预测旧力与历史保存旧力相对差4.90207×10⁻¹⁶。这是强定量支持的源码重建推断，**旧m和旧theta并未直接保存**。跨rank guard从315步起已被覆盖，但此前可比轨迹相同；没有证据把455步力差归给guard。之后排序、断裂接触及阻力反馈共同传播差异。
23. **断裂演化**：累计bulk少963条、interface少75条；总删除11451条、末态存活146186条。首次bulk/interface时间相同，首次游离推迟9 µs。每步速率、断裂参考中点位置、动能增量分别保存在`FRACTURE_RATE_AND_KINETIC_CHANGES.csv`和`FRACTURE_LOCATION_AND_BURST_AUDIT.json`。最大单步断键38条，524步，约占初始本体键的0.046%；无NaN/Inf或原生崩溃。有限数值不等于已证明时间收敛。
24. **旧细碎屑是否主要由bug造成**：证据不足。原生bug确实改变力和后续轨迹，但修复后末态仍以1–2点碎屑为主，不能把细碎化主因全部归给它。
25. **瞬时更大组件的含义**：13→20点是数值轨迹改变；仍未达到200 µm，也没有新增材料实验或断裂能依据，不能据此宣称物理模型改善或真实大块断裂已验证。

### 科学限制（26–30）

26. **theta限制**：仍为EXPLICIT_SCHEME_LIMITATION，未改算法顺序。{ns}个关键完整精度状态中，原生力与“断前theta＋断后存活键＋实际接触”预测最大相对差{a['max_force_relative_error_key_states']:.6g}；这些关键状态的滞后力与同坐标断后梯度力最大差1.83192 µN，**不是全4661步该误差的全局上界**。每步动能变化另存，不把拓扑释放或接触冲量全部算作theta误差。
27. **实际覆盖机制**：断后粒子接触从436步开始，最多36对；独立读取真实邻居表和special mask，经原Pair::single重建接触力，并与原生总力核对。所有199个实际脱离ID都在脱离当步开启背景阻力；附着点误加游离阻力、游离点误加附着阻尼、无效CFD游离点记录均为0。未触发管壁罚接触；V3数组增长、restart读回与5 s长期交接未执行或未触发，标NOT_EXERCISED。终态完整精度x/v和原交接源逻辑保留。
28. **Gc**：仍未标定。独立LPS势能梯度核查通过；静态native PE×V仍比相容势能低约2.69276%，不能用于Gc分子。PVCE的G_drive继续只表示试样驱动力。数值阻尼/阻力功的时间积分不是生物耗散标定，也不是完整能量守恒证明。
29. **是否改s00或本构**：没有足够新证据。本轮保持G=940 Pa、K=33333.333333333336 Pa、bulk s00=0.15、interface s00=0.140624859375、alpha=0、dx=50 µm、horizon=150.5 µm；没有为增大片径调参。
30. **生产保护**：生产lib SHA、env.sh及2703个历史/生产保护文件共3.979 GB全部匹配；原Y背景流/壁面资产也逐SHA匹配既有注册清单。没有git push、覆盖历史、永久环境切换或终止其他任务。

### 交付与下一步（31–37）

31. **实际新增计算次数**：1次正式4660步V3，另有同一微型静态检查2次调用（首次数值已计算、报告整理TypeError，修正后重复）。共3次原生Python/LAMMPS实例；没有第二次完整V3、原版重跑、dt/2 Y算例或参数扫描。只读后处理不计为原生计算。
32. **资源**：LAMMPS循环187.013 s；启动到退出189.286 s。两rank CPU合计约377.96 s；进程树RSS采样峰值1,912,688,640字节（1.913 GB），rank0/1自身峰值分别1,391,972/452,816 KiB。成功离线主分析约46.8 s；图像约10 s。初步总输出约289 MB，最终精确文件规模见manifest。没有新增近GB逐节点审计；原callback的两个大文件sink为本轮/dev/null，其原audit_rows计数表示尝试写入数，不能误称实际保存数。完整原始轨迹/重启在服务器，WSL为含精度不变NPZ和输入压缩文件的紧凑交付。
33. **可视化**：4组实际数据图：m/owner/排序迁移身份；6个严格同一步状态的旧/新对照；断键/游离/速度曲线；终态真实游离拓扑。没有生成连续医学壳、合并组件或捏造未执行时刻。静态OPEN_RESULTS.html提供查看与数据链接。
34. **服务器原始数据**：`/workspace/y_clot/project/verification/{C['task']}/{C['run']}/`；图片在对应`results/{C['task']}/{C['run']}/`。运行子目录`C1_final`，代码在`analysis/final_v3_v1`。
35. **WSL文档**：已实际修改`/home/lzy/projects/clot_rupture.md`，读回一致，保留①–⑨标题和历史217/208数据，并列新结果。
36. **备份/diff/SHA**：本目录`document_provenance/`及`DOCUMENT_UPDATE_AUDIT.json`。修改前SHA `{doc['before_sha256']}`；修改后SHA `{doc['after_sha256']}`。实际写前再次核对SHA，未发现并发修改。
37. **下一轮唯一优先项**：在冻结最终双修复库的受控小试样上，独立审查断键时间层与时间步敏感性，量化正常拓扑释放、显式theta滞后与接触反馈的区别。本轮停止，不自动启动下一轮或部署诊断库。

## 关键SHA256

| 对象 | SHA256 |
|---|---|
| 最终诊断liblammps.so | `{p['final_lib_sha256']}` |
| 既有最终三文件补丁 | `{p['patch_sha256']}` |
| 未修改生产库 | `d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a` |
{sha_rows}

实际源码/构建目录来自上轮真实清单（`source_owner_fix`、`build_owner_fix`）；本轮没有制作空的新补丁。过程中的只读观察器有独立源码、编译命令和SHA，不混作力学修复。

## 分项状态

| 项目 | 状态 |
|---|---|
'''+''.join(f'| {k} | {v} |\n' for k,v in status.items())+'''
## 数据索引和复现

- `DESIGN_FREEZE.json`、`V3_INPUT_IDENTITY_AUDIT.json`：实际原C1输入与冻结条件；`inputs/body.data.gz`解压后SHA与原输入相同。
- `RUN_COMMAND.json`、`SOURCE_AND_LIBRARY_IDENTITY.json`：MPI命令、源码及实际加载库；各rank实际退出文件在`C1_final`。
- `WVOLUME_RUNTIME_IDENTITY_AUDIT.json`、`MPI_SORT_AND_MIGRATION_EVENTS.json`、`CROSS_RANK_BREAK_AUDIT.json`：逐ID、排序、迁移、ghost和两端拓扑证据。
- `NATIVE_FORCE_AND_TOPOLOGY_AUDIT.json`、`KEY_STATE_FORCE_ENERGY_CHECKS.csv`：静态与关键状态的力/能量一致性和显式时间层限制。
- `EARLY_FLOW_COUPLING_AUDIT.json`、`FRAGMENT_ATTACHMENT_AWARE_V2.json`、`DYNAMICS_AUDIT.json`：真实附着、游离、阻力及每步动力学。
- `V3_OLD_NEW_COMPARISON.csv`、`FIRST_MECHANICAL_DIFFERENCE_AUDIT.json`、`FIRST_MIGRATION_COUNTERFACTUAL.npz`：严格步/ID对比和标明推断口径的旧m重建。
- `C1_final/accepted_step_scalars.csv`、`break_events.csv`、`identity_events_rank*.jsonl`、`states/`：轻量每步记录与无损FP64关键状态。
- `REPRODUCE.txt`：只读后处理命令和未来原生复验的独立目录要求。当前目录有原子创建的一次运行锁，禁止重复启动。
- `SERVER_OUTPUT_SHA256.json`、`LOCAL_DELIVERY_SHA256.json`及`DELIVERY_VERIFICATION.json`：完整服务器与紧凑WSL交付清单；清单不自我散列，最终校验记录单独存放。
'''
(O/rn).write_text(text);(F/rn).write_text(text)
vr=f'../../../verification/{C["task"]}/{C["run"]}/';figs=['Fig1_native_identity','Fig2_old_new_evolution','Fig3_fracture_dynamics','Fig4_final_fragments'];titles=['原生m、排序与迁移身份','严格同一步的历史/修复后演化','断键、游离拓扑与每步速度','末态真实碎片与质量守恒'];datalinks=['FINAL_STATUS.json','V3_OLD_NEW_COMPARISON.csv','FIRST_MECHANICAL_DIFFERENCE_AUDIT.json','WVOLUME_RUNTIME_IDENTITY_AUDIT.json','CROSS_RANK_BREAK_AUDIT.json','NATIVE_FORCE_AND_TOPOLOGY_AUDIT.json','FRAGMENT_ATTACHMENT_AWARE_V2.json','DOCUMENT_UPDATE_AUDIT.json','REPRODUCE.txt']
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>最终双修复完整V3回归</title><style>body{max-width:1500px;margin:32px auto;padding:0 22px;color:#20303c;font:17px/1.7 system-ui,sans-serif;background:#f8fafc}h1{font-size:30px}h2{font-size:24px}a{color:#24658e}section{background:white;border:1px solid #d8e1e8;border-radius:8px;padding:24px;margin:24px 0}img{width:100%;height:auto}li{margin:6px 0}code{overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:15px/1.7 system-ui,sans-serif}.status{padding:18px;background:#e8f3ef;border-left:5px solid #399777}</style>'''
page+=f'<h1>最终双修复完整 V3 回归</h1><p>Run {C["run"]} · 原 C1 · 2 MPI · 4660步 · 2.33 ms</p><div class="status"><strong>完整原生回归通过。</strong>生产库保持未切换。末态199个游离点 / 190个组件；最大4点、98.47 µm；全过程最大20点、168.39 µm，≥200 µm为0。显式theta、原生PE和材料Gc限制保留。</div><p><a href="{rn}">完整37项中文审查报告</a> · <a href="{vr}V3_OLD_NEW_COMPARISON.csv">下载新旧对照CSV</a></p>'
for title,name in zip(titles,figs):page+=f'<section><h2>{title}</h2><a href="{name}.png"><img loading="lazy" src="{name}.png" alt="{title}"></a><p><a href="{name}.png">打开原分辨率PNG</a>'+ (f' · <a href="{name}.svg">SVG</a>' if (F/(name+'.svg')).exists() else '')+'</p></section>'
page+='<section><h2>审计数据</h2><ul>'+''.join(f'<li><a href="{vr}{v}">{v}</a></li>' for v in datalinks)+'</ul><details><summary>展开37项完整审查报告</summary><pre>'+html.escape(text)+'</pre></details></section></html>'
(F/'OPEN_RESULTS.html').write_text(page);(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>V3回归图与数据</title><h1>最终双修复完整V3回归</h1><p><a href="../../../results/{C["task"]}/{C["run"]}/OPEN_RESULTS.html">打开四组审查图和37项报告</a></p><p><a href="{rn}">完整报告</a></p></html>')
(O/'REPRODUCE.txt').write_text(f'''只读复现（不启动新的LAMMPS运行）：
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project
python analysis/final_v3_v1/analyze.py
python analysis/final_v3_v1/figures.py

以上只重建本轮派生审计/图片；不得将O/P指向历史接受输出目录。
原始依据：C1_final/in.lammps,body.data,request.json,LOAD_PROVIDER.json；输入SHA见V3_INPUT_IDENTITY_AUDIT.json。
WSL紧凑包使用inputs/body.data.gz，解压后核对原SHA；原生文本全轨迹和final.restart保留服务器。
已执行的精确MPI argv、cwd和环境见RUN_COMMAND.json；两个实际进程映射见C1_final/LOADED_LIBRARY_rank*.json。
最终库复用上轮build_owner_fix，禁止用中间build或生产库冒充。source env.sh之后，只有当前子进程LD_LIBRARY_PATH以最终build目录优先；不永久修改env.sh。
源码桥observer.cpp仅增加END_OF_STEP只读回调、状态读取和contact_readout（原Pair::single只读计算）；不改力或拓扑。编译命令见OBSERVER_BUILD.json。
未来如获授权进行新的原生复验，必须换新run_id、新输出目录、重新核对所有冻结SHA和资源；不要删除FORMAL_V3_LAUNCHED.json来重复当前运行。start_managed.py只管理专用Supervisor任务，autorestart=false。
原静态报告脚本的float索引错误和离线JSON整数序列化问题均在科学报告中计数；没有隐藏失败或第二次正式V3。
本轮真实WSL文档更新由任务工作目录update_document.py执行，备份及SHA在document_provenance和DOCUMENT_UPDATE_AUDIT.json。
''')
print('REPORT_AND_INDEX_WRITTEN',len(text))
