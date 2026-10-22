from pathlib import Path
import json,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import poisson
root=Path(__file__).resolve().parent
c=json.loads((root/'FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json').read_text())
v=root/'validation/server'
def read(k):return json.loads((v/(k+'.json')).read_text())
a=read('A_B');ce=read('C_D_E_H');f=read('F');g=read('G');mark=read('F_mark_batches');boundary=read('E_boundary');real=json.loads((root/'validation/real_source_only/SUMMARY.json').read_text())
figdir=root/'figures';figdir.mkdir(exist_ok=True)
fig,axs=plt.subplots(2,1,figsize=(12,6),sharex=True)
axs[0].step([0,c['T_infusion_s'],c['T_infusion_s'],c['T_end_s']],[1,1,0,0],where='post',linewidth=2.5);axs[0].set_yticks([0,1],['OFF','ON']);axs[0].set_ylabel('Infusion')
for birth in ce['result']['source_proposals']:axs[1].plot([birth['birth_time_s'],c['T_end_s']],[birth['source_event_id']]*2,linewidth=2,marker='o',markevery=[0]);axs[1].scatter([c['T_end_s']],[birth['source_event_id']],marker='x',s=80,zorder=4,color='red')
for ax in axs:
 ax.axvline(c['T_infusion_s'],ls='--',color='#ba652c');ax.axvline(c['T_end_s'],ls='--',color='#b12c45');ax.set_xlim(0,c['T_end_s']+.18);ax.grid(alpha=.15)
axs[0].set_title('Fixed infusion + 10 s observation; endpoints independent of outcomes')
axs[0].text(c['T_infusion_s']+.12,.82,'INFUSION STOP\n1.240815660 s',color='#ba652c')
axs[0].text(c['T_end_s']-.12,.82,'OBSERVATION END\n11.240815660 s',ha='right',color='#b12c45')
axs[1].set_ylabel('Synthetic particle ID');axs[1].set_xlabel('Global physical time (s)');axs[1].set_title('Actual small CPU/GPU validation: two particles right-censored at fixed end')
fig.text(.02,.015,'Validation only: constant synthetic flow, dt = 0.37 s, source rate = 2/s. Production contract keeps dt = 0.0005 s and rate = 1200/s.',fontsize=9)
fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(figdir/'fixed_window_timeline.png',dpi=300);fig.savefig(figdir/'fixed_window_timeline.pdf');plt.close(fig)
x=np.loadtxt(v/'G_source_only_counts.csv',delimiter=',',skiprows=1)[:,1];fig,ax=plt.subplots(figsize=(10,5));ax.hist(x,bins=np.arange(x.min()-.5,x.max()+5,5),density=True,alpha=.6,label='1024 source-only realizations');k=np.arange(int(x.min())-5,int(x.max())+6);ax.plot(k,poisson.pmf(k,g['expected']),label='Poisson theoretical PMF');ax.set_xlabel('Proposed source events in fixed infusion window');ax.set_ylabel('Probability per event count');ax.legend();ax.set_title(f"Expected {g['expected']:.3f}; sample mean {g['mean']:.3f}; sample variance {g['sample_variance']:.3f}");fig.tight_layout();fig.savefig(figdir/'source_only_poisson.png',dpi=300);fig.savefig(figdir/'source_only_poisson.pdf');plt.close(fig)
text=f'''# 固定灌注与停灌后10秒观察：开发与验证报告

正式协议为 **FIXED_WINDOW_POPULATION_EXPERIMENT_V1**。以后所有由当前统一生产入口启动的 realization，都使用合同规定的固定输入时间窗，不再根据出口数量决定停灌时刻。此次仅小规模验证、真实源事件重放和零源 CLI 检查，没有启动正式大批量粒子实验或任何 CFD。

## 1. 已冻结的参数

|参数|正式值|
|---|---:|
|T_infusion|{c['T_infusion_s']:.9f} s|
|T_observation，明确指停灌之后|{c['T_observation_s']:.1f} s|
|T_end|{c['T_end_s']:.9f} s|
|浓度 C|{c['concentration_m3']:.6g} m⁻³|
|入口 Q|{c['Q_in_m3_s']:.16g} m³/s = 18 mL/min|
|λ=CQ|{c['lambda_source_s_inv']:.16g} s⁻¹|
|E[N_source]=λT_infusion|{c['lambda_source_s_inv']*c['T_infusion_s']:.9f}|
|生产名义 dt|{c['nominal_dt_s']:.7f} s = 0.5 ms|
|默认 master_seed|{c['master_seed']}|
|默认 realization_id|{c['realization_id']}|
|flow SHA256|`{c['flow_sha256']}`|
|SonoVue source distribution SHA256|`{c['source_distribution_sha256']}`|

本轮检查发现旧正式配置只有 coverage stop，没有显式固定灌注时长。因此采用用户指定的前次实际时长 **1.240815660 s（9位小数）** 作为新的输入参数。旧实际停灌为 1.2408156603855842 s；不把那次出口覆盖过程再运行一遍去决定新 T_infusion。新合同的终点按十进制定义求和后转为 float64；积分器必须到达这一约定值，不多走一个 dt。

合同 `FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json` 显式包含 stopping_policy、outcome_dependent_stop=false、right_censor_policy=ACTIVE_AT_T_END、源分布哈希及流场哈希。允许通过统一入口改变 seed/realization/C/dt/backend；修改后的实际参数写入每次新 run 的 resolved contract，时间不由 shell 参数覆盖。

## 2. 实际代码更改

源目录：`/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/src/particle_3d/population_v1/`。

- `run.py`：原先实际出口覆盖触发停止和“池为空就提前结束”均已从当前生产入口移除；兼容模块只转发到新入口，旧 --minimum-exits、--post-infusion、--fixed-infusion 参数被拒绝。
- `run_fixed_window_population.py:main`：唯一正式入口，强制读合同、检查与 runtime flow/source 身份一致，只允许新目录，保存实际合同、源哈希、源码逐文件哈希、设备信息和 previous/current protocol。
- `fixed_window.py:run_fixed_window`：下一个时刻只取 nominal tick、下一个合法源事件、T_infusion、T_end 的最小值；下游终止/出口计数不参与该选择。池为空也把全球时钟推进到固定端点。
- `fixed_window.py:censor_remaining`：恰好 T_end 后保存仍活动的微泡原始最后位置/速度/单元/壁面间隙和现时几何的 pair 状态，记录 RIGHT_CENSORED_OBSERVATION_END，再从活动池移除，不继续积分或接到出口。
- `fixed_window.py:FixedWindowOutput`：另存 `source_proposals.jsonl.gz`，保留原有带准入结果的 `source.jsonl.gz`；所有终态增加互斥 terminal_category。
- `scripts/launch_production.py`：不再带旧阈值或12秒时长，改为指定正式合同，默认目录改成 `realization_0000_fixedwindow`。

核心时序：

```python
birth = row['birth_time_s'] if row is not None and row['birth_time_s'] < tin else float('inf')
target = min(tick if engine.pool.n else float('inf'), birth, tin if source_on else end, end)
engine.advance(float(target))
if source_on and target == tin:
    source.stop(tin)
# 仅 birth < tin 才进入 admission；最终 engine.time == end 后执行 censor_remaining。
```

出生严格属于 **[0,T_infusion)**。预先抽到的未来 mark 不是已实现的事件；恰好发生于 T_infusion 的 mark 被排除，不算 rejection。源事件不取整到网格时刻。

**源提议与准入结果不同：** event-addressed RNG 使 arrival/直径/位置/朝向等源提议不依赖下游轨迹；但旧有 finite-size/pair admission 仍可能因为入口附近的活动微泡不同而接受/拒绝不同。要求“源过程不受 downstream 控制”不等于要求“改变 downstream 后准入决定也完全相同”。测试分别保存这两项，未修改准入或 MB–MB 模型。

终态类别为 EXITED_O1/O2/O3、MODEL_SUPPORTED_STATIONARY_TERMINATION、NUMERICAL_FAILURE、RIGHT_CENSORED_OBSERVATION_END；底层原状态名作为记录保留，writer 增加规范化 category。恒等式为 accepted = exits + supported stationary + numerical failure + right-censored。右删失不计作错误、出口或完成轨迹，不剔除它后把出口分数冒充总体概率。

最后速度是前一个已接受积分区间的 held velocity，不在终点偷偷重算/外推新速度。pair 状态使用最后位置重新查询精确球几何和原 V1 阈值，避免误用上一段开始时刻的 pair 集合。空列表表示此刻无活动近场 pair，并非漏填字段。

## 3. 实际完成的验证

|要求|实际执行与结果|证据|
|---|---|---|
|A 提前覆盖不改变停灌|调度 fixture 三出口各20次在 {a['old_20_each_outlet_coverage_time_s']:.9f} s 达成，继续提议至固定 T_infusion；实际提议 {a['proposal_count']}|validation/server/A_B.json|
|B 无结果依赖|“三个出口都很快覆盖”与“只有一个出口有流出”使用同源，同一半开时间窗内源提议逐字节一致|A_B.json；此项是调度 fixture，不是新流场求解|
|C 固定终点删失|实际 CPU 及 CUDA 合成常速场各有2泡保留到终点，均在 11.240815660 s 右删失，记录完整末态|C_D_E_H.json / F.json|
|D 精确切分|小验证 dt=0.37 s，不能整除停灌/终点，仍精确达到两个规定时刻，没有多走一步|C_D_E_H.json；生产仍 dt=0.0005 s|
|E 无停灌后出生|已实现的所有 birth<T_infusion；专门把首事件设为恰好停灌，事件数=0，空池仍走到固定终点|E_boundary.json|
|F CPU/GPU与batch可复现|实际 CPU/CUDA 2泡验证的源提议与准入 ledger 位相同；输出 chunk 1/1000一致，位置最大差 {f['position_error_m']:.6g} m；另 source marks 的 batch 1/17/256、批内逆序调度恢复 event-ID 次序后逐字节相同|F.json / F_mark_batches.json|
|G 固定窗Poisson统计|1024个 **arrival-only** realization，无粒子积分，mean={g['mean']:.6f}，variance={g['sample_variance']:.6f}，理论均值/方差={g['expected']:.6f}|G.json / G_source_only_counts.csv|
|H 互斥终态及计数|actual=accepted+rejected；accepted=全部互斥终态。测试中2 accepted=2 censored；另检验6类名称互斥和旧CLI参数拒绝|C_D_E_H.json / H_categories_CLI_guard.json|
|真实源定义重放|正式 BraVa flow/source/seed 下提议 {real['actual_source_count']} 个，最后出生 {real['last_birth_s']:.12f} s；与历史同种子源 mark 一致。未运行准入或轨迹|validation/real_source_only/|
|真实环境空源CLI|真实环境、CUDA后端、验证用 C=0，没有活动泡，仍精确经过停灌和观察终点；完整合同/manifest/时间线已保存|validation/empty_real_environment_fixedwindow_current/|

统计检验的判据在测试代码中设为6个标准误：均值 z={g['mean_z']:.6f}，方差 z={g['variance_z']:.6f}。样本方差比理论低约 {(1-g['sample_variance']/g['expected'])*100:.2f}%，本样本通过的是这个明确的宽判据，不能称方差已精确等于理论值。arrival-only 用真实 ARRIVAL role 与 poisson_interval，不包含直径/位置准入或1024次动力学；完整 mark 的一致性由 F 和真实源重放另行检查。

CPU/GPU 一致性只说明同一离散场与模型的执行一致，**不证明 CFD 物理正确**。metadata 的 no_numerical_failure_reported 仅表示该次运行没有记录数值失败；cfd_physical_accuracy_validated_by_protocol=false 明确排除物理正确性声明。

旧 production 的1565个接受泡、141个旧删失及2次数值失败没有迁移成新实验结果；没有回填新停止规则或缩短旧动画。

## 4. 图与保存的验证 run

![固定时间线](figures/fixed_window_timeline.png)

![源计数统计](figures/source_only_poisson.png)

`validation/small_synthetic_runs/run0/1/2` 是实际测试的2泡 CPU、CUDA、CUDA不同输出chunk副本；保存 IDENTITY、MANIFEST、提议/准入/终态 ledger 和接受时间状态。它们是明确标识的合成验证，不能作为 BraVa 正式生产结果。真实环境 CLI 是零浓度验证，亦不能替代含微泡正式运行。

## 5. 后续统一运行方式（本轮未启动）

服务器：

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
export PYTHONPATH="$PWD/source/particle_3d/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env/bin/python -m particle_3d.population_v1.run_fixed_window_population \\
  --environment data/environment.json \\
  --contract fixed_window_protocol/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json \\
  --output runs/realization_0001_fixedwindow \\
  --backend cuda --realization-id 1 --git-commit ACTUAL_CODE_COMMIT
```

每次选择新的目录及 realization_id；时间来自同一正式合同。代码不自动提交大批量 realization。若要在本地CPU运行，使用相同模块、`data/environment.local.json` 和 `--backend cpu`，Python 环境需要项目既有依赖。

验证命令：

```bash
FIXED_VALIDATION_DIR="$PWD/fixed_window_protocol/validation/server" \\
  env/bin/python -m pytest -q source/particle_3d/tests/population_v1/test_fixed_window.py
env/bin/python fixed_window_protocol/validate_real_source.py
```

本地实际代码与服务器源码同步，具体版本及新源码 SHA 见交付清单。旧入口的只读存档在 `migration/legacy_coverage_run.py.txt` 和 `migration/legacy_launcher.py.txt`；它们不是当前可运行生产入口。历史重现应使用历史记录的完整代码 commit，而不是把旧阈值重新加入当前协议。
'''
(root/'FIXED_WINDOW_PROTOCOL_REPORT_ZH.md').write_text(text)
(root/'README.md').write_text('''# 固定时间 population 实验\n\n主报告：[FIXED_WINDOW_PROTOCOL_REPORT_ZH.md](FIXED_WINDOW_PROTOCOL_REPORT_ZH.md)。\n\n正式合同：[FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json](FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json)。\n\nT_infusion=1.240815660 s；停灌后观察10.0 s；T_end=11.240815660 s；出生区间[0,T_infusion)。停止规则不读取出口覆盖结果。0.5 ms是正式名义步长；少量合成测试使用0.37 s以验证非整除时间窗。\n\n唯一正式入口：`python -m particle_3d.population_v1.run_fixed_window_population`，参数及服务器复现命令见主报告。当前兼容入口`population_v1.run`同样强制固定合同。没有运行正式大批量实验，没有CFD计算。\n\n源提议ledger与准入ledger分开。`RIGHT_CENSORED_OBSERVATION_END`不是失败或出口；accepted births及全部终态始终保留。\n\n迁移来源：旧production无显式固定灌注时长，因此按用户要求固化前次1.240815660 s。旧history的精确停止时间1.2408156603855842和停灌后12s观察保留原样，不能改称新10s协议。旧结果hash不变性见validation/history_immutability.json。\n\nCPU/GPU只做已有模型执行一致性检查；不替代近壁CFD审计或物理精度验证。\n''')
(root/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><title>固定实验协议</title><style>body{font:17px system-ui;max-width:1200px;margin:30px auto}img{width:100%}</style><h1>固定灌注 + 停灌后10秒观察</h1><p><a href="FIXED_WINDOW_PROTOCOL_REPORT_ZH.md">中文报告</a> · <a href="FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json">正式合同</a> · <a href="../cfd_wall_audit/OPEN_RESULTS.html">近壁CFD审计</a></p><p>未启动正式大批量实验。小规模合成CPU/CUDA回归和真实源重放证据见报告。</p><img src="figures/fixed_window_timeline.png"><img src="figures/source_only_poisson.png">')
print(root/'FIXED_WINDOW_PROTOCOL_REPORT_ZH.md')
