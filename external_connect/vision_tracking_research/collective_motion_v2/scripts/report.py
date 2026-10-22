"""Human-readable delivery derived from frozen numeric records."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from vcommon import ROOT,OUT,DATA,LIMITS,NAMES,P,dump,identity

def read(p):return json.loads((OUT/p).read_text())
def write(p,t):(OUT/p).write_text(t.strip()+'\n')
def nums(a,d=3):return ', '.join('null' if v is None else f'{v:.{d}f}' for v in a)
def pct(x):return 'undefined (0 accepted)' if x is None else f'{100*x:.2f}%'
def main():
 a=read('results/analysis.json');ch=a['champions'];m=a['methods'];f=read('results/model_freeze.json');suite=read('training/suite.json');ev=read('results/final_evaluation_complete.json');rt={n:read(f'results/real_replay_{n}.json') for n in ['PIV_SINGLE','direct_gray','direct_fused','geo_fused']};diag=read('results/posthoc_collapse_diagnostic.json')
 names={'PIV_SINGLE':'PIV single pass','direct_gray':'Direct6 gray','direct_fused':'Direct6 + heat','direct_heat':'Direct6 heat only','direct_fused_robust':'Direct6 + heat robust','geo_fused':'GeoMotion','direct_fused_2frame':'Direct6 two-frame','geo_no_local_supervision':'GeoMotion no local loss'}
 table='| Method | MAE denominator | Six-parameter MAE (vx, vy, div, vort, s1, s2) | Accepted / 1536 | Failed accepted |\n|---|---:|---|---:|---:|\n'
 for key in ['PIV_SINGLE','direct_gray','direct_fused','direct_heat','direct_fused_robust','geo_fused']:
  r=m[key if key=='PIV_SINGLE' else ch[key]];table+=f'| {names[key]} | {r["n_mae_denominator"]} | {nums(r["mae_all_finite_global"])} | {r["n_accepted"]} | {r["n_accepted_failed"]} |\n'
 risk='| Method | Error risk at 25% diagnostic coverage (384/1536) |\n|---|---:|\n'
 for key in ['PIV_SINGLE','direct_gray','direct_fused','direct_heat','direct_fused_robust','geo_fused']:
  k=key if key=='PIV_SINGLE' else ch[key];risk+=f'| {names[key]} | {pct(a["matched_coverage"]["0.25"][k]["risk"])} |\n'
 runtable='| Method | Pure model p50/p95 ms | Full pipeline p50/p95/max ms | Consumed /120 | Fresh six-parameter Hz |\n|---|---:|---:|---:|---:|\n'
 for n,r in rt.items():
  pure=r['pure_model_ms'];runtable+=f'| {names[n]} | {"N/A" if pure is None else nums([pure["p50"],pure["p95"]],2)} | {nums([r["full_pipeline_ms"][x] for x in ["p50","p95","max"]],2)} | {r["n_consumed"]} | {r["fresh_six_hz"]:.2f} |\n'
 traintable='| Configuration / seed | Epochs / best | Optimizer steps | Wall s | Stop reason |\n|---|---:|---:|---:|---|\n'
 for r in suite['runs']:traintable+=f'| {r["config"]["name"]} / {r["config"]["seed"]} | {r["epochs_completed"]} / {r["best_epoch"]} | {r["optimizer_steps"]} | {r["resource"]["elapsed_s"]:.2f} | {r["stop_reason"]} |\n'
 timed={'training_suite':suite['elapsed_s'],'validation_and_calibration_predictions':sum(r['validation_timing']['wall_s']+r['calibration_timing']['wall_s'] for r in f['models'].values()),'final_predictions_including_H3':sum(r['timing']['wall_s'] for r in ev['evaluations']),'PIV_CPD_accuracy_runs':sum(read(f'results/{n}_final_test.json')['wall_s'] for n in ['PIV_SINGLE','PIV_MULTIPASS','CPD']),'causal_replays_monitored':sum(r['resource']['elapsed_s'] for r in rt.values()),'server_data_generation':sum(read(f'evidence/server_data_manifests/{s}/manifest.json')['wall_s'] for s in P['source_splits']),'supplemental_demo':read('results/synthetic_demo/manifest.json')['wall_s']}
 cold={n:read(f'results/cold_process_{n}.json') for n in ['direct_gray','direct_fused','geo_fused'] if (OUT/f'results/cold_process_{n}.json').exists()}
 if cold:timed['new_process_cold_probes']=sum(r['total_since_script_start_s'] for r in cold.values())
 coldtext='; '.join(f'{names[n]}: {r["total_since_script_start_s"]:.3f}s (imports {r["imports_s"]:.3f}s, checkpoint/CUDA {r["checkpoint_and_cuda_initialization_s"]:.3f}s, first inference {r["first_inference_with_transfer_s"]:.3f}s)' for n,r in cold.items())
 batchtext='; '.join(f'{names[n]}: {1536/next(r["timing"]["wall_s"] for r in ev["evaluations"] if r["model"]==ch[n] and r["corruption"] is None):.1f} triples/s' for n in ['direct_gray','direct_fused','geo_fused'])
 gpu_step_wall=sum(r['step_s']['mean']*r['step_s']['n'] for r in suite['runs'])
 # The server manifest copies are local audit evidence; add audited video rendering separately.
 if (OUT/'videos/manifest.json').exists():timed['video_render_full_decode_local']=read('videos/manifest.json')['wall_s']
 gpu=[x.get('gpu_util_percent',np.nan) for r in suite['runs'] for x in r['resource']['samples']];cpu=[r['resource']['cpu_core_equivalents'] for r in suite['runs']]
 dump(OUT/'results/runtime_ledger.json',dict(timed_stage_wall_s=timed,sum_timed_stage_s=sum(timed.values()),training_gpu_util_sample_mean=float(np.nanmean(gpu)),training_parent_cpu_core_range=[min(cpu),max(cpu)],peak_training_allocated_mib=max(r['peak_cuda_allocated_mb'] for r in suite['runs']),scope='sum of explicitly timed stages, not CUDA kernel time or whole development elapsed time; excludes install/download, smoke warmup, some setup and review; CPU excludes child DataLoader CPU'))
 write('PIV_AI_COMPARISON.md',f'''# PIV versus AI: retain PIV

Frozen final test: 1,536 triples = 1,320 VALID_GLOBAL_MOTION + 54 NON_AFFINE_MOTION + 162 NOT_OBSERVABLE. Single source recording. MAE is over all finite estimates in class0, including rejected predictions; unavailable estimates are counted separately and as failure in overall global recovery. PIV has 1,114 finite class0 estimates, whereas every AI has 1,320. Units: px/s, px/s, 1/s, 1/s, 1/s, 1/s.

{table}
PIV accepted 1,111: 1,110 class0 (84 over tolerance) plus one wrongly accepted non-affine case. Thus joint accepted risk is 85/1111 = 7.65%; class0 accepted error risk is 84/1110 = 7.57%. This is the unmodified PIV quality gate, not a 5%-risk guarantee. AI had no feasible validation threshold with at least 20 accepted samples and <=5% joint error; all formal AI fields are rejected. **0/0 is undefined risk, not zero-error success.**

{risk}
These curves retrospectively compare fixed rankings at matched coverage; they do not authorize deployment thresholds selected from test. PIV ranks within its existing accepted quality set; AI ranks algebra-valid predicted-global candidates before the selected reject-all threshold. Test rankings demonstrate that merely relaxing AI rejection would admit many errors.

Two-pass PIV and CPD were rerun only on the first 112 preregistered test samples (four material frames ×28 cases), not all1536: multipass86 accepted/7 failed; CPD41 accepted/10 failed. These subset results are never substituted for full-test denominators. F10U historical216-window results remain unchanged and separate.

{runtable}
Identical raw ROI, actual PTS, hardware and software environment; each pipeline has its own preprocessing. AI internally resizes to128² and maps the full affine generator back to raw pixels; PIV uses raw resolution. The timing comparison includes these costs. It is a short same-recording causal replay, not a repeated-session benchmark. True camera end-to-end latency is NOT_MEASURED.
''')
 write('HEATMAP_ABLATION.md',f'''# Heatmap ablation

H0 gray, H1 gray+actual detector heat, H2 heat-only, H4 separately trained robust gray+heat each have two seeds and the same data/cap/early-stop policy. H3 applies fixed dropout/false-center/jitter/severe corruption to the validation-selected H1, H2 and H4 checkpoints without retraining or recalibration. Underlying independent test gray images and truth are identical. Full numeric tables are in `results/analysis.json`; failures remain in `results/all_failures.jsonl`.

{table}
No reliable six-parameter advantage is established by adding heatmaps. H1 slightly changes shear error but its overall recovery remains far below PIV. H2 is especially sensitive: vx MAE rises from {m[ch['direct_heat']]['mae_all_finite_global'][0]:.2f} to {m[ch['direct_heat']+'_severe']['mae_all_finite_global'][0]:.2f} px/s under severe heat corruption. H4 is less sensitive, but low sensitivity is not evidence of useful robustness when outputs can collapse toward a prior.

The predicted vx standard deviation / truth standard deviation is {diag['models'][ch['direct_fused']]['std_ratio'][0]:.4f} for H1 and {diag['models'][ch['geo_fused']]['std_ratio'][0]:.4f} for GeoMotion. This post-hoc failure diagnostic was not used to retrain or select a new model. Three-frame versus two-frame and local-loss ablations also fail the joint risk gate; no temporal or geometric reliability improvement is established.
''')
 ci='| Model | Nominal | Actual six marginal coverages (%) | Six mean full interval widths |\n|---|---:|---|---|\n'
 for n in ['direct_gray','direct_fused','geo_fused']:
  for level in ['0.9','0.95']:
   r=m[ch[n]]['intervals'][level];ci+=f'| {names[n]} | {level} | {nums(np.array(r["coverage"])*100,2)} | {nums(r["mean_width"])} |\n'
 write('UNCERTAINTY_CALIBRATION.md',f'''# Independent synthetic calibration

Weights and validation rejection thresholds are frozen before final-test access. The independent calibration pool has 1,024 triples and96 class0 source-frame ×case groups. Each group score is its maximum standardized absolute residual, separately per parameter. The quantile uses ceil((G+1)*nominal) in sorted group maxima. This is an explicit finite-sample rank construction; no unverified MAPIE API is substituted. Source frames share a recording and cases reuse material, so exchangeability and independent-experiment coverage are not established. This procedure is an empirical conservative group calibration, not a real-domain guarantee.

{ci}
Coverage is often above nominal because group maxima create broad intervals. Broad intervals cannot rescue an inaccurate point estimate or the failed risk gate. GeoMotion's nominal90% divergence coverage is89.92%, while most95% marginals are above98%. Mean widths are provided to expose informativeness, and condition-specific coverage is in `analysis.json` and Figure5. No accepted-AI interval coverage exists because the operational accepted denominator is zero.

Source-frame bootstrap has only four source frames and one recording. Its intervals are descriptive, not independent physical-experiment confidence intervals. Real microscopy uncertainty remains UNCALIBRATED; no90/95% real guarantee is asserted.
''')
 write('GENERALIZATION_AND_FAILURES.md',f'''# Generalization and failures

**SINGLE_RECORDING_GENERALIZATION_NOT_VERIFIED.** Four material regions/frame pools and unique synthetic seeds prevent exact source-patch reuse across splits, but cannot establish cross-recording generalization. F10T human detection review is reused only as historical Figure1 evidence, not as six-parameter motion truth or network training masks.

All trained models fail the validation 5% joint-risk admission. The fallback rejects all final-test estimates and all real-replay AI estimates. That also rejects all216 non-global/unobservable samples, but simultaneously discards all1320 observable-global samples: this is not successful semantic classification. GeoMotion's raw classifier predicts global for all1536 final samples. Direct6+heat recognizes118/162 not-observable examples, wrongly calls54 global samples unobservable, and misses all54 non-affine examples. PIV's quality residual detects51/54 non-affine examples, wrongly accepts one and rejects two for insufficient evidence.

The most serious errors are translation and vorticity. Low variation in model outputs, visible in Figure3, shows mean/prior collapse rather than learned reliable motion. Gray-only, heat-only, fusion and robust fusion differ, but none resolves it. Finite gradient tests, independent checkpoint loading, and exp/log coordinate tests passed; these engineering checks cannot prove that the trained features recover motion.

A plausible mechanism is loss of small subpixel-motion evidence under downsampling and spatially diffuse local matching, combined with background/material shortcuts and validation class loss dominating checkpoint choice. This is a hypothesis, not a proven causal diagnosis. GeoMotion checkpoints selected at epoch1 are explicitly disclosed; later training worsened the preregistered validation score. We did not change the score after seeing final-test results.

Stationary source dots, finite feathered patches, limited real appearance diversity, per-frame blur/noise, dropout, repeated patterns and absent depth variation limit synthetic realism. Severe heat errors particularly harm heat-only input. `all_failures.jsonl`, confusion matrices, source/case breakdowns, and posthoc collapse diagnostics preserve the negative evidence. Optional NeuFlow/SynthPix are NOT_RUN and not ranked.

Next round: use new development material and a new final holdout; validate an established local matching/flow initialization, texture-supported motion losses and coarse-to-fine motion range; require visible improvement on held-out development material before expanding training. Collect independent recording sessions and known stage motion for real six-parameter truth. Keep PIV observation-only and retain its failure/rejection checks.
''')
 write('RUNTIME_AND_RESOURCES.md',f'''# Runtime and resources

Server: AMD Ryzen9 7950X3D,16physical/32logical cores; one RTX4090 24GiB. NVIDIA550.120, PyTorch2.14.0+cu126 / TorchVision0.29.0+cu126, BF16 operations and gradients verified. Actual environment, free RAM/disk, driver query and bounded64MiB fsync/warm-read probe are in `evidence/environment.json`. Other Supervisor jobs were not stopped. Isolated env passed dependency checks; a missing setuptools dependency was copied into the new env from an existing installed wheel. Slow package download and one interrupted parallel transfer were resolved; transfer hashes were verified.

Formal suite: {len(suite['runs'])} trained runs, {sum(r['epochs_completed'] for r in suite['runs'])} epochs and {sum(r['optimizer_steps'] for r in suite['runs'])} optimizer steps; suite wall {suite['elapsed_s']:.2f}s ({suite['elapsed_s']/60:.2f}min). The individual run total is {sum(r['resource']['elapsed_s'] for r in suite['runs']):.2f}s. These are substantive repeated full-dataset passes, not smoke batches. Each run early-stopped by the declared validation rule, with best/last/periodic checkpoints saved.

{traintable}
Training GPU utilization (1s device samples) averaged {np.nanmean(gpu):.1f}%; peak allocated CUDA memory {max(r['peak_cuda_allocated_mb'] for r in suite['runs']):.1f}MiB. Parent CPU core-equivalent range {min(cpu):.2f}–{max(cpu):.2f}; DataLoader child CPU is not included in that measure. Device-wide utilization includes any concurrent device work; GPU process audit showed no other compute job at launch. No unsupported kernel-time utilization total is claimed.

{runtable}
Cold new-process probes: {coldtext}. These include Python module imports and first CUDA inference, excluding interpreter spawn/SSH/video decode. Offline batch128 cached-input throughput: {batchtext}; this excludes detection and is not streaming fresh-state Hz. CUDA-synchronized training-step wall time sums to {gpu_step_wall:.2f}s (includes host/optimizer overhead; GPU kernel time was not separately measured).

Model timing uses60 synchronized warmed repetitions; cold checkpoint load is separate. Full pipeline includes detection, heat creation if used, resize/transfer/model, geometric conversion, distribution, schema and consumer send. Processing is paced to actual recording PTS with a single latest-frame slot. PIV consumed114/120, dropped6, and formed100 valid fresh windows; AI consumed120/120, formed118 windows and had zero valid outputs due to frozen rejection. Empty counts are not presented as low-risk successes. History spans66.6ms, separate from software age; interval reference age includes half that span. No-new-frame TTL checks passed.

Explicitly timed stage wall seconds: `{json.dumps(timed)}`. Their sum is {sum(timed.values()):.2f}s ({sum(timed.values())/60:.2f}min). This ledger excludes installation/download, some setup/smoke and manual review, and must not be described as whole engineering elapsed time or pure GPU kernel time. Local rendering is listed separately from server computation.

Local independent reproduction uses PyTorch2.11.0+cu128 / TorchVision0.26.0+cu128 and RTX4060. Gray/heat images and detector centers reproduce byte-for-byte; label differences are <=2.71e-12. Checkpoint inference differs by at most0.087 frozen error-tolerance units in the tested128 validation triples (BF16/hardware/version variation), not bit-identical. Optimizer state restores and changes weights on a temporary in-memory train step; formal files remain unchanged.
''')
 answers=[
 '直接复用F10U检测、分布统计、OpenPIV单/双遍、CPD、齐次生成元、SciPy矩阵对数、因果窗口、旧状态构造和MockConsumer；历史相机、AFG和物理代码未修改。',
 '没有重写CNN主干、PIV、CPD、线性求解器或视频编解码器；新增的是三帧数据扩展、模型任务头、局部匹配与仿射适配、校准和科研评价。',
 '仅来自2026-10-15的同一段Chronos录像；train帧100–240，validation380–440，calibration520–580，final_test780/800/820/860，具体离散帧和四个不相交ROI见协议。没有跨录像验证。',
 '存在合成域限制：有限大小羽化图块、固定源图中的暗点、缺少真实三维出平面变化。图2公开这些痕迹；没有把阈值轮廓当成精确泡膜。',
 '源帧、空间素材区域和随机种子明确隔离；训练8192、验证1024、校准1024、最终测试1536。共同录像仍使独立实验泛化未验证。',
 '没有证明中心热图带来可用的六参数收益。少数剪切误差变化不足以通过联合风险门槛，H1输出还出现明显收缩到先验。',
 f'热图错误会传播：纯热图模型vx MAE从{m[ch["direct_heat"]]["mae_all_finite_global"][0]:.2f}增至严重扰动下{m[ch["direct_heat"]+"_severe"]["mae_all_finite_global"][0]:.2f}px/s。H4较稳定，但在预测塌缩下不能称为已解决鲁棒性。',
 f'纯灰度网络未可靠恢复：六参数MAE为[{nums(m[ch["direct_gray"]]["mae_all_finite_global"])}]，单位依次px/s、px/s和四项1/s。',
 'GeoMotion的若干平均误差略小，但预测变化几乎消失、分类器未区分非仿射/不可观测，联合风险仍失败；没有证明几何网络提高了可用可靠性。',
 f'GeoMotion8979参数，Direct6融合635271参数。纯模型p95分别{rt["geo_fused"]["pure_model_ms"]["p95"]:.2f}和{rt["direct_fused"]["pure_model_ms"]["p95"]:.2f}ms；完整管线p95分别{rt["geo_fused"]["full_pipeline_ms"]["p95"]:.2f}和{rt["direct_fused"]["full_pipeline_ms"]["p95"]:.2f}ms。小参数量不保证完整管线更快。',
 '没有。PIV接受1111/1536，其中84个全局案例超容差、1个非仿射误接受；AI无法在验证集找到满足5%风险且至少20接受样本的阈值，按协议全部拒判。',
 f'在25%诊断覆盖率384/1536下，PIV风险0.52%，灰度{pct(a["matched_coverage"]["0.25"][ch["direct_gray"]]["risk"]) }，融合{pct(a["matched_coverage"]["0.25"][ch["direct_fused"]]["risk"])}，GeoMotion{pct(a["matched_coverage"]["0.25"][ch["geo_fused"]]["risk"])}。曲线不是在测试集选部署阈值。',
 '平移和涡量最容易失败；预测接近小运动先验。小的散度/剪切MAE部分来自这些真值本来较小，不能单独证明恢复了形变。',
 '对称/重复分布使局部匹配有歧义；失焦降低位移证据；固定背景会与移动暗点竞争。条件误差和拒判分母均已保存，未把这些困难统统标成不可观测。',
 '正式AI因全拒判没有错误接受分母，不能据此宣称安全；诊断风险曲线显示大量错误。PIV仍有84个超容差接受和1个非仿射误接受，全部逐样本保留。',
 '策略拒绝全部162个不可观测样本，也拒绝全部1320个可观测全局样本，因此不是有效分类成功。融合网络原始分类识别118/162个不可观测样本；GeoMotion原始分类把全部样本判为全局。',
 '54个平滑非线性位移案例按预登记网格仿射残差>1.5px确定为非仿射；不是因为看起来相反运动就标非仿射。PIV判出51个、误接受1个、另2个因证据不足拒绝。',
 '逐参数实际覆盖率、平均宽度和困难类别已列于UNCERTAINTY_CALIBRATION.md。GeoMotion90%为['+nums(np.array(m[ch['geo_fused']]['intervals']['0.9']['coverage'])*100,2)+']%，95%为['+nums(np.array(m[ch['geo_fused']]['intervals']['0.95']['coverage'])*100,2)+']%。区间很宽，不能外推真实实验。',
 '已完成真实Chronos录像因果回放，但没有真实六参数真值，准确率仍NOT_VERIFIED。AI全部拒判，没有证据表明合成训练已可靠迁移。',
 'SynthPix预训练NOT_RUN。官方主要接口是图像对；本轮未验证连续三帧身份接口和独立JAX环境，优先完成核心12项训练及严格测试，不宣称有收益。',
 'NeuFlow v2对照NOT_RUN；本轮未下载/使用其权重，不用论文中的自然图像指标替代微泡验证。',
 f'正式训练12项、258epochs、11094步，墙钟{suite["elapsed_s"]:.2f}s；训练设备GPU采样平均{np.nanmean(gpu):.1f}%，峰值约4.97GiB。完整计时范围和缺项见RUNTIME_AND_RESOURCES.md。',
 f'是，融合AI完整管线p95={rt["direct_fused"]["full_pipeline_ms"]["p95"]:.2f}ms，小于同机单遍PIV{rt["PIV_SINGLE"]["full_pipeline_ms"]["p95"]:.2f}ms；这是短片软件回放，不是相机端到端测量。',
 f'没有。PIV新鲜有效六参数约{rt["PIV_SINGLE"]["fresh_six_hz"]:.2f}Hz，三种AI均0Hz。处理帧多不等于有效观测多。',
 '已保持兼容：新2.0Schema、显式旧1.0映射、逐字段null/valid/reasons、断帧/旧帧/TTL及真实本地ZeroMQ MockConsumer测试通过。动作始终null。',
 '七组英文PNG/PDF及三段实际数据MP4：真实检测热图、合成真值、六参数恢复、热图/风险覆盖、区间校准、资源性能、困难失败；视频为真实因果回放、合成真值对照和困难/热图消融。',
 'training/保存每项best/last权重、配置和曲线；results/保存逐样本与失败；logs/保留执行日志；README和ACTIVE_CODE_INDEX给出复现命令。大数据及周期checkpoint保留路径和SHA，不无选择提交。',
 '未解决可靠局部运动学习、真实非仿射判别、跨录像泛化、紧致且可靠的不确定性、真实像素/时间标定及硬件端到端延迟。PIV也未获通用物理准确准入。',
 '需要多个独立显微录像、已知载物台平移/旋转与可控二维形变参照、人工少量困难片段审核；随后只读Chronos在线验证。任何PZT响应实验都应另行预登记并校准声压与时间。',
 '没有充分证据推荐替换。保留PIV观测主线，保留本轮失败权重与接口；下一轮先在新开发素材验证可靠局部位移，再开新盲测，禁止继续据本轮final_test调参后重称独立测试。'
 ]
 questions=['F10U直接复用','是否重复造轮子','训练素材','合成伪影','数据隔离','热图收益','热图错误传播','纯CNN六参数恢复','几何网络收益','参数量与延迟','是否优于PIV','相同覆盖率风险','易失败参数','困难成像影响','错误接受','不可观测判断','非仿射识别','90/95%区间','真实域迁移','SynthPix','NeuFlow v2','资源和时间','完整管线速度','新鲜有效更新','M4兼容','图和视频','权重日志复现','未解决问题','下一轮实验','是否替换PIV']
 report='# F10V-Vision 中文科研结论\n\n**结论：已完成实际训练和独立验证；学习型运动观测器未通过科学准入，保留PIV主线。**\n\n'+'\n\n'.join(f'{i+1}. **{q}**\n\n{ans}' for i,(q,ans) in enumerate(zip(questions,answers)))
 write('FINAL_RECOMMENDATION.md',report)
 write('F10V_DEVELOPMENT_REPORT.md',f'''# F10V-Vision 科研开发报告

本轮完成了可执行的学习型群体运动观测原型、12项实际训练、验证选模、独立校准、1536样本最终测试和真实录像因果回放。**科研结论为不替换PIV**。代码运行完成与科学准确准入严格分开。

{table}
{risk}
{runtable}
训练共258个epoch、11094优化步、{suite['elapsed_s']/60:.2f}分钟；不是仅烟雾测试。最佳权重按冻结验证目标选择，两个GeoMotion主模型最佳点均在epoch1，后续验证恶化的证据没有删除。最终测试没有用于重新训练、调阈值或挑选随机种子。

所有AI因未找到合法验证风险阈值而拒判；零接受不是零错误成功。PIV的1111个接受中85个失败也完整保留。七组图、三段视频展示了实际数据与失败。只有一个录制会话、没有真实六参数真值；EXPLORATORY，PHYSICS_UNCALIBRATED，training_ready=false。相机、PZT、AFG及物理动作世界模型均未操作或训练。

30项逐问回答见[FINAL_RECOMMENDATION.md](FINAL_RECOMMENDATION.md)，分类/条件/种子/素材统计见[results/analysis.json](results/analysis.json)，源码入口与恢复方法见[ACTIVE_CODE_INDEX.md](ACTIVE_CODE_INDEX.md)和[README.md](README.md)。
''')
 write('README.md',f'''# F10V-Vision — Geometry-Aware Collective Microbubble Motion

Actual trained prototype; **not scientifically admitted to replace PIV**. Read [Chinese findings](FINAL_RECOMMENDATION.md) and [development report](F10V_DEVELOPMENT_REPORT.md). No physical device actions. `training_ready=false` concerns future physical world-model data, not whether these vision models were trained.

Code baseline:9cb22d729fece92fe160082edf09222a18068d74. Training snapshot7aed42436; frozen evaluationed36aa804; continuous demoa5a4ef53d. Full identities and source hashes accompany artifacts. Branch:dev/f10v-geometry-aware-motion-20261017. Local commit only; no F10V GitHub push.

Local module: `{ROOT}`. Remote: `/workspace/f10v_vision_20261017`, env `.venv`, data `data/`, output `runs/prototype_7aed/`, immutable code `snapshots/`. These folders are independent of the physics and old vision environments.

Quick review: `figures/Fig_F10V_1...7.*`, `videos/Video_F10V_1_real_causal_replay.mp4`, `Video_F10V_2_synthetic_truth.mp4`, `Video_F10V_3_difficult_ablation.mp4`. All plots are actual data; rejected raw diagnostics are clearly labeled.

Run in this module after activating the isolated environment:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python scripts/figures.py
.venv/bin/python scripts/videos.py
```

Full reproduction should use a **new output and cache directory** so frozen results are preserved. Restore the real-patch bank to the new cache; it is only about1.2MB and has per-file hashes. The raw video is needed only to re-extract banks or replay real microscopy. Regeneration without raw video:

```bash
export F10V_DATA_ROOT=/absolute/new/cache
export F10V_OUTPUT_ROOT=/absolute/new/run
# Copy data/banks and bank_manifest.json into F10V_DATA_ROOT first.
.venv/bin/python scripts/prepare_data.py --splits train validation calibration final_test --workers 8
.venv/bin/python scripts/train.py --all --batch 192
.venv/bin/python scripts/evaluate.py freeze
.venv/bin/python scripts/evaluate.py test
.venv/bin/python scripts/baselines.py --split final_test --workers 4
.venv/bin/python scripts/analyze.py
.venv/bin/python scripts/replay_runtime.py --video '/absolute/path/source.mp4'
```

The training script has a bounded budget and `--resume`; use the server Supervisor configuration retained in `evidence/server_execution/` for long runs. Never rerun `test` to tune a model; the script refuses to overwrite a completed test. Reproduction in a fresh output directory is an audit, not a new independent test. Canonical server labels and source hashes are retained; floating differences across environments are explicitly measured.

Offline inference:

```bash
.venv/bin/python scripts/infer.py --images frame0.png frame1.png frame2.png --pts 0 0.0333 0.0666 --model geo_fused --output /tmp/f10v_state.json
```

Current frozen thresholds reject all AI outputs. This is intentional evidence-based behavior; do not bypass it for hardware control. `--legacy` exports an explicit1.0-compatible message. All input images must be same ROI and use actual increasing timestamps.

Weights are retained locally under `training/*/best.pt` and `last.pt`; periodic checkpoints remain on the server. Weight hashes are in run summaries and `INPUT_OUTPUT_SHA256.json`. Virtual environments, full caches, third-party downloads and weights are excluded from Git; `ASSET_RESTORE.md` documents recovery. Some large text records are stored as `.json.gz` in Git and automatically restorable with the provided script.
''')
 write('ACTIVE_CODE_INDEX.md','''# Active code index

| Area | Entry | Purpose |
|---|---|---|
| Source reuse | src/legacy.py | Read-only F10U adapters |
| Data | src/dataset.py; scripts/prepare_data.py | Real patches, continuous triples, frozen split |
| Heatmap | src/heatmaps.py; src/learning.py | Actual detector input and H3/H4 perturbations |
| Geometry | src/geometry.py; tests/test_math.py | exp/log, fixed center, full resize conversion |
| Networks | src/networks.py | Official MobileNet blocks, Direct6, GeoMotion |
| Training | scripts/train.py; scripts/smoke_profile.py | Bounded full training and hardware profile |
| Selection | scripts/evaluate.py | Validation gate, independent calibration, locked test |
| Comparators | scripts/baselines.py | Frozen single/multipass OpenPIV and CPD |
| Metrics | src/metrics.py; scripts/analyze.py | Explicit denominators, risk curves, grouped analysis |
| Failure analysis | scripts/posthoc_diagnostic.py | Post-hoc prior-collapse description only |
| State | src/state_v2.py; scripts/infer.py | New schema and old-consumer adapter |
| Runtime | scripts/replay_runtime.py | Actual-PTS causal latest-slot replay |
| Visualization | scripts/figures.py; scripts/synthetic_demo.py; scripts/videos.py | Real figures and chronological MP4 |
| Reproduction | scripts/independent_check.py; tests/ | Independent checkpoint, gradients, math, MockConsumer |
| Audit | scripts/dependency_audit.py; scripts/report.py | Dependencies and result-derived reports |

Historical collective_motion_v1 and smart_online_tracking_v2 remain frozen. No Chronos receiver, PZT or AFG implementation is modified.
''')
if __name__=='__main__':main()
