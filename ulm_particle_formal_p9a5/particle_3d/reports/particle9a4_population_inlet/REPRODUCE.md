# P9-A.4 reproduce / review

本轮只在新worktree开发，没有push/merge。原worktrees只读。所有下列命令从 `/home/lzy/projects/ulm_particle_population_inlet_p9a4` 执行。Python 3.13.11、NumPy 2.5.3、SciPy 1.18.1、PyVista 0.49.0；服务器Python3.13.15、相同NumPy/SciPy/PyVista。记录源/输入SHA而不假设跨任意软件版本bitwise一致。

```bash
cd /home/lzy/projects/ulm_particle_population_inlet_p9a4
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONPATH="$PWD/particle_3d/src"
export SONOVUE_ROOT="$PWD/sonovue_size_distribution_v0"
P9A4_PY=/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python
P9A4_CHECKDIR=$(mktemp -d /tmp/p9a4-review-XXXXXX)
"$P9A4_PY" sync_metadata/network_h0_particle_20260925/run_checks.py --output "$P9A4_CHECKDIR/portable"
"$P9A4_PY" -m pytest -q -p no:cacheprovider --junitxml="$P9A4_CHECKDIR/new_legacy.xml" particle_3d/tests/particle9a4_population_inlet particle_3d/tests/particle9a2_inlet particle_3d/tests/particle8_2a/test_open_inlet_not_solid_wall.py
"$P9A4_PY" sync_metadata/network_h0_particle_20260925/verify_snapshot.py --scientific-inputs
```

`run_checks.py --output` 必须绝对路径；其子进程工作目录为临时副本。该旧工具将5份路径manifest只在临时副本重定位，并只在临时源码保护合同中排除92条生成pyc，所有原科学源SHA仍强制匹配。新tests在独立包命名空间，避免多conftest重名。

## 入口生成与checkpoint

本次已经完成：`scripts/run_inlet_audit.py prepare`（独立20k Qacc + immutable contract）、`generate --workers 6 --count 100000 --label inlet100k`、`replay --workers 1`、`replay --workers 3`。不要在已有合同上重跑prepare，不要覆盖现有receipt/ledger。重新生成入口审计使用新label，例如：

```bash
"$P9A4_PY" particle_3d/reports/particle9a4_population_inlet/scripts/run_inlet_audit.py generate --workers 6 --count 100000 --label inlet100k_review_repeat
```

只生成proposals与accepted births，不跑轨迹。新label如果存在会显式拒绝。需要检验1/3/6而不写既有receipt，可在Python中加载同一冻结contract并调用 `make_source(load_new_environment(root), contract)`、`generate(source,100000,workers=n)`，对每行 `canonical_bytes` 的串接计算SHA，与 `data/worker_replay_1.json` 比较。source ledger SHA预期 `b89b019874b8d142f31f227ce45ed20f2d0b5301595b1b8dc3dd92cdad3c8cd7`。生产源身份包含两份新源码SHA；源码若改变不能假装原ledger仍是同一科学版本。

`PopulationLedger.restore(folder, source.identity)` 读取完整manifest、ledger和state，验证后得到同样next IDs/time。`generate(source, extra_count, ledger=restored, workers=n)` 从下一个source ID延续，保存到全新checkpoint目录。永久测试覆盖53个事件（包括拒绝）checkpoint后75个事件续跑，对照128个不中断事件。缺失/篡改/identity不同的checkpoint拒绝，不能删除拒绝行来减小文件。

## 实际输入与保护

NEW输入：`particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle/inputs/NEW.vtu`。
SHA：`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。
source合同：`particle_3d/contracts/P9A4_CONTINUOUS_INFUSION_V1.json`。
原科学源码保护：`particle_3d/reports/network_derived_flow_mb_validation_v1/data/code_snapshot_manifest.json`。

两份新增科学模块直接使用原 full-positive-flux sampler、原 conditional SonoVue sample、原 admission.check。原 legacy B/C、STREAMS、FrozenFEMField、WALL builder、P9-A.1 equations、OLD contract均未编辑。旧轨迹checkpoint由原版本重放，不为新增module放宽其科学身份。

## 审核脚本与输出

`gate_inlet.py` 读取data中的既定宽容差、100k及独立控制样本。真实审计容差在gate执行前写入，但不是在所有初步aggregate计算前盲定；synthetic永久测试容差预先固定。`compare_legacy.py` 的500表示入口请求，不是trajectory；比较accepted分布时应保留full vs conditional source与分母差异。

`run_smoke30.py` 只允许 `/workspace/particle9a4_population_inlet_<UTC>`，核对deployment manifest、gate、固定30cohort与NEW身份后创建独占输出目录。它复用原Network mb_job；原积分器启动时再次验证既定birth一致性，不重抽source。`deploy_smoke.py` 本轮创建的remote_root见data/server_deployment.json。已完成30 smoke，不要为审核图表重复运行或把它扩大为正式500。服务器操作遵循已保存的logs/server_agent_guide.txt。

`analyze_smoke.py` 只读保存结果，复用原analyze.one并额外独立壁距/入口相交审核，不推进轨迹。`render_figures.py inlet/smoke` 只读数据，英文白底300dpi PNG+PDF；会拒绝覆盖已有图件。`build_reports.py` 从审计产物构建报告，也拒绝覆盖最终交付；重建草稿需另存目录/版本。原阶段输出永不覆盖。

完整100k ledger约88 MiB，是重要科学证据；全部reject行保留。服务器原始30轨迹、metadata、压缩audit、point结果均已回传并保存SHA。报告支持的结论来自这些完整数据。图件目录figures/，主报告和机器摘要见OPEN_RESULTS.html；不需要旧可视化目录参与本轮。

## 本轮执行结果

baseline115 passed（87.810 s）；smoke前portable115、新35、legacy15通过；最终portable115 + new35 + legacy15=165通过、0fail/skip。完整100k worker1/3/6逐字节一致；30 smoke28completed/2supported stationary/0solverfail。46,670旧文件SHA不变。失败的开发测试/辅助脚本日志保留于logs，原因见BUG_FIX_LEDGER_ZH.md。没有正式500、没有push或merge。
