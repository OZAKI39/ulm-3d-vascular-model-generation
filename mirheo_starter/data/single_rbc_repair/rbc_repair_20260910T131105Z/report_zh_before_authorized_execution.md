# HemoCell / Mirheo 单红细胞剪切流修复核查

本轮未完成同等质量、同一终点的完整比较，qualified_speedup=null。已完成 CPU 排查和最小修复候选；按用户第 2、5 节要求，新的 GPU 与隔离编译等待集中授权。

- 直接退出原因已确认：粗候选超过 6400。写入有边界保护，不等于越界写。候选为什么增长到这个数量仍缺出错一步的完整状态。
- 更早的异常已核实：默认主运行分别在准备第 2500、2000 步出现成员不符；半步在准备第 5000 步出现。独立检查确认距离约 0.1633、0.04188、0.3460，均远大于 1e-4 近膜带。这排除了阶段交接作为全部异常的唯一原因。
- 受限队列诊断：Γ=2.2 首个坏存帧有 21 对真实自交、成员不符，最大边长为 WLC 上限的 1.652 倍；前一正常存帧 Γ=2.1。不能凭这两个存帧判断具体哪个异常先发生。
- 调度候选：local 与 halo 使用同一 bouncer 暂存表而缺少执行依赖，空 halo 也会清零缓冲。源码确认风险；实际失败因果尚待 GPU 对照。隔离补丁增加依赖并跳过空 halo；旧库未替换。
- 材料尚不匹配：初始力 RMS 为 HemoCell 3.865e-10、Mirheo 1.768。独立 Kantor 弯曲能量梯度重现 Mirheo 初始总力，误差 0.14%。原 17.82% 伸长差、约 7.7% 准备形状差、6.47% 黏度半步差仍是失败证据，未用新阈值抹去。弯曲、释放、惯性和共同准备质量尚待新检查。
- 本轮新 GPU 求解、CPU 求解和编译均为 0；没有新的完整耗时。历史 HemoCell 58.17194、60.78252 秒，均完成 Γ=4；历史 Mirheo 四次失败进程分别花费 147.76706、620.42071、355.33590、786.97269 秒，不能据此计算合格速度比。

第一批真实待执行队列：四个各最多 90 秒的零剪切短对照。隔离编译按需在第三个对照前进行。--execute 只执行这一冻结短队列；之后必须先审查实测数据，才能准备后续更长测试。没有伪装成已实现的无人干预全流程。

集中申请上限：GPU 8500 秒，A=1600 / B=1900 / C=5000；CPU 求解 600 秒；隔离编译 1800 秒。单 GPU、单任务最长 2400 秒，失败计费、无自动重试，前关失败阻止后关昂贵运行。旧剩余 788.296746 秒保持旧账本。估算只用于是否可启动的预算判断，不是完成 Γ=4 的承诺。

实际命令：

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.repair_single_rbc_benchmark --help
.venv/bin/python -B -m py_scripts.repair_single_rbc_benchmark --config py_scripts/single_rbc_benchmark_repaired.yaml --preflight-only
.venv/bin/python -B -m py_scripts.repair_single_rbc_benchmark --config py_scripts/single_rbc_benchmark_repaired.yaml --review
# 仅在新授权正式记录后；执行冻结的短诊断队列：
.venv/bin/python -B -m py_scripts.repair_single_rbc_benchmark --config py_scripts/single_rbc_benchmark_repaired.yaml --execute
```

离线页面：../../../test_code/outputs/single_rbc_repair/rbc_repair_20260910T131105Z/comparison_review.html。页面上部为本轮 CPU 结论，下部明确标为旧 campaign 的完整离线页面，复用旧真实膜动画、D/倾角/A/V、计时和交互代码；不伪装成修复后的新轨迹。

原始 CSV、NPZ 和日志继续保留在 runs/single_rbc_benchmark/rbc_shear_20260910/solver/{main_mirheo_1,main_mirheo_2,strict_mirheo,diagnostic_mirheo_full}；精确路径及哈希在 failure_timeline.json。本轮派生 CSV、首个三角面对及粒子坐标、CPU 能量核查和补丁均在当前 data 目录。新 runs 目录尚无原生求解记录。

CPU 检查、浏览器检查、人工验收分别记录；人工保持 PENDING。本轮未提交或推送 Git。
