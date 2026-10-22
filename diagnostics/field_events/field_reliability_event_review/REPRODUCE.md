# 原事件区域背景场诊断复核

本目录是 2026-10-05 的独立诊断交付。先打开 `OPEN_RESULTS.html`、`FIELD_RELIABILITY_REPORT_ZH.md` 和 `FIELD_RELIABILITY_RESULTS.json`。这里没有 CFD 求解、网格生成或微泡积分入口。参考细场不是解析真值；结论是受控数值差异及有限特征支持。

## 输入与执行位置

本地工程根目录为 `/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin`；现有服务器工程根目录为 `/workspace/brava_flow_roi_18mlmin_20260928`。服务器连接使用既有私有配置，本目录不保存连接信息。服务器操作前读取 `/etc/vast-agents-guide.md` 并重新核查实际资源。

原输入相对工程根目录：

- `cases/balanced_pressure_final/frozen_flow/steady_flow.vtu` 与 `mesh/SV_MESH/mesh-complete.mesh.vtu`。
- `population_v1/analysis/speed_event_cause_review/` 中四个 `original_<ID>.npz`、事件边界及原始证据。
- 两个参考网格：`mesh_particle_validation_v2/02_mesh_candidates/M1_fixedplc_{medium,fine}/mesh/SV_MESH/`。
- 两个参考解：`mesh_particle_validation_v2/03_cfd_runs/M1_fixedplc_{medium,fine}/frozen_flow/steady_flow.vtu`，及其真实原生终态、身份、求解器输入和稳态记录。

逐文件身份见 `INPUT_SHA256.json`，全部边界、物性、历史 BC、求解实现、稳态和实际点覆盖见 `FIELD_COMPARABILITY.json`。大体积参考原始载荷保留在现有服务器位置；本地缺失不能通过忽略检查或用图像替代载荷解决。

依赖现有 `population_v1/env/bin/python`，以及原生产 `population_v1/source/particle_3d/src/particle_3d` 的 CUDA FP64 查询和既有扩展缓存。使用 numpy/scipy/pyvista/VTK/torch/matplotlib/Pillow/pytest/psutil，不安装新环境、不修改生产源代码。CPU BLAS 线程由既有监督器限制为 1；逐场加载，未复制全网格并行竞争内存。

## 冻结协议和重做规则

本次协议 SHA256：`9e1aa92438e86dadaddc79c5d837d73fe558af447cf4d986d56bbdc2d54a5aea`。它在读取参考速度载荷前冻结。本次交付不能原地重跑覆盖。若今后明确授权复核，创建本目录的独立同级复核目录，仅复制脚本及冻结协议/协议源 SHA；不要复制分段缓存后冒充重新枚举。保持目录深度以让 `common.py` 正确定位工程根目录。`freeze_protocol.py` 留作本次协议生成过程的证据，不应重新挑选已冻结事件。

在服务器工程根目录、已创建的独立复核目录中，使用以下顺序。这里的 `review_dir` 需指向新的同级目录，每条命令以前一条正常完成和检查输出为前提；不得改为无人监管自动重试队列。

```bash
review_python="$PWD/population_v1/env/bin/python"
review_dir="$PWD/population_v1/analysis/field_reliability_event_review_recheck"
"$review_python" "$review_dir/supervise.py" admit_fields
"$review_python" "$review_dir/supervise.py" prepare_coordinates
"$review_python" "$review_dir/supervise.py" sample_fields normals
"$review_python" "$review_dir/supervise.py" sample_fields M0
"$review_python" "$review_dir/supervise.py" sample_fields M1_fixedplc_medium
"$review_python" "$review_dir/supervise.py" sample_fields M1_fixedplc_fine
"$review_python" "$review_dir/supervise.py" analyze
"$review_python" "$review_dir/supervise.py" spatial_features
"$review_python" "$review_dir/supervise.py" plot_results
"$review_python" "$review_dir/supervise.py" run_tests
"$review_python" "$review_dir/supervise.py" write_report
```

`supervise.py` 复用现有 `production_workflow/scripts/managed_job.py`，创建本目录专有 `monitoring/<stage>/` 的日志、LIVE 和 EXECUTION。累计诊断预算 5400 s，单段不超过 1800 s，进程树 RSS 10 GiB、cgroup 非可回收内存估计 34 GiB；采样结束还核查 CUDA 张量分配峰值 <8 GiB，该值不是设备总占用或连续 GPU 强制限额。停止应向本目录监督器发送终止信号，由它终止所观察的项目进程树，不按全局进程名杀进程。

本次失败尝试及修正包含在 `EXECUTION_RECEIPT.json`。VTK 包围盒查询三次停滞后改用保守 cKDTree 候选加原重心裁剪；未改几何。早期 M0 控制体出现未处理的严格边界包含，退出 0 也未作为完成结果；最终使用原面精确 P1 trace 和不移点的扩展精度包含核查。空间特征初次键类型错误另行修正。全部失败耗时计入预算。

## 数据字典和图像复核

- `path_<ID>.npz`：全部接受点、原时间、弧长；各场过面分段的并集；Gauss 1/2/4 阶共点及弧长/时间权重。`partitions_<field>_<ID>.json` 保留各自单元内分段，不跨网格匹配 cell ID。
- `plane/profiles/rise_<ID>.npz`：实际 SI 坐标、固定局部基底及有效范围。法向线遇原有限面片投影无效即截断。三例近壁事件的辅助 A–B 直弦也保留，但主结论使用原折线；822 的 A–B 本身是一条真实接受线段。
- `<field>_<kind>_<ID>.npz`：`valid`、原坐标处所属 cell、`flow` 和 `cell_geometry`。`flow` 列为 u[0:3]、p[3]、G 行优先[4:13]、curl[13:16]、strain[16:25]。几何列为最短/最长边、最小/最大面高、体积。全部 SI；NPZ 的无效数值配合显式 mask，JSON 用 null，禁止将域外值外推。
- `wall_<kind>_<ID>.npz`：原有限壁面最近特征、中心距 d、原法向和歧义；球面 gap=h=d−a。`WALL_NORMAL_TIES.json` 补充接受点和法向线的全部并列面片；其它采样点未逐个完整枚举所有并列面片。
- `cv_<ID>_L<0/1/2>.npz`：相同闭合物理表面的所有三角形、原节点、外法向、面积、壁面标记及体积；`<field>_BALANCE.json` 保存逐面通量和积分细化历史。throughput=0.5Σ|原三角面的通量|，保留面分组，不加 epsilon。
- `FIELD_METRICS.json`、`SPATIAL_FEATURES.json`、`representative_values.csv`：图表源值、共点各阶段统计、独立切面统计、逐单元有符号增速贡献、822 固定阈值穿越。
- `vtk/`：原壁面、事件 ROI、实际场查询邻域和闭合面，可直接用 ParaView 打开。u 为原节点 P1，G/div/curl 为原单元量，不做显示节点平均后冒充粒子输入。
- `FIGURE_METADATA.json`：切面、原点平移、相机、共同色标、箭头比例及数据源。五幅 PNG 均至少 300 dpi，英文图内文字、中文报告；`VISUAL_QC.json` 记录人工看图核查，不把图像平滑当科学门槛。

`TEST_RESULTS.json` 是实际测试输出；测试范围是诊断实现。`FINAL_INPUT_INTEGRITY.json` 同时核对本轮登记输入和上轮 195 项原输入。`OUTPUT_SHA256.json` 可复核交付文件运输完整性（不含它自身、监督运行状态和缓存）。最终 `finalize_delivery.py` 仅在所有受监督阶段结束后运行，检查活动任务、汇总真实 receipt、生成哈希，不做新科学计算；图像有变动时必须重新人工复核后更新 `VISUAL_QC.json`。

所有新任务现已正常结束或明确停止；没有自动下一档。建议仅为取得同几何/同历史边界、局部误差受控的事件区域参考，本轮没有启动该后续。
