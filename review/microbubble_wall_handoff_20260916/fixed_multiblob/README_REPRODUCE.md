# 重现与核查入口

本目录是独立科研审查归档，未替换正式 wall V0。首先阅读主报告中的失败范围、未运行范围和参考资格限制。原始矩阵与失败记录都保留。

## 只重算结果（不运行水动力求解）

环境版本在 `source/requirements-lock.txt`。归档未捆绑 Python 虚拟环境；工作目录中的隔离环境可直接使用。

```bash
cd /home/lzy/projects/compre_output/fixed_multiblob_wall_feasibility_audit/20260916_143602
sha256sum -c SHA256SUMS
cd remote_raw
sha256sum -c SHA256SUMS
cd ..
OPENBLAS_NUM_THREADS=1 /home/lzy/projects/fixed_multiblob_wall_feasibility_audit_20260916_143602/env/bin/python scripts/finalize_fixed_multiblob_wall_audit.py
/home/lzy/projects/fixed_multiblob_wall_feasibility_audit_20260916_143602/env/bin/python scripts/visualize_fixed_multiblob_wall.py
```

重算会更新派生 H5/CSV/图；原始 `remote_raw/` 不变。重算后应生成新的派生结果清单，不能把旧清单继续称为当前文件的校验结果。finalizer 不导入主求解器，不读取求解器 PASS 标签作为数值依据。

## 源码与输入

- `source/upstream_stokesian-dynamics/`：固定提交的上游原始源码与 MIT license。
- `source/pecnut_runtime/`：相同数学源码，加入原生生成器生成的所需尺寸比表。
- `src/pecnut_fixed_wall_adapter.py`：规定墙运动、矩阵缩放、标准 Schur 消元、缓存和确定性几何。
- `scripts/run_remote_audit.py`：有界静态参数扫描和原始 H5 写入。
- `scripts/remote_control/`：本轮实际 SSH、supervisor 与预算控制脚本。端点是本轮实例，不能假设以后仍存在；程序拒绝覆盖已有运行。
- `contracts/`：冻结科学定义、主资源上限、系数生成计划及 2,000 球实测内存预检依据。
- `input/`：冻结 RMBW 表、原 STL、此前 local-plane 指标与最终生成的系数表。
- `raw/SCALAR_GENERATION_*_FIRST.npz`：每次新增尺寸比的原生 near/mid 输出及合并表；继承列没有重加权或插值。
- `provenance/RUNNER_INITIAL_SOURCE.py`：最早的调度版本；最终脚本添加了输入身份记录与小规模原生一致性门槛，水动力核未改变。

当前生成表对应的本轮 β 为 .5、.25、.125；不要直接调用未生成尺寸比的原生 lookup。更细分辨率已被资源预检排除。

新运行必须另建输出目录，先验证输入 SHA、可用内存、源码身份与小规模原生一致性；不允许覆盖本归档。低精度替代、简化边界、拟合墙高或缩小 patch 来冒充同尺度细化都不属于本 contract。

运行入口采用 `AUDIT_ROOT`（含 src/scripts/contracts）、`AUDIT_RESULT`（新的输出目录）、`PECNUT_CODE`（含 settings.py 的原生运行目录），CPU BLAS/OMP/Numba 线程均为 4。单个文件级 runner 有资源检查；跨阶段总时长由本轮 supervisor 外层控制，不能把每个文件的上限各自重复消费当成无限预算。

## 时间、矩阵和可视化含义

矩阵、动作向量与耗散使用 work-conjugate bulk-drag 缩放，DOF 顺序为 Vx,Vy,Vz,Ωx,Ωy,Ωz。H5 中另存冻结三种半径对应的 SI 转换。有效查询保存目标球力/力矩响应；不声称输出了完整墙面 traction 场。

`PERFORMANCE_SCALING.csv` 的计时来自实际 CPU，包含未执行和查询失败行。2,000 球装配完成而查询失败；首个保守预检拒绝也保留。查询失败时不填写成功吞吐。内存字段是进程峰值 RSS，包含装配临时量，不等于单个矩阵大小。原生 pair solver 主对照采用一个活跃近场对；cutoff 边界控制结果不用于主速度比。

泄漏 CSV 同时保留 signed transmission 与 amplitude ratio；方向反转不是“负渗透率”，也不能利用负号通过层数改善门槛。

`PLANAR_WALL_BLOBS.vtp` 坐标及 blob_radius 为米，按 d50 尺寸显示。object_type=0 表示固定墙球；bubble_id=-1 表示墙球；未定义的单球 gap/matrix_error 用 NaN 标记，validity_flags=0 表示未获生产资格。其余三个 VTP 为空并带 NOT_RUN 状态，因为曲壁/真实 STL 门槛未通过。不得把空文件解释为实际几何建模失败或零水动力误差。

HUMAN_VISUAL_REVIEW=PENDING。NEXT_STAGE 仅为建议；本目录不启动其他科研阶段。
