# TopBrain → BraVa MCA 语义迁移

本轮已实现并运行语义点提取、Open3D 粗配准、POT Partial FGW、等权 donor 集成及真实 Pilot。**当前状态为 `TOPBRAIN_LABEL_TRANSFER_NOT_VALIDATED`，尚未进入 BG001 ROI 或 VascularMD 建模。**

## 方法与范围

TopBrain 严格掩膜存在数字拓扑环，未找到有依据的无损树转换；旧路线及专用环境已删除，原因见[简要历史](docs/research_history/RETIRED_TOPBRAIN_TREE_PIPELINE.md)。独立 CFD 表面扩展仍使用 VMTK，保留原状。

新路线只研究双侧 MCA 的 M1、M2、M3、UNKNOWN。TopBrain 提供原生解剖语义；BraVa 将提供树连接、坐标、半径及建模几何。配准只变换 donor，不应移动 BraVa。ACA/PCA 不做分段推断，原生标签注册表继续保留。

TopBrain M1 OR M2 OR M3 使用 Lee thinning 形成语义采样点，不建立图或输出 SWC。每侧最多 275 个代表点；物理网格内选择最接近均值的真实采样点，直接继承该点的原生标签。每个原始点保留映射，网格原始点数形成运输质量。抽样选择不读取类别，目标病例隐藏标签不会通过抽样进入拟合。

Open3D 仅做 similarity ICP，估计统一缩放、正旋转及平移。初始化使用 MCA/ICA 原生接触附近的根、90% 径向尺度、近端方向和四个固定轴向角度。无镜像、非刚性配准或内部 M1/M2 边界辅助。低质量、非有限或塌缩的 ICP 结果被拒绝。规则在配准 JSON 中固定记录。

POT 使用原生 `ot.gromov.partial_fused_gromov_wasserstein`。C1、C2 和配准后跨集合距离 M 均为欧氏距离，分别除以正距离中位数。运输矩阵列中按 donor 原生标签累加并归一化为三类概率。列质量低于正列质量中位数的 0.25 倍时无有效支持。最多选择 5 个配准合格 donor，各点对提供有效支持的 donor 等权平均；有效支持 donor 少于入选数一半时标记 UNKNOWN。没有按分支代次或远端位置补标签。

## 环境与运行

从项目根目录执行：

```bash
PY=/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python
"$PY" tools/setup_topbrain_brava_transfer.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$PY" tools/validate_topbrain_transfer.py \
  --topbrain-root data/TopBrain --pilot
```

固定 Open3D 0.20.0、POT 0.9.7.post1；本机使用主环境，已有包版本未变。Open3D 缺失的 libusb 使用已校验的 conda-forge 二进制，放在项目输出目录的 `runtime/lib`，由 `transfer_runtime` 在导入 Open3D 前加载；没有 sudo、系统 Python 修改或源码编译。部署脚本已重复执行验证，不重复安装正确版本。环境记录见 [environment_summary.json](outputs/topbrain_brava_transfer/environment_summary.json)。

`--cache-only` 只建立语义缓存，`--debug-save-matrices` 才额外保存 C1/C2/M/T。`--full-validation` 是显式可选任务，本轮未运行。现有结果可以直接重新判读并生成点云 VTP，无需重跑配准与运输：

```bash
"$PY" tools/validate_topbrain_transfer.py --evaluate-existing
```

退出码 2 表示验证未获准进入生产，不表示包安装失败。

## 真实 Pilot 结果

目标为 MRA001–005 的双侧，共 10 个侧别。缓存来自 25 个病例的 48 个可用侧别；010/R 和 014/R 缺原生分段。每次排除目标病例，严格使用同侧 donors，入选数为 3–5。评价在全部原始语义采样点进行，各侧 F1 等权平均；UNKNOWN 计入对应真实类别的假阴性。

| 方法/参数 | M1 F1 | M2 F1 | M3 F1 | macro F1 | UNKNOWN |
|---|---:|---:|---:|---:|---:|
| Open3D + 最近邻基线 | 0.7413 | 0.8012 | 0.7820 | 0.7748 | 0.00% |
| 默认 α=0.5，m=0.8 | 0.7250 | 0.8315 | 0.6301 | 0.7289 | 17.54% |
| α=0.4，m=0.75 | 0.7030 | 0.8111 | 0.5790 | 0.6977 | 22.55% |
| α=0.4，m=0.90 | 0.7378 | 0.8377 | 0.6895 | 0.7550 | 8.85% |
| α=0.7，m=0.75 | 0.7409 | 0.8225 | 0.5605 | 0.7080 | 21.78% |
| 最佳 α=0.7，m=0.90 | 0.7668 | 0.8434 | 0.7031 | 0.7711 | 7.45% |

最佳结果的 M1/M2、M2/M3 边界中位误差为 **3.527 mm、3.690 mm**，两项均在 10 个侧别可计算。这是相同局部几何近邻图上真实/预测界面的对称最近距离：每侧先取中位数，再在侧别间取中位数，不等价于体素表面的精确分段误差。

最佳参数满足三项绝对 F1 门槛，但 macro F1 比基线低 **0.00373**，M3 F1 也低于基线。按本次要求第 100 节“不优于 baseline 就停止”，最终门禁要求 FGW macro F1 严格高于基线。最初运行日志使用了允许 0.02 差值的宽松判据；该判据已经纠正，最终报告仅重新评价原有数值，没有重跑或改变预测。最终以 [pilot_summary.json](outputs/topbrain_brava_transfer/pilot_summary.json) 和 [final_gate_review.log](outputs/topbrain_brava_transfer/final_gate_review.log) 为准。

该小网格在同一 Pilot 上选择参数，没有独立测试集，因此即使通过也只能视为后续概念验证的入口。本轮在上述停止点结束，没有增加参数组合或引入其他算法。

## 查看结果与后续导出边界

每个目标侧别的 `outputs/topbrain_brava_transfer/pilot/MRA<case>_<side>/semantic_validation_points.vtp` 可在 ParaView/PyVista 中打开，切换 `predicted_label`、`truth_for_evaluation_only`、`p_M1/p_M2/p_M3`、`confidence` 或 `support_mass` 着色。它是语义点云，没有中心线树连接；真值字段仅用于拟合后的评价。完整原始点概率保存在对应 `*_predictions.npz`。

BG001 的 transfer CLI、分支聚合、新 `--source topbrain-brava` UI 模式、ROI SWC 导出和 VascularMD 调用没有进入开发或执行阶段。没有生成伪造成功的 VTP、概率 CSV、SWC、VTK 或 STL。既有 `s1-3_swc_roi_generate_MeVO.py` 仍可检查原生 MCA 掩膜，退休的 TopBrain→SWC 选项已移除，成熟渲染器及 s1-2 保持不变。

只有验证门禁通过后，才应实现 BraVa 原始点/分支映射、按弧长聚合、标签逆序检查、UNKNOWN 间断与 ROI 连通性检查；合格 ROI 仅能保留原始子图，之后通过现有 `vascular_processing.pipeline.process_file` 调用 VascularMD。当前没有可供调用的真实合格 ROI，因此这些步骤均为未执行，不能称为已建模或已制造。

## 验证与完整性

清理后稳定回归 **141 项通过**，覆盖原始数据、标签映射、affine、BraVa 解析/拓扑、SWC、VascularMD 适配、界面和渲染辅助。新增核心测试 **8 项通过**，覆盖粗配准恢复、禁反射/无效变换、软运输概率、低支持 UNKNOWN、等权集成、远端额外结构、抽样映射和保护校验、基线停止规则。新 BraVa 聚合测试属于未进入的后续阶段。

268 个受保护文件无变化，包括 100 个 NIfTI 和 123 个原始 BraVa SWC；VascularMD upstream、mouse 和成熟 renderer 未改。s1-2 SHA256 仍为 `3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729`。

删除清单、环境、逐目标结果及完整交付状态见 `outputs/topbrain_brava_transfer/`。官方接口依据为 [Open3D PointToPoint](https://www.open3d.org/docs/release/python_api/open3d.pipelines.registration.TransformationEstimationPointToPoint.html) 和 [POT Partial FGW](https://pythonot.github.io/gen_modules/ot.gromov.html#ot.gromov.partial_fused_gromov_wasserstein)，实际 API 签名与版本另已从本机安装包记录。
