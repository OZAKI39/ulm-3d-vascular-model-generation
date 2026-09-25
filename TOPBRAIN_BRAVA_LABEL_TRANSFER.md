# TopBrain → BraVa 语义迁移

正式方法为原 Open3D similarity registration + 最多五个同侧 donor 的 NN 投票。旧 MeVO 比较结论固定为 `MEVO_BINARY_NN_PREFERRED`，NN F1=0.977445；旧比较 JSON/CSV/Markdown 均保留。Partial FGW 实现、参数搜索、专用测试和 POT 依赖声明已经退休，详见 `docs/research_history/RETIRED_PARTIAL_FGW_TRANSFER.md`。现有环境无需安装或卸载软件。

十侧 NN 硬预测及其逐 donor 投票、距离、注册标识已冻结在 `outputs/topbrain_brava_transfer/frozen_nn_pilot_predictions/`。四组 UNKNOWN 校准已完成且缓存。没有组合同时满足降低 M1 错纳率和 F1 降幅不超过 0.01，因此选择 agreement ≥0.60、距离参考分布第 99 百分位、至少三个合格 donor，状态为 `OOD_SUPPORT_GATE_ONLY`。这是支持范围拒绝，不宣称改善分类。

用户已明确选用官方 BG001 长度表推断的对应关系：**TYPE 3 = LMCA，TYPE 4 = RMCA**。`configs/brava_major_arteries.json` 记录 `USER_ACCEPTED_OFFICIAL_TABLE_DERIVED`，仅适用于所绑定的 BG001 两个源文件。原有来源冲突作为证据保留，不宣称获得了官方直接 TYPE 字典。来源证据见 `references/文本记录 - BraVa ColorCoded TYPE 来源追踪.md`。

离线复核入口：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tools/trace_brava_type_mapping.py
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tools/validate_topbrain_transfer.py --freeze-nn --calibrate-unknown
```

第二条命令在现有缓存下只读验证，不重复配准、NN 查询或四组校准。`tools/setup_topbrain_brava_transfer.py` 现在仅离线检查已安装的 Open3D 0.20.0 和本地 libusb，不下载或安装。

真实 BG001 结果：LMCA 的 25 个 donor 配准均未通过原有 QC，626 个点、39 条分支全部记录为 UNKNOWN，状态为 `REGISTRATION_UNSUPPORTED`；RMCA 的 5 个 donor 通过检查，得到 44 条 MeVO、8 条 M1、3 条 UNKNOWN 分支。44 条 MeVO 分支构成 5 个独立 ROI，总计 687 个原始节点。

RMCA 的 part01、02、04、05 均通过原生 VascularMD `model_network` 并输出 smooth SWC、VTK 和 STL。part03 的根为分叉且没有相邻上游 M1 context，原有输入检查拒绝建模；严格 SWC 仍保留，未改造入口。几何是 `ANATOMICAL_TRANSFER_CANDIDATE`，制造状态为 `MANUFACTURING_NOT_VALIDATED`。

生产与显示入口：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tools/transfer_topbrain_to_brava.py
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python s1-3_swc_roi_generate_MeVO.py --source topbrain-brava
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python s1-3_swc_roi_generate_MeVO.py --source topbrain-brava --brava-view surface
```

默认右视窗读取 strict 原始 ROI；surface 模式只列出成功建模组件，表面包含已注明的真实上游 M1 context。左右箭头切换 ROI，A/R/S/C 保留框选显示，F12 截图。UI 不做配准或标签查询。s1-2、共享 renderer、原始数据和冻结的配准实现保持不变。输出位于 `outputs/topbrain_brava_transfer/nn_production/BG001/`，原始输入和 donor 缓存由 `input_provenance.json` 绑定，完成的结果可直接复用。
