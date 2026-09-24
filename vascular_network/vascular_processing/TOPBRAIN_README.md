# TopBrain 原生数据与 MCA 显示

保留官方 2025 v2 数据发现、ITK-SNAP 标签解析、完整 NIfTI affine、来源哈希、严格 M2 OR M3 掩膜及双视窗。TopBrain 不再导出中心线树或 SWC。ACA/PCA 仅保留原生标签名称注册，不再执行合并标签拆分。

数据入口为 `data/TopBrain`；原始数据、License.txt 和 `outputs/topbrain_validation/dataset_provenance.json` 均保留。来源为 [Zenodo 16878417](https://zenodo.org/records/16878417)，DOI 10.5281/zenodo.16878417，压缩包 MD5 为 b703ea31cd1f0e7115a5d3e6e61f59b3。

```bash
python s1-3_swc_roi_generate_MeVO.py --discover
python s1-3_swc_roi_generate_MeVO.py --case-id 001 --modality mr --roi LMCA_M2M3
python s1-3_swc_roi_generate_MeVO.py --case-id 001 --off-screen
```

`--no-gui` 保存原生掩膜与表面，`--show-native-labels` 叠加原始分段边界，`--smoke-gui-seconds` 执行定时界面检查。原来的布局、相机、Trackball 和快捷键继续复用共享显示模块。缓存校验原始数据、标签表、程序及产物哈希；`--recompute` 创建新的结果目录。

退休原因见 [简要历史](../docs/research_history/RETIRED_TOPBRAIN_TREE_PIPELINE.md)。新的生产候选路线使用 TopBrain 语义迁移至 BraVa 原生树，数据与显示工具不参与改变 BraVa 几何。
