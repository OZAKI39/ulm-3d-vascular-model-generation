# 3D 脑血管打印模型：当前代码与结果快照

本分支同步截至 2026-09-25 的 TopBrain → BraVa BG001 RMCA → MeVO 语义迁移、ROI 精修、打印候选、端口与浇注盒生成代码，以及关键数据、日志、QC 和测试记录。源工作区的完整业务代码与所选结果逐字节保留；同步不改变血管模型、半径或可视化风格。

**当前限制：O3 接缝的光顺性仍待进一步修复。** `print_fixture_design_aligned` 已验证端面轮廓与轴向对齐，但这些检查不证明侧壁切向连续。已有 3MF 中血管芯与盒体仍为两个实体；一体式 union 阶段尚未实现。本快照不宣称模型已经通过打印或制造验收。

## 代码与关键结果

| 用途 | 入口／目录 |
| --- | --- |
| 小鼠 SWC、官方逐分支 AIC 派生数据 | `s1-1_swc_roi_generate_mouse.py`、`tools/vascularmd_prepare_mouse.py` |
| 人脑 BraVa SWC 与 ROI | `s1-2_swc_roi_generate_human.py` |
| MeVO／打印 ROI 双视窗 | `s1-3_swc_roi_generate_MeVO.py` |
| TopBrain 数据准备、语义迁移 | `tools/topbrain_prepare_dataset.py`、`tools/transfer_topbrain_to_brava.py` |
| RMCA 精修、制造候选、紧凑候选 | `tools/refine_bg001_rmca.py`、`tools/manufacture_bg001_rmca.py`、`tools/compact_bg001_rmca.py` |
| 端口与浇注盒 | `s1-4_sacrificial_box_and_ports.py`、`config/sacrificial_box_BG001.yaml` |
| SWC 表面重建与共享实现 | `s2_swc_stl_model_generate.py`、`vascular_processing/`、`utils/` |
| 第三方实现及许可 | `third_party/vascularmd/`、`Ultraliser/` |
| 真实数据、运行日志与 QC | `outputs/topbrain_brava_transfer/`、`outputs/*vmd*/`、`outputs/topbrain_*/` |
| 输入与备份 | `vessel_model/T - Brava/`、小鼠单样本原始／AIC 派生目录 |
| 中文方法说明及历史记录 | `references/`、`docs/research_history/` |
| 同步证据、排除清单、测试 | [`reports/github_print_sync/`](reports/github_print_sync/) |

打印候选的共同目录为：

```text
outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/
├── refined_roi/                         # 精修 SWC、映射与 QC
├── manufacturing_roi/                   # 完整制造 ROI、STL、方向与切片记录
└── compact_manufacturing_roi/
    ├── candidates/{MINI,BALANCED,RICH}/  # 3 个紧凑候选与 VascularMD 表面
    ├── BG001_RMCA_BALANCED_print_candidate.stl
    ├── print_fixture_design/            # 原始严格间距版本
    ├── print_fixture_design_adaptive/   # 间距回退版本
    └── print_fixture_design_aligned/    # 最新端面对齐版本；光顺性仍待修复
        ├── core/BG001_RMCA_BALANCED_core_with_ports.stl
        ├── box/BG001_RMCA_BALANCED_casting_box.stl
        ├── box/BG001_RMCA_BALANCED_casting_box.step
        ├── assembly/BG001_RMCA_BALANCED_assembly.3mf
        ├── assembly/assembly_preview.vtm
        ├── QC/
        ├── tables/
        ├── geometry_qc.json
        ├── exported_attachment_alignment.json
        ├── validation_results.json
        └── design_report.md
```

保留历史报告的原始结论；判断本轮状态时以本文开头和 [同步验证说明](reports/github_print_sync/SYNC_REPORT_ZH.md) 为准。`prior_union_printer_probe.json` 仅是先前打印机发现记录，不代表 union 或切片已经完成。

## 获取与检查

```bash
git clone --single-branch --branch sync/vascular-print-models-20260925 \
  https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git ulm_3D_vascular
cd ulm_3D_vascular
python3 tools/verify_print_snapshot.py
```

本分支直接保存实际文件，不依赖 Git LFS。单文件上限 50 MiB；完整原始影像数据集、虚拟环境、嵌套 Git 历史和编译缓存不纳入快照。[排除清单](reports/github_print_sync/excluded_files.csv) 记录了文件或目录及原因。源文件逐项 SHA-256 与同步来源见 [source_file_manifest.json](reports/github_print_sync/source_file_manifest.json)。

## 环境与运行条件

本次验证使用 Python 3.13.11，项目原解释器为 `/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python`。打印／语义模块的实际依赖版本记录在 `requirements-print-snapshot.txt`；它是独立的验证环境记录，不应与旧 `requirements.txt` 的 PyVista 上界混合安装。完整 Ultraliser／HDF5 流程还有自己的依赖，详见 `docs/ULTRALISER_PIPELINE.md` 和环境 JSON 中的未安装条目。本次不宣称已重新编译或验证全部旧后端。

VascularMD 代码和本地兼容补丁已纳入分支；其来源检查会执行 `git rev-parse HEAD`，因此运行需要该目录自己的 Git 元数据。新克隆可以按以下命令恢复固定提交的元数据，保留已随本分支提供的兼容源码；**不要执行 `reset --hard`，不要重复应用补丁**：

```bash
git -C third_party/vascularmd init
git -C third_party/vascularmd fetch --depth=1 https://github.com/megdec/vascularmd.git \
  770feb8fbfb591d6d43f966db57b1465010b9a13
git -C third_party/vascularmd reset --mixed FETCH_HEAD
python3 tools/verify_print_snapshot.py
```

历史清单及保护记录保留 `/home/lzy/projects/ulm_3D_vascular/...` 等原始绝对路径和哈希。本次异地验证实际发现：部分缓存审阅／制造测试在不同项目路径下会触发来源绑定失败。因此这些历史缓存入口应在原项目路径使用；跨机器迁移需单独完成路径迁移及哈希链校验，不能把清单批量替换后仍宣称历史验收原样有效。只读文件校验和 STL／STEP／3MF 读取不依赖该原路径。

在原项目路径查看最新 ROI（沿用现有 UI）：

```bash
python s1-3_swc_roi_generate_MeVO.py --source compact-brava
```

使用当前代码重新生成端口与盒体时，应选择新的输出目录，保留本快照的所有模型和日志：

```bash
python s1-4_sacrificial_box_and_ports.py \
  --config config/sacrificial_box_BG001.yaml \
  --output-root outputs/print_fixture_rerun
```

该命令重放当前算法，不构成尚未完成的 O3 光顺性修复。Windows Bambu Studio、字体及 VMTK 路径见配置；外部可执行程序没有打包。

## 数据来源

- TopBrain Challenge Data Release：University Hospital of Zurich, Department of Neurology；Zenodo record **16878417**；DOI **10.5281/zenodo.16878417**；[官方来源](https://zenodo.org/records/16878417)。原压缩包 MD5：`b703ea31cd1f0e7115a5d3e6e61f59b3`。
- `data/topbrain_provenance/` 保留原始 `License.txt`、下载清单、ITK-SNAP 标签表及 25 份 MRA 标签；不包含完整影像，不能视作完整训练数据目录。完整数据可由 `python tools/topbrain_prepare_dataset.py --download` 获取。
- TopBrain 许可允许非商业使用并要求注明来源；商业用途须取得数据所有者许可。以所附原始许可为准。
- BraVa 原始 SWC、ColorCoded 数据、当前 BG001 派生 SWC、相关论文及小鼠原始备份均按现有项目内容保存。派生语义标签不等同于人工确认的 BraVa 解剖标注。
- VascularMD 的 GPL-3.0 许可与 Ultraliser 原许可随各自源码保留。第三方来源提交记录在同步清单中。
