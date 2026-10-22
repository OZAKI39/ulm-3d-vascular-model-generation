# 血管打印模型工作区

> **当前通道交付（2026-10-15）：** Y 形、三端口直管和双上端口直管均采用 ABS 血管芯、PLA 夹具，交付未切片工程。双上端口芯另有约 45° 少支撑及水平离板约 1.20 mm 的新姿态，见 [最新姿态说明](docs/Straight_2Top少支撑与水平离板姿态.md) 和 [开发与使用指南](开发与使用指南.md)。下文保留 BraVa 与旧制造阶段记录。2026-10-22 GitHub 快照的保留范围、排除清单及环境复现方法见 [同步说明](docs/GitHub同步说明_20261022.md)。

**新服务器入口（2026-10-09）：** `95.3.33.46:15283` 的 `/workspace/vascular_printing` 使用独立 `.venv` 和 Linux Bambu Studio。连接、迁移范围、验证记录及继续开发命令见 [新服务器使用说明](docs/新服务器迁移与使用_20261009.md)。

入口 `./run_python.sh s1-11_vascular_core_support.py --printer p1s-0.4` 当前针对 **Bambu P1S / 0.4 mm / ABS**，只重新生成冻结血管芯及显式支撑，不生成外模或 C 插片。默认打印机已改为 P1S；独立输出在 `vascular_core_support_production/P1S_04_ABS_20261009/`。使用 `--stage review` 复核已有交付，不重新切片。**真实切片成功与制造验收是两项不同检查，状态以 `regeneration_summary.json` 为准；`NEEDS_ADJUSTMENT` 文件仅供诊断。** 模型、3MF、检查图和日志路径见 [P1S 0.4 mm 使用说明](docs/S1_11_P1S_04_ABS.md)。历史 X1 Carbon 输出 `final_practical_v3/` 原位保留，复核须显式使用 `--printer x1c-0.2`，见 [历史实用验收](docs/S1_11最终实用验收.md)。

**六组件一起打印：** [打开最新同板 3MF](outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/all_abs_sacrificial_casting_system/compact_c_retainer_v2/mold_2mm_clearance_5mm/print_plate_all_components/BG001_RMCA_ALL_ABS_X1C_02mm_single_plate.3mf)，在 Bambu Studio 中选择“打开项目”，保留 X1 Carbon / 0.2 mm / ABS 与逐对象支撑设置。2 mm 外模、血管芯及四个 C 插片已独立排在同一打印板上并通过实际切片。复现命令：`./run_python.sh s1-10_all_abs_sacrificial_casting_system.py --print-plate`。布局与检查记录见 [全 ABS 使用说明](docs/ALL_ABS牺牲浇注系统使用说明.md#六组件同板打印)。

当前入口仍为 `./run_python.sh s1-10_all_abs_sacrificial_casting_system.py`。已将全 ABS 盒体侧壁和底板统一为 **2 mm**，完整血管到内底面、顶部开口的净距均为 **5 mm**；继续使用 90° 紧凑 C 插片，血管和插片几何不变。当前输出在 `R/all_abs_sacrificial_casting_system/compact_c_retainer_v2/mold_2mm_clearance_5mm/`，`R=outputs/topbrain_brava_transfer/nn_production/BG001/RMCA`。模型、装配坐标、实际 X1 Carbon / 0.2 mm / ABS 打印文件及日志见 [全 ABS 系统使用说明](docs/ALL_ABS牺牲浇注系统使用说明.md)。旧模型原位保留。

**s1-9 PLA 夹具（含 V2）已标记为 RETIRED / RESEARCH HISTORY**。原代码、最终报告、数据及源哈希保持原样；[BCF V2 说明](docs/BCF_V2内墙支撑与切割使用说明.md) 仅用于历史研究。s1-8 及下文 candidate 0/15、一体盒工艺也属于历史路线。

这里集中存放血管数据处理、MeVO 区域提取、血管实体建模、四端口衔接、ABS 倒模盒生成及 Bambu 切片的代码、数据和日志。上游路径是 **TopBrain → BraVa BG001 RMCA → 血管 ROI → 四端口**；当前制造阶段使用 **独立冻结 ABS 芯 + 全 ABS 牺牲外模**。旧一体式倒模件保留作历史。

## 历史 candidate 0 / 15 模型

| 候选 | 实际姿态 | 文件目录 | 整套压缩包 | 预计打印时间 / 总耗材 |
|---|---|---|---|---|
| Candidate 0 | 底板贴打印板，顶部开口朝上 | [candidate_0](deliverables/candidate_0/) | [candidate_0.zip](deliverables/candidate_0.zip) | 约 111 分钟 / 55.77 g |
| Candidate 15 | +X 侧壁贴打印板，开口朝侧方 | [candidate_15](deliverables/candidate_15/) | [candidate_15.zip](deliverables/candidate_15.zip) | 约 191 分钟 / 70.40 g |

两个目录均有 **STL、开启树状支撑的真实 Bambu 3MF、G-code、姿态预览、变换和检查记录**。在 Bambu Studio 中打开 `_support_ON.3mf` 即可看对应切片；不要把两个候选叠加到同一个打印任务。当前文件对应 P1S、0.4 mm 喷嘴、ABS、0.20 mm Standard 工艺。

血管、四个端口和五面盒体已经合为一个实体；四周和底板封闭，顶部开放。当前结果可供人工审核，尚未完成实物打印、支撑拆除和灌注验证。Candidate 0 的支撑清理通路代理存在风险标记；Candidate 15 通过当前代理检查，但支撑更多。这些检查没有证明实际支撑一定取不出或一定能取出。

## 在哪里看说明

- [开发与使用指南](开发与使用指南.md)：后续开发接手入口，说明整体逻辑、目录职责、各阶段输入输出、运行依赖和修改边界。
- [使用说明](docs/使用说明.md)：查看模型、运行可视化、生成新倒模件、重新导出两种候选。
- [目录与数据说明](docs/目录与数据说明.md)：代码、原始数据、模型、日志分别在哪里。
- [历史一体式倒模评估报告](deliverables/casting_mold_report.md) 与 [检查图](deliverables/QC/)：几何、O3 间距、支撑和真实切片结果。
- [迁移与清理记录](docs/迁移与清理记录.md)：迁移范围、兼容链接、删除清单和验证结果。

## 运行环境

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-3_swc_roi_generate_MeVO.py --source compact-brava --no-gui
```

`run_python.sh` 自动进入本目录，使用已有共享虚拟环境。`.venv` 是指向 `/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv` 的链接；共享环境已于 2026-09-27 迁入 temp_storage，环境链接及启动器已更新；既有依赖版本未升级。已验证版本记录在 [requirements-printing.txt](requirements-printing.txt)。

迁移保留了历史报告及其原始哈希，因而部分历史记录仍写着旧路径。旧工作区中的兼容链接会把它们指向本目录；新开发、查看和运行均从这里开始。CFD、Network H0 和流场求解仍属于原工作区，没有混入本项目的打印入口。
