# 血管打印模型工作区

这里集中存放血管数据处理、MeVO 区域提取、血管实体建模、四端口衔接、ABS 倒模盒生成及 Bambu 切片的代码、数据和日志。主要路径是 **TopBrain → BraVa BG001 RMCA → BALANCED 血管 → 四端口 → 一体式 ABS 倒模件**。

## 先拿到当前模型

| 候选 | 实际姿态 | 文件目录 | 整套压缩包 | 预计打印时间 / 总耗材 |
|---|---|---|---|---|
| Candidate 0 | 底板贴打印板，顶部开口朝上 | [candidate_0](deliverables/candidate_0/) | [candidate_0.zip](deliverables/candidate_0.zip) | 约 111 分钟 / 55.77 g |
| Candidate 15 | +X 侧壁贴打印板，开口朝侧方 | [candidate_15](deliverables/candidate_15/) | [candidate_15.zip](deliverables/candidate_15.zip) | 约 191 分钟 / 70.40 g |

两个目录均有 **STL、开启树状支撑的真实 Bambu 3MF、G-code、姿态预览、变换和检查记录**。在 Bambu Studio 中打开 `_support_ON.3mf` 即可看对应切片；不要把两个候选叠加到同一个打印任务。当前文件对应 P1S、0.4 mm 喷嘴、ABS、0.20 mm Standard 工艺。

血管、四个端口和五面盒体已经合为一个实体；四周和底板封闭，顶部开放。当前结果可供人工审核，尚未完成实物打印、支撑拆除和灌注验证。Candidate 0 的支撑清理通路代理存在风险标记；Candidate 15 通过当前代理检查，但支撑更多。这些检查没有证明实际支撑一定取不出或一定能取出。

## 在哪里看说明

- [使用说明](docs/使用说明.md)：查看模型、运行可视化、生成新倒模件、重新导出两种候选。
- [目录与数据说明](docs/目录与数据说明.md)：代码、原始数据、模型、日志分别在哪里。
- [当前倒模评估报告](deliverables/casting_mold_report.md) 与 [检查图](deliverables/QC/)：几何、O3 间距、支撑和真实切片结果。
- [迁移与清理记录](docs/迁移与清理记录.md)：迁移范围、兼容链接、删除清单和验证结果。

## 运行环境

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-3_swc_roi_generate_MeVO.py --source compact-brava --no-gui
```

`run_python.sh` 自动进入本目录，使用已有共享虚拟环境。`.venv` 是指向 `/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv` 的链接；共享环境已于 2026-09-27 迁入 temp_storage，环境链接及启动器已更新；既有依赖版本未升级。已验证版本记录在 [requirements-printing.txt](requirements-printing.txt)。

迁移保留了历史报告及其原始哈希，因而部分历史记录仍写着旧路径。旧工作区中的兼容链接会把它们指向本目录；新开发、查看和运行均从这里开始。CFD、Network H0 和流场求解仍属于原工作区，没有混入本项目的打印入口。
