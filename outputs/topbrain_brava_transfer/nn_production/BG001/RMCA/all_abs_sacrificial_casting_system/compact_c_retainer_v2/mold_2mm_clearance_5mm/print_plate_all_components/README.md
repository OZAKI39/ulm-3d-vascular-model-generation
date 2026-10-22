# 六组件同板打印

主文件：[BG001_RMCA_ALL_ABS_X1C_02mm_single_plate.3mf](BG001_RMCA_ALL_ABS_X1C_02mm_single_plate.3mf)

在 Bambu Studio 中选择“文件 → 打开项目”，载入上面的 3MF 并保留项目设置。项目有一块打印板和六个独立对象，切片结果已包含在文件中。请保留“按层打印”和逐对象支撑配置；组件打印后再装配。

| 对象 | 位置及打印姿态 | 自动支撑 |
|---|---|---|
| 01_ABS_mold_2mm | 前部；底板朝下、顶部开口朝上；侧壁与底板均 2 mm | 关闭 |
| 02_ABS_vascular_core | 后部；沿用已检查的血管芯打印姿态 | 树状支撑 |
| 03_C_retainer_I1 | 右侧第 1 个；沿用原插片轴向端面落床姿态 | 关闭 |
| 04_C_retainer_O1 | 右侧第 2 个 | 关闭 |
| 05_C_retainer_O2 | 右侧第 3 个 | 关闭 |
| 06_C_retainer_O3 | 右侧第 4 个 | 关闭 |

X1 Carbon，0.2 mm 喷嘴，ABS，Textured PEI Plate，0.10 mm 层高，首层 0.12 mm。布局仅平移组件，没有缩放或改变其打印姿态；装配后的血管上下净距仍为 5 mm。各组件保留 10 mm brim 范围，血管另预留 6 mm 支撑范围，最小分隔约 8.57 mm。

实际 Bambu Studio 02.07.01.57 切片成功，估计耗时 12 小时 48 分钟、ABS 32.93 g。模型、支撑均在可用打印板内，六个组件的模型挤出路径与各自参考切片在 0.2 mm 容差下的采样覆盖率均为 100%。该检查是软件几何与切片核验，尚非实物打印或材料参数校准结果。

- [布局预览](single_plate_layout.png)
- [总检查记录](single_plate_summary.json)、[布局检查](layout_audit.json)、[实际刀路检查](actual_toolpath_audit.json)
- [Bambu 原生预览](bambu_all_abs_plate/native_preview.png)
- [实际切片命令](bambu_all_abs_plate/command.json)、[输出日志](bambu_all_abs_plate/stdout.log)、[错误日志](bambu_all_abs_plate/stderr.log)
- [回归测试日志](pytest.log)

`input_project.3mf` 是切片前的中间工程；`plate_layout_review.stl` 是审阅网格。直接使用本页顶部的原生 3MF，即可保留独立对象、打印预设和切片。

重新生成：

```bash
cd /home/lzy/projects/vascular_printing
./run_python.sh s1-10_all_abs_sacrificial_casting_system.py --print-plate
```

源码：`vascular_processing/all_abs_print_plate.py`。可用 `--output <新目录>` 另存。当前模型来自父目录 `mold_2mm_clearance_5mm/`；历史模型和源文件保持不变。程序只在本地生成和切片文件，不上传或启动打印机。
