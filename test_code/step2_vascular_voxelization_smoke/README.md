# Step 2：真实血管闭合体素化与四端口几何映射

这是独立、仅处理几何的测试代码。实际运行记录位于：

`/home/lzy/projects/compre_output/step2/20260912_225418/`

结果入口是该目录中的 `STEP2_REPORT.md`、`PARAVIEW_REVIEW.md` 和 `FINAL_STATUS.txt`。
当前结果仅自动核验通过，人工 ParaView 检查仍待完成。未运行时间步、细胞、粒子、BC 或 GPU。

## 文件职责

- `prepare_port_contract.py`：验证冻结 meter STL，读取 Step 1 合同及旧 CFD 的 mesh.dx_m；从完整 STL 恢复真实 cap 三角形、比较面积和法向。输出数值不写死在源码中。
- `closedVascularVoxelizer.cpp`：真实 Palabos TriangleSet → DEFscaledMesh → TriangleBoundary3D → VoxelizedDomain3D → MultiScalarField3D。仿照官方默认参数及 0.001 LU inflate，在 X 两端复制之前停止。
- `CMakeLists.txt`：独立编译单一测试程序，只读链接已有 `build/libhemocell.a`。无顶层工程写入。
- `map_ports_and_diagnose.py`：证明物理/格点映射、检查闭合连通性，按真实 cap 的单格相交邻接建立四个外侧小端口。输出标签、逐体素 witness、VTK。
- `verify_artifacts.py`：独立 VTK 文件回读、原始 STL 包含性检查、开始/结束源码及冻结输入审计。
- `write_reports.py`：读取实测结果生成中文报告、人工核验说明、最终字段输出。
- `run_step2.py`：可选的新运行入口。始终创建新时间戳输出目录；每阶段一次，失败停止，无自动加密或重新运行。

## 本次执行与依赖

使用 `/usr/bin/mpicxx`（系统 g++）、`/usr/bin/cmake`、`/usr/bin/make`、`/usr/bin/mpirun -np 1`。
Python 为 `/usr/bin/python3 -B`，使用现有 NumPy、SciPy、PyYAML、VTK；未激活 Conda，未安装依赖。
构建、运行命令及返回码保存在 run 的 `logs/`。编译并行度1。

本次实际各阶段由逐段命令执行；`run_step2.py` 是把相同最终步骤归并的后续复核入口，本轮未额外执行它，以免重复几何运行。
首次 appended VTI 导出曾被本机 VTK reader 拒绝，保留原文件和日志，随后仅重导出为压缩 inline binary；所有 raw 数组哈希未改变。
当前代码默认 inline binary，所有最终 VTI 已成功回读。

## 后续需要重做时

仅在继续授权几何烟雾测试的任务内执行：

```bash
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
/usr/bin/python3 -B /home/lzy/projects/hemocell_starter/test_code/step2_vascular_voxelization_smoke/run_step2.py
```

入口不修改旧 run、Step 1、旧血管工程或既有 HemoCell 文件；已有静态库缺失则停止。
它不会运行 MPI2、细化网格、流体时间步或 Step3，也不会提交 Git。

本次内存峰值约1.88 GiB、几何运行约50秒，仅供准备机器资源，不能当作血流性能基准。
VTI 文件压缩后很小，但展开有约6860万个单元；ParaView 手册解释如何先 Threshold 选择流体。

## 固定算法范围

合同端口保持 ASSUMED 生理身份。dx=旧 CFD 0.2 µm 只是参考 smoke 分辨率，有效 dx 由最长 Y 轴和整数 refDirN=494决定。
cap 不被圆盘/凸包替代；每个新增外侧体素的法向投影必须在实际 cap 三角形并集内，并有单格线段穿过真实 cap 三角形、直接连到原 lumen 的证据。
标签只标新增层，不改变 HemoCell core flag 语义。实际 BC、生产分辨率、MPI2、人类几何评审均不在本轮验证范围。
