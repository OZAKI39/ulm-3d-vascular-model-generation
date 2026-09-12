# 官方 voxelizer 的本地源码行为

本次逐段检查 `/home/lzy/projects/hemocell_starter/helper/voxelizeDomain.cpp` 及 `.h`。未调用其 getFlagMatrixFromSTL 最终接口。

`voxelizeDomain.cpp:76` 指针版本设 extraLayer=0、borderWidth=1、margin=1；以 DBL 精度读 TriangleSet。
构造 DEFscaledMesh(refDirLength, refDir, margin, extraLayer)、TriangleBoundary3D(*defMesh,false)，调用默认 inflate()。
随后构造 VoxelizedDomain3D，并把 voxelFlag::inside=3、innerBorder=4 转为 flagMatrix=1。

同一函数末尾（搜索 `Since the domain is closed`）继续运行两次 CopyFromNeighbor：

- 目标 x=0..1，源偏移 (+1,0,0)。
- 目标 x=nx-2..nx-1，源偏移 (-1,0,0)。

CopyFromNeighbor::process 按 x 从小到大原位复制。这是管道端面假设，无法表达本合同四个不同方向的 cap。
本轮 test-local C++ 在二值 flagMatrix 完成后结束构造流程；没有调用上述复制，也没有创建流体求解 lattice。

反事实诊断只对导出数据模拟全局单 rank 的原位复制，不调用官方 helper、不改正式输出。其四片可能新增流体数为
[0, 18, 90, 97]，合计 205。
该反事实不声称重现所有未来 MPI 分块/envelope 同步行为。
实测 CLOSED 六个包围盒边界面流体数均为 0。OPENED 的 X-max 有 6 个流体体素，全部属于真实 outlet_03（label=4）的一层 cap opening；其余边界面为 0。

`examples/pipeflow/CMakeLists.txt` 将可执行文件写到现有 example 源码目录；其 compile.sh 使用已有 ../../build 和顶层目标，因此没有用它们编译。
本测试独立 CMake 只编译一个新增 C++ 单元，并只读链接已有 `build/libhemocell.a`；归档哈希见 diagnostics/archive_link_provenance.json。
没有修改顶层 CMake、cmake/setup_googletest.cmake 既有 patch、HemoCell core 或 Palabos。

源码证据：

- `/home/lzy/projects/hemocell_starter/helper/voxelizeDomain.cpp`
- `/home/lzy/projects/hemocell_starter/helper/voxelizeDomain.h`
- `/home/lzy/projects/hemocell_starter/palabos/src/offLattice/voxelizer.h:46`（flag 语义）
- `/home/lzy/projects/hemocell_starter/test_code/step2_vascular_voxelization_smoke/closedVascularVoxelizer.cpp`（本轮实现）
