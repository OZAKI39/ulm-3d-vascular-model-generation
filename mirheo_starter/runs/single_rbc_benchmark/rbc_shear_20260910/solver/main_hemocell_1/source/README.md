# 原生单红细胞剪切适配

派生自固定 HemoCell `examples/oneCellShear`、`examples/stretchCell` 和本地纯流体计时小模块，携带 COPYING（AGPL-3.0）。新目标仅链接现有 `build/upstream/libhemocell.a`，不重编译依赖或覆盖示例。

`benchmark.cpp` 真正初始化细胞场、加载一颗 `RbcHighOrderModel`、逐步执行 `hemocell.iterate()`；材料和 IBM 速度更新时间均为 1。X/Y 周期，Z 速度壁，z=0 是 −U，z=(nz−1)dx 是 +U。输出中周期副本 ID 用 HemoCell `base_cell_id` 映射，顶点按唯一所有权收集；丢失或重复顶点使任务失败。

静态探针直接调用原生 `ParticleMechanics`，无代理材料力。`reference.off` 来自实际加载后的旋转与平移，输出使用共同数值单位；不混用官方日志中的 deformation index。

配置 XML **必须有 XML declaration**。固定版本 Config::load 使用 `FirstChild()->NextSibling()` 决定 checkpoint 标志，缺少声明会留下未初始化状态；本案例生成器补全声明，不修改上游库。

由 Mirheo starter 中的 `py_scripts.benchmark_single_rbc --prepare-cpu` 管理独立输出、限制和账本。完整说明见 `/home/lzy/projects/mirheo_starter/test_code/README_single_rbc_benchmark.md`。
