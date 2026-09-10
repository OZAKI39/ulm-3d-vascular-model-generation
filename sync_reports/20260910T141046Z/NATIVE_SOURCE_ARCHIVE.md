# 原生源码归档与重建边界

Mirheo upstream commit 为 `8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf`，原安装源来自 `https://github.com/cselab/Mirheo.git`。HemoCell commit 为 `5a410848bd5c57d5ae1c171112e78eab4a82e650`，来源与 Palabos archive/官方补丁身份见[既有 sources.json](../../hemocell_starter/metadata/sources.json)。两边本轮均没有新编译。

Mirheo 完整第三方树不上传。归档包含[原 CPU 审计引用的 11 个原始源码文件](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/original_source)、[隔离候选的 5 个文本文件](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source)、[新增完整 repair_trace.h](../../mirheo_starter/py_scripts/single_rbc_repair/repair_trace.h)、[许可证](../../mirheo_starter/py_scripts/single_rbc_repair/MIRHEO_LICENSE)。`native/source` 是明确的子集，不能直接作为完整编译树。

复原顺序是检出上述 Mirheo upstream commit，应用[既有 CUDA include 兼容补丁](../../mirheo_starter/metadata/compat_patches/20260908_141405_699555/cuda_include_fix.patch)，再应用[local_before_halo_v2.patch](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/local_before_halo_v2.patch)。第二份补丁包含新头文件，不重复包含兼容补丁。三个修复产物与兼容头可用[冻结计划](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/frozen_benchmark_plan.json)和[清单](SYNC_MANIFEST.json)核对。归档校验只在同步目录的私有临时副本检查补丁应用与字节身份，没有编译或运行。

[新 HemoCell 材料案例](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/hemocell_case/benchmark.cpp)、[CMake](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/hemocell_case/CMakeLists.txt)、[许可证](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/hemocell_case/COPYING)及[与旧案例的补丁](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/hemocell_material_probe.patch)均保留。原生库和系统依赖不上传，既有库的 SHA-256 见[环境快照](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/environment.json)。没有任何编译成功或运行成功的新声明。

原始 JSON 中绝对路径保持不变；[PATH_MAPPING.json](PATH_MAPPING.json)列出精确映射。源码证据里不带目录前缀的 `src/...:line` 指向当时的原安装源，可在 `native/original_source` 核对；被修改的同名文件在 `native/source` 单独保留，二者不能混作一个版本。
