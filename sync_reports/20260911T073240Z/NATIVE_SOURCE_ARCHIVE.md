# 原生源码与构建身份

Mirheo upstream 为 `https://github.com/cselab/Mirheo.git` / `8fa67b9aaa7f04c9de2d74a335c9c8c4665068cf`。完整第三方树未重复归档；[当前身份](NATIVE_SOURCE_IDENTITY.json)、[全部 tracked 差异](native_current_tracked_changes.patch)、[实际新头文件](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source/src/mirheo/core/bouncers/repair_trace.h)、[原候选完整补丁](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/local_before_halo_v2.patch)及[许可](../../mirheo_starter/py_scripts/single_rbc_repair/MIRHEO_LICENSE)均保留。tracked 差异不包含新头文件，因此单独保存该原件；不要把两种修复补丁重复应用。native/source 是明确的归档子集，不是可直接构建的完整源树。

[原始 CPU 审计代码](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/original_source)从旧分支继承；[实际修改后的 native/source](../../mirheo_starter/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source)单独保留。[本地已有 CMakeCache](native_build_evidence/CMakeCache.txt)、[隔离构建账本](../../mirheo_starter/runs/single_rbc_repair/rbc_repair_20260910T131105Z/build/budget_ledger.json)、[云端构建身份](../../cloud_results/mirheo-sm120-20731713865ae510/compile_identity.json)和[云端兼容补丁](../../cloud_results/mirheo-sm120-20731713865ae510/cloud_build_compatibility.patch)记录实际环境，没有重新编译。

云端实际库 SHA-256 为 `d45b4fd1b4498b6f365a012e758eec35bc58a7844dabc66f22e617a7b8cf84ee`，原修复树摘要 `2b30d14be357b1bcd7a63805cc164f4907d291b38fd8ae4d4db44980b3f0ac90`。sm_120 的编译/配置兼容变化与物理修复分开记录，库本体不上传。

HemoCell / Palabos 的[源版本与许可](../../hemocell_starter/metadata/sources.json)及自定义案例继承并校对；本次没有运行它们。[路径映射](PATH_MAPPING.json)说明原始 WSL、云端上传副本和仓库位置。
