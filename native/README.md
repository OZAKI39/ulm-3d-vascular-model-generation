# Native solver

唯一自定义实现位于 `src/`：

| 文件 | 职责 |
|---|---|
| `bond_bpm_peri_variational.cpp/.h` | 非均匀 endpoint-star 能量/力、原始 ID、逐键阈值、κ/d、损伤耗散、拓扑同步、V5 重启序列化、C API |
| `fix_lps_variational_dynamic.cpp/.h` | Verlet 漂移后一次损伤提交、质量/回调顺序检查；V4 中材料/接触力已是 N，不再乘统一体积 |
| `pair_bpm_peri.cpp/.h` | 非存活参考键对的线性排斥，平均节点体积平方的对称系数 |
| `mapped_geometry.cpp` | 映射胞射线、三角面遮挡与参考键穿面查询；几何核，不是材料求解器 |

前六份文件逐字复制当前有效 V4 的生成结果；没有重新执行历史字符串补丁链。`mapped_geometry.cpp` 也保留原始字节。哈希见 `archive/reports/MIGRATION_ASSETS.json` 和构建后的 `build/BUILD.json`。

`vendor/lammps_core_20260930.tar.gz` 是实际历史原生构建使用的必要核心、BPM、MOLECULE、fmt/nlohmann、CMake 和匹配 Python 包的源码快照，约 4.4 MiB。它不是整个旧科研工程，不包含旧 Python 血栓模型、计算结果或对象文件。自定义六个文件从包中排除，由 `src/` 覆盖进构建暂存。版本头为 `30 Sep 2026`；不能把这个历史修改快照称为未经修改的官方发布版。内容锁定在 `vendor/MANIFEST.json`，原迁移仓库的提交为 `519a5a049285af2bb6b99e95506e3df6a6ebbb3c`，实际构建快照另以内容 SHA 为准。

```bash
source .venv/bin/activate
python native/build.py --jobs 4
source native/env.sh
python -m pytest -q
```

构建使用 C++17、Release、shared library、MPI、OpenMP、BPM 和 MOLECULE。`native/build.py` 安装匹配的正常 `lammps` Python 包，再完整编译共享库和映射几何核；没有历史 `.o` 链接、动态 Python 源码执行或运行时自动编译。`native/work/` 可删除并重建，不能在那里维护自定义模型。

外部依赖：系统 OpenMPI、libstdc++、libgomp、zlib 和 libc；版本/链接信息在验证报告中。本次 WSL 使用系统 OpenMPI，未依赖服务器 `/workspace/y_clot/deps`。

原已接受库 SHA：`4510aed9195a30b03a4895435be05eca005d462507a12c4423fe8d2a854481e2`。本地重新编译产生不同二进制 SHA，使用新的 `BUILD.json` 记录；没有改写旧冻结记录。是否可用以当前回归测试的实际覆盖范围为准，不能用“相同 C++ 文件”替代数值验证。

LAMMPS 依赖按所附 `LICENSE`（GPL-2.0）保留。历史 C++ 注释中出现的 V3/力密度措辞属于来源记录；V4 的实际 SI 约定以代码、Python 校验和测试为准。
