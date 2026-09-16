# Microbubble 与墙模型科研交接 — 2026-09-16

入口：[NEW_CHAT_CONTEXT.md](NEW_CHAT_CONTEXT.md) → [CURRENT_PROJECT_STATUS.md](CURRENT_PROJECT_STATUS.md) → [CURRENT_STATE.json](CURRENT_STATE.json)。这是独立归档，不是生产模型升级。

本目录补齐上次 HemoCell Stage4/RBC Stage1 交接后 9 个阶段的项目代码、冻结合同、输入、报告、验证表、图及适合 Git 保存的原始数据。所有复制文件保持原字节；SOURCE_MANIFEST.tsv 逐文件记录来源与 SHA256。既有仓库内容未改。未创建 PR 或合并 main。

## 当前结论

当前墙 V0 被真实几何的局部平面有效性门槛阻塞；最新固定硬球墙审查为 FAIL_WALL_REPRESENTATION。最优已测诊断配置仍有 88.678% 核心最坏矩阵作用误差，超过 5% 门槛。不能将失败写成 production-ready，也不能外推为所有 regularized blobs 不可行。下一阶段仅建议 BEM_CURVED_WALL_FEASIBILITY_AUDIT，未启动。

## 数据与完整性

在本目录执行：

```bash
sha256sum -c SHA256SUMS
```

这个顶层清单覆盖本次 Git 交接内容。各阶段子目录原有 SHA256SUMS、FROZEN_INPUT_SHA256SUMS 及 provenance 清单是完整原始归档的历史证据，可能引用本次明确未上传的文件；不能把它们误认为这个裁剪副本的清单。9 个完整原归档的清单在复制前均实际通过，记录在 provenance/ORIGINAL_ARCHIVE_SHA_VERIFICATION.json。

最新 fixed_multiblob/remote_raw 保持完整，可另执行：

```bash
cd fixed_multiblob/remote_raw
sha256sum -c SHA256SUMS
```

固定小球审查的全部远端原始 HDF5、有效/失败查询、独立分析脚本、输入和全部 12 张 PNG/4 个 VTP 均在仓库内。未运行的曲壁/真实 STL 图及空 VTP 有明确 NOT_RUN 标签，不能当作运行结果。主派生 FIXED_MULTIBLOB_WALL_AUDIT.h5（569742131 bytes）留在原本地归档。

OMITTED_ARTIFACTS.tsv 列出所有排除项：二进制、构建缓存、临时服务器调度信息、外部文献完整副本、大型或重复轨迹与主派生 HDF5；每项保留原路径、大小、SHA256 和保留位置。原本地文件没有删除或修改。Git 数据不是完整运行环境；没有捆绑编译器或 Python venv。

## 重算最新静态审查

在单独工作副本进行，避免修改已冻结的 Git 交接目录。用 fixed_multiblob/source/requirements-lock.txt 建立匹配环境后，运行：

```bash
cd fixed_multiblob
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python scripts/finalize_fixed_multiblob_wall_audit.py
python scripts/visualize_fixed_multiblob_wall.py
```

finalizer 从已保存矩阵重算，不启动水动力求解；会写出约 543 MiB 的派生 HDF5 并更新表格。HDF5 容器元数据可能使文件级 SHA 不同，应比较数值、合同 SHA、查询数与分类。历史报告中的原机路径是 provenance，不是新机器安装路径；本次未重跑任何模拟。

## 许可证和使用边界

项目自定义代码保留来源记录。Pecnut 为 MIT，RMBW 与 RigidBodyIB 的许可证及固定提交记录随相关阶段归档；其原有许可继续适用，本目录不改变上游许可。外部论文使用报告中的引用与 URL，完整论文副本未上传。服务器地址、私钥和凭据不属于公开交付。
