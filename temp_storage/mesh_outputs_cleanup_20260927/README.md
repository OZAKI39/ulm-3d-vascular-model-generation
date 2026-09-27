# 当前网格输出目录清理回执

范围：`formal_3D_flow_solver/FEM_SimVascular/outputs/mesh_and_flow/`。

- 删除 `vascular_flow/`：旧未验收瞬态场，质量守恒失败且未达到稳态。
- 删除 `official_fluid_smoke/`：早期官方示例测试输出。
- 删除 `plot_cache/`、`api_documentation.json`：可重新生成的缓存和接口说明快照。
- 保留 `model/`、`mesh_generation/primary/`、`solver_mesh/`：当前表面导入、网格生成与检查所需数据。

`before.json` 仅记录删除/保留文件的路径、哈希和大小，没有复制被删除的数据。`result.json` 记录实际结果：删除 21 个文件路径（19 个独立文件，包含 2 组硬链接），释放文件数据块 14,909,440 字节，约 14.22 MiB；保留 19 个数据文件，哈希全部不变。

保留求解网格的 9 个文件与当前 H0 算例一致，当前冻结流场哈希保持不变。服务器未修改。
