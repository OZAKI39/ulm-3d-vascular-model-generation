# 小鼠微血管、3D 流场、微泡与 RBC 当前工作区

这是 **2026-09-27 当前开发内容的同步快照**，分支 `sync/current-h0-dt1ms-wss-audit-20260927`。目录对应清理后的 `/home/lzy/projects/` 工作区；旧分支仍保留其历史。本次只整理和同步文件，未调整物理参数、边界条件或科学计算代码。

## 当前入口

| 内容 | 仓库内路径 |
|---|---|
| 血管 A、ROI、1D/0D Network-H0 | [ulm_3D_vascular/](ulm_3D_vascular/README.md) |
| 当前几何的原始 SWC、影像/掩膜及共享预处理 | [vascular_printing/](vascular_printing/)；只收录当前 A 的数据与所需运行结果 |
| 表面导入、四面体网格及检查 | [formal_3D_flow_solver/](formal_3D_flow_solver/README.md) |
| **当前 H0 边界条件的正式 FEM 算例** | [mean-2p0-mmps-A-H0-pressure-v1/](ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/) |
| 当前微泡：生成、积分、500 条轨迹、dt=1 ms | [ulm_particle_formal_p9a5/](ulm_particle_formal_p9a5/README.md)、[结果索引](ulm_particle_formal_p9a5/CURRENT_RESULTS.md) |
| RBC 几何、轨迹及独立展示 | [RBC 结果索引](ulm_particle_formal_p9a5/CURRENT_RESULTS.md) |
| 流线、速度矢量、压力、WSS 可视化 | [OPEN_RESULTS.html](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/OPEN_RESULTS.html)（下载后用浏览器打开） |
| **最新中文 WSS 审计及数值证据** | [WSS_AUDIT_REPORT.md](formal_3D_flow_solver/FEM_SimVascular/wss_audit/WSS_AUDIT_REPORT.md) |
| 新流场残差表格 | [A_Global_Transient_Evolution_A_H0.xlsx](A_Global_Transient_Evolution_A_H0.xlsx) |
| 服务器当前路径 | [CURRENT_SERVER_PATHS.md](CURRENT_SERVER_PATHS.md) |
| 本次同步范围、验证、排除项 | [SYNC_REPORT_ZH.md](sync_metadata/current_20260927/SYNC_REPORT_ZH.md) |

当前冻结流场为 `steady_flow_mean_2p0_mmps_A_H0.vtu`，SHA256：
`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。

500 条正式微泡轨迹结果为 O1=0、O2=139、O3=291、接触支持静止=70，计算失败=0。正式动画及最终交付核验仍有未完成项，不能把较早的 Particle HTML 页面当成全部完成证明。RBC 展示和 HemoCell 留存资料使用各自的参考场/模型，**不表示已完成当前 H0 场中的 RBC 耦合生产计算**。

## 使用与复核

本快照保留原科学代码和数据字节，部分入口仍使用 `/home/lzy/projects/` 绝对路径，Python/MPI/PETSc/SimVascular 运行环境没有打包。因此这不是在任意目录克隆后即可完整运行的安装包。复现全工作流时应按各工作区 README 配置依赖和原目录布局；不要用历史 `frozen_reference` 替代当前 H0 场，也不要批量替换受保护代码中的路径而不更新来源契约。

`wss_audit/` 可以单独下载复核；其 [REPRODUCE.md](formal_3D_flow_solver/FEM_SimVascular/wss_audit/REPRODUCE.md) 记录依赖与命令。核心命令：

```bash
python -B formal_3D_flow_solver/FEM_SimVascular/wss_audit/scripts/verify_bundle.py
python -B ulm_particle_formal_p9a5/scripts/verify_current_data.py
```

其他测试的确切执行命令、日志、SHA256 清单及服务器文件映射保存在 `sync_metadata/current_20260927/`。服务器缺少的当前证据按内容去重补入 `server_evidence/current_20260927/`；原目录结构可由映射表还原，它不是另一份应直接覆盖当前代码的部署包。
