# 排除项与恢复方式

本次保留正式网格、校准/正式CFD的原始结果与日志、失败记录、两姿态10段流场MP4、4K PNG/PDF、微泡样本与所有已生成的基准/未完成轨迹。既有小鼠1500条、RBC及WSS结果不删除。

未纳入的内容：

- Python/Conda环境、Git数据库、安装包、缓存、套接字/PID、可重建共享库/对象文件/求解器二进制。求解器完整修复源码、补丁、构建命令和二进制SHA256保留。
- BraVa `microbubble/data/gpu_mesh_input.npz`，109,638,745字节。它是完整网格/节点场/梯度/壁面数组的派生输入，超过GitHub单文件限制；正式VTU、NPZ场和壁面数据保留，CUDA计算记录也保留。原副本仍在WSL及服务器。
- 以前已经排除的用户取消细档大网格、第三方示例输出和上游未下载的LFS真实载荷。LFS指针文本记录于`server_inventory.json`，不能视为真实测试数据。
- 与本轮流场无关的大型TopBrain原始影像、全数据集、切片安装包、G-code/3MF/ZIP和重复完整盒体；保留所用BG001 SWC、refined ROI、BALANCED候选、四端口芯、姿态清单、配置与开发源码。其它打印候选及打印专用完整输出仍在原`vascular_printing/`工作区。

逐项排除路径、原因、大小和可用SHA256见 `excluded_files.json`、`server_inventory.json`。目录级运行环境不会逐文件扫描私有配置。它们是排除清单，不表示源文件已删除。

仅恢复现有GPU输入文件（不会启动计算）：

```bash
scp -p vast4090:/workspace/brava_flow_roi_18mlmin_20260928/microbubble/data/gpu_mesh_input.npz brava_flow_roi_18mlmin/microbubble/data/
sha256sum brava_flow_roi_18mlmin/microbubble/data/gpu_mesh_input.npz
```

校验值请与清单对应记录核对。不要通过运行整个`campaign.py prepare`恢复缓存，以免新建或覆盖样本契约。当前微泡取消标记应保留，用户明确授权恢复之前不执行计算命令。

同步源/配置/结果保持原始字节；导航文档变换和Git规则归档另有记录。使用者仍需按原路径契约配置科学环境；当前快照不是一键部署包。
