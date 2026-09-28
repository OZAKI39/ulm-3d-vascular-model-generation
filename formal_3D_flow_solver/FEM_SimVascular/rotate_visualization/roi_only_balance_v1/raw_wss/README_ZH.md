# 原始壁面三角形 WSS 旋转动画

[打开预览页面](OPEN_RESULTS.html) · [MP4](animations/wss_raw_full_vessel_enlarged.mp4) · [4K PNG](figures/wss_raw_overview_4k.png) · [PDF](figures/wss_raw_overview_4k.pdf)

本动画显示当前 ROI-only-balanced-pressure-v1 流场的 `WSS_raw_Pa`，每个血管壁三角形直接使用自己的原始 WSS 模长，单位 Pa。入口/出口人工截面不属于这些壁面面片。只读取现有数据，没有重算 WSS、CFD 或轨迹。

## 与既有 WSS 动画的区别

- 既有版本读取节点字段 `WSS_display_Pa`，由周围壁面面片模长做面积加权节点平均。
- 本版本读取单元字段 `WSS_raw_Pa`。VTK mapper 明确使用 `UseCellFieldData`，选择 `WSS_raw_Pa`，关闭标量插值；没有 cell-to-point 转换。不同三角形可以显示不同的原始值及面片间跳变。
- 完整保留既有相机路径、观察范围、固定 Z 轴、6° 仰角、旁置箭头及文字渐变、背景、配色、照明和字体。仅调整标题/脚注以明确原始面片含义。
- 保留原来的法向照明（smooth normals），这只改变光照，不平均 WSS、不移动网格。光照、色表离散和 MP4 压缩会影响屏幕像素，不能从单个像素反推精确应力值；精确数值以 VTP 中的 cell data 为准。
- 色标仍为 0–55 Pa、线性。当前全部原始值都在范围内，无数据裁剪或归一化。
- 1920×1080、24 fps、432 帧、18 秒。18 秒为稳态场的显示旋转时长，不是物理时间演化。

面片更明显、局部颜色跳变更清楚，是取消节点平均后的显示差异；本次渲染没有改变物理量，也不能据此单独判断 WSS 精度或网格无关性。

## 数据和验证

数据源为父目录 `input_data/field_diagnostics/data/wall_wss_si.vtp`。
该文件保存 45,221 个壁面三角形；精确统计、输入哈希及 GPU 后端记录在 `MEDIA_VALIDATION.json`，完整解码和本地核验见 `LOCAL_VERIFICATION.json`。

`render_raw_wss.py` 继承原 `render_surface_fields.SurfaceView`，复用原文字、相机、渲染、视频解码及标签检查函数，仅显式选择原始单元字段。执行时检查 mapper 输入数组逐值相等，面片连接与坐标保持不变，相机与标签文件同原动画逐字节相等。源输入、原渲染代码、原节点平均动画均由 `audit/source_lock.json` 保护。

服务器使用 NVIDIA EGL/RTX 4090 渲染，CPU libx264 编码。执行命令和日志见 `supervisord.conf` 与 `logs/`。该专用任务 autostart=false、autorestart=false，不包含求解任务；渲染正常退出后，专用 supervisor 已关闭。

本地核验命令：
```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -B /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1/raw_wss/verify_raw_wss.py
```

服务器目录：`/workspace/roi_only_flow_visualization_20260928/raw_wss/`。
复现时在服务器通过 supervisor 执行 `/root/particle8_2_runs/env/bin/python -B render_raw_wss.py --output <新的空输出目录>`，保留当前结果。脚本拒绝覆盖已完成 MP4。
