# 当前新流场：应变率与等效剪切率旋转可视化

[旋转动画](animations/shear_rate_full_vessel_enlarged.mp4) · [预览入口](OPEN_RESULTS.html) · [4K 图](figures/shear_rate_overview_4k.png) · [PDF](figures/shear_rate_overview_4k.pdf)

本次只读取当前 ROI-only-balanced-pressure-v1 最终第 71 步的真实 FEM 速度场，计算运动学派生量和绘图。没有启动 CFD，没有改变黏度、边界条件、网格、轨迹或已有 WSS。

## 定义

令速度梯度为 $G_{ij}=\partial u_i/\partial x_j$。应变率张量是

$$\mathbf D=\frac{\mathbf G+\mathbf G^T}{2}.$$

因此用户截图中的 $\mathbf G+\mathbf G^T$ 是两倍的应变率张量。张量没有唯一的“大小”命名，本次分别保存并明确区分：

$$\|\mathbf D\|_F=\sqrt{\mathbf D:\mathbf D},\qquad
\dot\gamma=\sqrt{2\mathbf D:\mathbf D}=\sqrt2\|\mathbf D\|_F.$$

冒号表示把对应位置的分量相乘后求和：$\mathbf D:\mathbf D=\sum_{i,j}D_{ij}^2$。动画选择后一个量，即等效剪切率 **Shear rate magnitude**，单位 **s⁻¹**。这一定义在简单剪切 $u_x=\kappa y$ 中给出 $|\kappa|$，刚体旋转给出 0。

[OpenFOAM 官方实现](https://cpp.openfoam.org/v13/strainRateViscosityModel_8C_source.html) 使用 `sqrt(2.0)*mag(symm(grad(U)))`，可作为定义核对；本项目实际求解仍是 svMultiPhysics，不是 OpenFOAM。

计算用原始米制坐标和 m/s 速度；不是对速度模长求导，不通过 WSS/黏度生成，也没有使用涡量代替变形率。梯度复用当前生产 WSS 的 `p1_gradients()`：线性四面体每个单元内为常量。

当前 CFD 仍为恒黏度牛顿流体，μ=0.00345312 Pa·s。剪切率也适用于牛顿流场的描述；绘制剪切率不意味着运行了非牛顿模型。这一标量包含剪切及拉伸变形，不能理解为只保留某个坐标方向的剪切分量。

采用完整 D，没有减去 trace(D)/3。当前离散速度的局部散度不严格为零，保留其数值贡献而没有人为修正；`VelocityDivergence_s_inv` 一并保存。WSS/μ 是壁面切向牵引对应的剪切率，在一般三维离散场中不能直接冒充本次等效剪切率。

## 显示方法与风格

- 先在全部 **371,402 个四面体**上计算，再将每个外表面三角形映射到其唯一相邻四面体，直接显示那个单元的剪切率。
- 动画展示 **45,704 个外表面三角形**，其中血管壁 45,221 个，其余为入口/出口人工截面。外表面显示不会呈现内部全部单元；完整内部场保存在 VTU 中。
- 与上一版原始 WSS 面片动画一致，直接使用 cell data；没有节点平均、标量插值、阈值筛选、对数变换或数值裁剪。入口/出口截面的值来自流体单元，不是边界指定的剪切率。
- 0–16,000 s⁻¹ 线性色标覆盖体网格和表面的全部值。沿用配色、黑色背景、英文图内文字、原放大尺度、固定 Z 轴旋转、6° 仰角、出口箭头及渐变标签；只更新物理量名称、单位和数值刻度。
- 保留原有表面法向光照与材质。这只影响光照，不改变单元标量或网格。像素受光照、色表和视频压缩影响，精确值应读取科学数据。
- 1920×1080、24 fps、432 帧、18 秒；18 秒只是稳态场的显示旋转时间。另存四个相差 90° 的 4K PNG，主视角同时提供 PDF。

## 实际数据统计

以下全部为 s⁻¹；分位数使用各区域的面积/体积权重。由于权重和区域不同，不能把全域体积均值当作壁面均值。

| 区域 | 最小 | 加权均值 | P5 | P50 | P95 | 最大 |
|---|---:|---:|---:|---:|---:|---:|
| 全部体网格（体积权重） | 11.334405 | 3025.705957 | 291.957251 | 2824.043450 | 7830.665099 | 15529.155286 |
| 全部外表面，含人工截面（面积权重） | 266.499700 | 4459.774721 | 476.915238 | 4039.714294 | 9604.263830 | 14000.387594 |
| 仅血管壁（面积权重） | 266.499700 | 4476.138586 | 476.087997 | 4052.100300 | 9606.027041 | 14000.387594 |

更多端口统计及应变率张量模长统计见 [rate_summary.csv](data/rate_summary.csv)。

## 可复核数据与计算检查

- [完整体网格 VTU](data/strain_shear_rate_volume_si.vtu)：保持原 SI 坐标、速度和压力，新增下列 cell data。
  - `StrainRateTensor_s_inv`：9 个分量，行优先次序 xx, xy, xz, yx, yy, yz, zx, zy, zz。
  - `StrainRateFrobenius_s_inv`：$\sqrt{D:D}$。
  - `ShearRate_s_inv`：$\sqrt{2D:D}$，动画所示物理量。
  - `VelocityDivergence_s_inv`、`TetraVolume_m3`、`OriginalTetraID_zero_based`。
- [外表面 VTP](data/shear_rate_exterior_si.vtp)：两个标量、相邻四面体编号、原边界标签、面片面积和外法向。
- `prepare_shear_rate.py`：从固定数据派生应变率，复用哈希一致的生产 P1 梯度核心；无求解器调用。
- `render_shear_rate.py`：继承现有表面渲染类，显式选择 `ShearRate_s_inv` 单元字段；相机和标签文件与当前 WSS 逐字节一致。
- `verify_shear_rate.py`：通过另一种边向量叉积/形函数梯度表达式核算全部单元，再检查全部表面对应及完整视频解码。
- `COMPUTE_VALIDATION.json`：实际运行了简单剪切 1000 s⁻¹、刚体旋转 0 s⁻¹、不可压缩单轴拉伸 √3×800 s⁻¹ 的解析场检查；同时核对展开公式和坐标旋转不变性。它们验证这次计算定义与实现，不是新增 CFD 算例或网格收敛证明。
- `MEDIA_VALIDATION.json` 和 `LOCAL_VERIFICATION.json`：映射模式、原始数据完整性、相机/箭头、GPU 设备、432 帧解码检查。

源流场 SHA256：`fb3c6ad0815d156b09b289fc347766486f24ddb521b8e51e8bbc4cf1fd2ace12`。
所有保护输入列在 `audit/source_lock.json`，生产梯度核心按原字节保存在 `audit/wss_core.py`。
当前梯度派生量仍受原 CFD 的空间离散及局部不可压缩误差影响，本次可视化不构成精度或网格无关证明。

## 运行位置与复现

本地：`/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1/strain_shear_rate/`

服务器：`/workspace/roi_only_flow_visualization_20260928/strain_shear_rate/`

服务器 NVIDIA EGL / RTX 4090 渲染，CPU libx264 编码。专用 supervisor 只包含绘图任务，autostart=false、autorestart=false；渲染正常退出后已关闭。

本地复核：
```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -B /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/roi_only_balance_v1/strain_shear_rate/verify_shear_rate.py
```

重新派生：`prepare_shear_rate.py --output <新的空数据输出目录>`；默认拒绝覆盖已完成数据。
重新渲染：在现有数据就绪后，服务器通过专用 supervisor 调用 `/root/particle8_2_runs/env/bin/python -B render_shear_rate.py --output <新的空图片输出目录>`，保留当前 MP4。
