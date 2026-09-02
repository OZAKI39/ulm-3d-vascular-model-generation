# 分子键合全局／局部动画

本流程使用真实的有限尺寸微泡轨迹与确定性平均场分子键模型。为了尽快看到键合，
靶区不采用原有疾病区域，而是放置在基线轨迹中第一个进入捕获壳的微泡—壁面位置。

## 1. 生成快速观测靶区

```powershell
python ulm_microbubble_traj_gen\figure_draw\prepare_fast_binding_target.py `
  --baseline-result ulm_microbubble_traj_gen\results\20260723_082654
```

默认选择规则和本次结果：

- 捕获距离：`0.50 µm`
- 首个可捕获观测：`t = 0.256 s`
- 微泡：`ID 4`
- 当时壁面间隙：`0.456 µm`
- 靶阳性连续壁面长度：约 `240 µm`

## 2. 生成含键合状态的轨迹

```powershell
python ulm_microbubble_traj_gen\generate_microbubble_trajectories.py `
  --config ulm_microbubble_traj_gen\configs\molecular_binding_visualization.yaml `
  --reuse-field-from ulm_microbubble_traj_gen\results\20260723_082654 `
  --skip-render
```

有效演示参数位于
`configs/molecular_binding_visualization.yaml`：

| 参数 | 数值 |
|---|---:|
| 靶分子面密度 | `270 molecules/µm²` |
| 配体面密度 | `500 molecules/µm²` |
| 有效捕获距离 | `0.50 µm` |
| 静息键长 | `0.45 µm` |
| 二维结合率 | `3.7037037037037035e-14 m²/(molecule·s)` |
| 零力解离率 | `0.06 s⁻¹` |
| 有效键刚度 | `100 pN/µm` |
| 反应柔顺长度 | `0.039 nm` |

这些数值用于稳定、快速地展示键合与滚动，不是文献参数拟合。尤其 `0.50 µm`
捕获距离相对微泡半径较大，程序会给出局部球冠近似适用性的提示。

## 3. 输出 2×2 三模式 GIF

```powershell
python ulm_microbubble_traj_gen\figure_draw\draw_typical_trajectory_modes_2x2_gif.py
```

脚本自动选择真实的代表性轨迹，并在同一时间轴上同步显示：

- 左上：完整血管网、三个移动放大框和通往局部面板的延伸线；
- 右上：正常流动分叉选择；
- 左下：非靶向贴壁运动；
- 右下：靶向分子键合。

默认的正常分叉和非靶向贴壁场景读取带完整旋转诊断的
`results/20260721_150913`，键合场景读取本流程生成的
`figure_draw/outputs/molecular_binding_demo_data`。

输出：

- `figure_draw/outputs/molecular_binding_global_local.gif`
- `figure_draw/outputs/molecular_binding_global_local_peak.png`

三个局部面板使用相同屏幕尺寸、相同圆形外观和相同旋转标记的微泡。管腔和外部
实体保持同一背景颜色；绘图仅使用连续壁面线段，不绘制栅格 lumen mask 或
wall mask。默认使用 `180` 个时间点和 `10 fps`，播放时间约 `18 s`。
