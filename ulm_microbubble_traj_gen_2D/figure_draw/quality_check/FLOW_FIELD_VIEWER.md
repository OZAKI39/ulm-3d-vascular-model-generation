# 流场结果交互检查器

`flow_field_quality_viewer.py` 读取主工作流 `runner.py` 保存的
`velocity_and_wall_shear_field.npz`，用于检查求解完成后的流场和数值诊断结果。

## 启动

从仓库根目录打开最新结果：

```powershell
python -m ulm_microbubble_traj_gen.figure_draw.quality_check.flow_field_quality_viewer
```

指定某次结果：

```powershell
python -m ulm_microbubble_traj_gen.figure_draw.quality_check.flow_field_quality_viewer `
  --field "ulm_microbubble_traj_gen/results/20260728_010604"
```

`--field` 可以指向结果目录，也可以直接指向
`velocity_and_wall_shear_field.npz`。

指定初始图层并保存图片：

```powershell
python -m ulm_microbubble_traj_gen.figure_draw.quality_check.flow_field_quality_viewer `
  --field "ulm_microbubble_traj_gen/results/20260728_010604" `
  --layer divergence_s_inv `
  --snapshot "ulm_microbubble_traj_gen/figure_draw/quality_check/divergence.png" `
  --no-show
```

加入 `--full-color-range` 可使用严格的最小值和最大值作为颜色范围；默认使用
0.5%–99.5% 稳健百分位范围，避免少量极端值掩盖主体结构。

## 图层

查看器根据 NPZ 中实际存在的字段显示图层：

- 最终速度大小、X/Z 速度分量；
- 压力和壁面剪切应力；
- 散度、穿壁速度、最终速度与初始速度之差；
- 初始速度大小及其 X/Z 分量；
- 边界规定速度和开放边界通量。

旧结果缺少某个可选字段时，对应图层不会出现，但其余图层仍可使用。

## 操作

- 右侧单选按钮：切换流场图层；
- 鼠标滚轮：以鼠标位置为中心缩放；
- Matplotlib 工具栏：平移和恢复视野；
- 鼠标悬停：显示单元坐标、速度、压力、壁面剪切应力和散度；
- `V`：显示或隐藏稀疏速度方向箭头；
- `G`：显示或隐藏笛卡尔单元边；
- `W`：显示或隐藏连续管壁；
- `O`：显示或隐藏解剖入口/出口；
- `C`：切换稳健/完整颜色范围；
- `R` 或 `Home`：恢复完整视野。

解剖开口中，蓝色表示入口，橙色表示出口。速度箭头只表达局部方向，
长度经过统一归一化，不用于比较速度大小；速度大小应读取底图颜色或悬停数值。
