# 混合方案网格划分检查

`hybrid_velocity_quality_viewer.py` 专门显示混合速度方案中的两套网格，以及它们相对于连续血管壁的划分情况。它不再混入血管编号、黏度、流向等与网格质量无关的场。

## 可查看的规则网格内容

- **Cell type**：规则方格与连续管腔边界的关系；
- **Velocity region**：该方格中心属于有限元区、过渡区还是规则网格区；
- **FEM blend weight**：有限元速度在混合速度中所占的比例；
- **Lumen coverage**：方格被真实管腔覆盖的面积比例；
- **Wall distance**：方格中心到连续固体壁面的距离。

规则网格单元类型为：

| 数值 | 单元类型 |
|---:|---|
| 0 | 血管外，方格没有被管腔覆盖 |
| 1 | 边界切割单元，方格中心在血管外 |
| 2 | 边界切割单元，方格中心在管腔内 |
| 3 | 方格完全位于管腔内 |

这里的“切割”依据连续几何计算出的 `lumen_fraction`，不是把弯曲壁面强行变成台阶。

## 可查看的 DOLFINx 三角网格内容

- **Triangle region**：三角形属于有限元区、跨越过渡区或规则网格区；
- **Relative size**：三角形相对于规则网格步长更细、相当或更粗；
- **Triangle area**：三角形面积；
- **Minimum/maximum edge**：三角形最短边和最长边；
- **Shape quality**：三角形是否过于狭长；
- **Centroid-wall distance**：三角形重心到连续壁面的距离。

三角形区域使用三个顶点到连续壁面的距离来判断：

- 三个顶点都不超过有限元区距离：有限元区；
- 三个顶点都不小于规则网格区距离：规则网格区；
- 其余情况：三角形跨越过渡区。

三角形形状质量为

\[
q=\frac{4\sqrt{3}A}{a^2+b^2+c^2},
\]

其中 \(A\) 是面积，\(a,b,c\) 是三条边。等边三角形有 \(q=1\)，三角形越扁长，\(q\) 越接近 0。

## 使用

打开最新结果：

```powershell
python -m ulm_microbubble_traj_gen.figure_draw.quality_check.hybrid_velocity_quality_viewer
```

指定结果目录：

```powershell
python -m ulm_microbubble_traj_gen.figure_draw.quality_check.hybrid_velocity_quality_viewer `
  --field ulm_microbubble_traj_gen/results/20260723_082654
```

直接保存某一种网格划分图：

```powershell
python -m ulm_microbubble_traj_gen.figure_draw.quality_check.hybrid_velocity_quality_viewer `
  --field ulm_microbubble_traj_gen/results/20260723_082654 `
  --layer fem_triangle_shape_quality `
  --snapshot ulm_microbubble_traj_gen/figure_draw/quality_check/triangle_quality.png `
  --no-show
```

操作：

- 鼠标滚轮：缩放；
- `M`：显示或隐藏 DOLFINx 三角形边；
- `G`：显示或隐藏规则方格边；
- `W`：显示或隐藏连续壁面；
- 鼠标悬停：查看规则单元或最近三角形的分类和质量数据。

查看器只接受包含混合方案网格、区域和连续壁面数据的新 NPZ 文件。
