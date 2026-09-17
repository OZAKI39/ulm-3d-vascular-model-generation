
`build_domain_from_vessels()` 的作用可以概括为一句话：

> 根据全部血管的位置，在 $X-Z$ 平面上铺一张足够大的方格纸，供后续血管绘制和血流计算使用。

虽然血管端点用三维坐标 $(x,y,z)$ 表示，但当前血流模拟主要在 $X-Z$ 平面中进行，因此最终生成的是二维网格。

## 1. 输入是什么

函数接收两个输入：

```python
build_domain_from_vessels(vessels, cfg)
```

其中：

- `vessels`：所有血管段。
- `cfg`：网格设置。

每条血管都有两个端点：

$$
\mathbf{x}_p=(x_p,y_p,z_p)
$$

和

$$
\mathbf{x}_d=(x_d,y_d,z_d).
$$

可以把一条血管理解为连接这两个点的一根线段。

配置中主要使用：

- `grid_spacing_um`：相邻网格点之间的距离。
- `padding_um`：血管网络外面额外保留的空白。
- `max_grid_cells`：允许生成的最大网格规模。

所有长度的单位都是微米 $\mu\mathrm{m}$。

## 2. 收集所有血管端点

假设共有 $N$ 条血管，每条血管有两个端点，那么一共有 $2N$ 个点：

$$
\left\{ \mathbf{x}_{p,1},\mathbf{x}_{d,1}, \mathbf{x}_{p,2},\mathbf{x}_{d,2}, \ldots, \mathbf{x}_{p,N},\mathbf{x}_{d,N} \right\}.
$$

代码把这些点放入一个数组：

```python
points = np.asarray(
    [point for vessel in vessels for point in (vessel.x_p, vessel.x_d)],
    dtype=float,
)
```

这样就可以统一寻找整个血管网络的最左、最右、最上和最下位置。

## 3. 找到血管网络的边界

设所有血管端点中的最小和最大 $X$ 坐标分别为：

$$
x_{\min}^{\mathrm{vessel}}, \qquad x_{\max}^{\mathrm{vessel}}.
$$

同样，$Z$ 方向的范围是：

$$
z_{\min}^{\mathrm{vessel}}, \qquad z_{\max}^{\mathrm{vessel}}.
$$

如果网格刚好贴着血管网络边缘，靠近边界的血管可能没有足够的计算空间。因此，代码在四周增加宽度为 $p$ 的空白，其中

$$
p=\texttt{padding\_um}.
$$

扩展后的范围为：

$$
x_{\min}=x_{\min}^{\mathrm{vessel}}-p,
$$

$$
x_{\max}=x_{\max}^{\mathrm{vessel}}+p,
$$

$$
z_{\min}=z_{\min}^{\mathrm{vessel}}-p,
$$

$$
z_{\max}=z_{\max}^{\mathrm{vessel}}+p.
$$

可以把它想成：先用一个长方形框住所有血管，然后把长方形的四条边向外移动一段距离。

需要注意，函数根据血管端点确定边界，并没有读取血管半径。因此 `padding_um` 必须足够大，才能让血管本身及其周围区域完整落入网格。

## 4. 计算需要多少个网格点

设相邻网格点之间的距离为：

$$
h=\texttt{grid\_spacing\_um}.
$$

$X$ 方向的总长度为：

$$
L_x=x_{\max}-x_{\min}.
$$

$Z$ 方向的总长度为：

$$
L_z=z_{\max}-z_{\min}.
$$

代码计算：

$$
n_x= \left\lceil \frac{L_x}{h} \right\rceil+1,
$$

$$
n_z= \left\lceil \frac{L_z}{h} \right\rceil+1.
$$

其中 $\lceil a\rceil$ 表示向上取整。

向上取整是为了保证网格一定能够覆盖完整区域。例如：

$$
\frac{L_x}{h}=10.2
$$

表示需要超过 10 个间隔。若向下取整为 10，网格就会少覆盖一部分，因此必须取：

$$
\left\lceil10.2\right\rceil=11.
$$

公式最后再加 $1$，是因为“间隔数”和“网格点数”不同。例如一条长度为 $10\,\mu\mathrm{m}$ 的线，每隔 $5\,\mu\mathrm{m}$ 放一个点：

$$
0,\quad5,\quad10
$$

共有 2 个间隔，但有 3 个点。

## 5. 生成每个网格点的坐标

$X$ 方向第 $i$ 个网格点的位置为：

$$
x_i=x_{\min}+ih, \qquad i=0,1,\ldots,n_x-1.
$$

代码对应：

```python
x_coordinates = min_xyz[0] + spacing * np.arange(nx)
```

$Z$ 方向第 $j$ 个网格点的位置为：

$$
z_j=z_{\min}+jh, \qquad j=0,1,\ldots,n_z-1.
$$

代码对应：

```python
z_coordinates = min_xyz[2] + spacing * np.arange(nz)
```

因此，每个二维网格位置都可以表示为：

$$
(x_i,z_j).
$$

整个网格共有：

$$
N_{\mathrm{grid}}=n_xn_z
$$

个位置。

## 6. 为什么还要计算一个固定的 $Y$ 坐标

原始血管坐标是三维的，但这里进行的是 $X-Z$ 二维计算。为了让输出仍然能够写成三维坐标，代码取所有血管端点 $Y$ 坐标的平均值：

$$
\bar y=\frac{1}{2N}\sum_{k=1}^{2N}y_k.
$$

代码是：

```python
fixed_y = points[:, 1].mean()
```

之后整个二维网格都放在：

$$
y=\bar y
$$

这个平面上。

因此，网格左下角的三维坐标为：

$$
\mathbf{x}_{\mathrm{origin}}=(x_{\min},\bar y,z_{\min}).
$$

这里的固定 $Y$ 只是为了把二维 $X-Z$ 计算结果放回三维坐标格式，并不表示程序在 $Y$ 方向也进行了网格计算。

## 7. 哪些输入会被直接拒绝

这些输入无法产生可靠网格，因此函数会直接报错。

网格间距必须满足：

$$
0<h<\infty.
$$

如果 $h\leq0$、`NaN` 或无穷大，就无法正常铺设网格。

边缘空白必须满足：

$$
0\leq p<\infty.
$$

负的边缘空白会缩小计算区域，可能切掉血管，因此不允许。

最大网格数量必须满足：

$$
N_{\max}>0.
$$

血管端点的所有坐标也必须是有限数值，不能包含 `NaN` 或无穷大。

此外，生成的网格至少要满足：

$$
n_x>1, \qquad n_z>1.
$$

如果

$$
n_xn_z>N_{\max},
$$

程序也会停止，因为网格过大会占用过多内存和计算时间。

## 8. 最终返回什么

函数最终返回一个 `GridDomain`，里面保存：

- 网格起点：
  $$
  (x_{\min},\bar y,z_{\min})
  $$
- 网格间距：
  $$
  h
  $$
- 网格形状：
  $$
  (n_x,n_z)
  $$
- 固定的 $Y$ 坐标：
  $$
  \bar y
  $$
- 所有 $X$ 坐标。
- 所有 $Z$ 坐标。

例如：

```text
origin_um       = [0, 5, 10]
spacing_um      = 2
shape           = (51, 41)
fixed_y_um      = 5
```

表示网格从三维位置

$$
(0,5,10)
$$

开始，在 $X$ 和 $Z$ 方向每隔

$$
2\,\mu\mathrm m
$$

放置一个网格点，共有

$$
51\times41=2091
$$

个网格位置。

最终可以把整个函数理解为：

$$
\boxed{ \text{血管端点} \longrightarrow \text{血管范围} \longrightarrow \text{增加边缘空白} \longrightarrow \text{计算网格数量} \longrightarrow \text{生成二维计算区域} }
$$

它本身不计算血液怎样流动，只负责准备后续血管形状绘制和血流求解所需的“方格纸”。
