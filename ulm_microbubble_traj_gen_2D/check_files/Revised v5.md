> From: https://chatgpt.com/c/6a54af07-a098-83ee-a4ff-ccb4ef2b52ad

# chatgpt response

# 整个血管内的迁移率矩阵建模方案

建议将当前的粒子推进方式

$$
\mathbf x_i^{\,n+1}
=
\mathbf x_i^{\,n}
+
\mathbf u(\mathbf x_i^{\,n})\Delta t
$$

替换为一个**位置相关、尺寸相关的全血管迁移率模型**：

$$
\boxed{
\text{微泡速度}
=
\text{背景血流产生的速度}
+
\text{迁移率矩阵}\times\text{微泡碰撞力}
}
$$

这里不是为整个血管只设置一个固定矩阵，而是对每个微泡、在每个时间步，根据它所在的位置、半径和离壁距离，重新构造一个局部迁移率矩阵：

$$
\boxed{
\mathbf M_i
=
\mathbf M
\left(
\mathbf x_i,R_i,\mu,h_i,\mathbf n_i
\right)
}
$$

你的现有程序已经在二维连续速度场中逐位置插值流速，并保存连续壁面距离和法向，因此不需要重新建立血管几何和流场，只需替换微泡运动计算部分。

---

## 一、先固定模型的物理含义

虽然轨迹只在 $X-Z$ 平面中计算，但建议把微泡理解为：

> 一个三维刚性球形微泡，其球心运动被限制在二维 $X-Z$ 平面内。

因此，微泡具有三个运动自由度：

$$
\boxed{
\mathbf U_i
=
\begin{pmatrix}
V_{x,i}\\
V_{z,i}\\
\Omega_{y,i}
\end{pmatrix}
}
$$

其中：

- $V_x,V_z$ 是二维平移速度；
- $\Omega_y$ 是绕垂直于 $X-Z$ 平面的 $Y$ 轴旋转的角速度。

对应的广义力为：

$$
\boxed{
\mathbf G_i
=
\begin{pmatrix}
F_{x,i}\\
F_{z,i}\\
T_{y,i}
\end{pmatrix}
}
$$

这种处理与目标论文中“刚性球形微泡—低雷诺数流动—迁移率矩阵”的基本假设一致。论文采用迁移率矩阵直接由力和力矩求平移、旋转速度，而不是长期积分微泡加速度。

当前阶段明确忽略：

$$
\mathbf F^{\mathrm{buoyancy}}=0,
\qquad
\mathbf F^{\mathrm{acoustic}}=0,
\qquad
\mathbf F^{\mathrm{adhesion}}=0.
$$

只保留：

1. 背景血流；
2. 微泡之间的碰撞；
3. 近壁水动力修正；
4. 现有的血管壁不可穿透约束。

为了单独验证迁移率矩阵，建议暂时不要同时修改连续灌注、黏附或其他模块，以免无法判断轨迹变化来自哪一部分。

---

## 二、整个模型的核心方程

对于第 $i$ 个微泡，推荐采用：

$$
\boxed{
\mathbf U_i
=
\mathbf U_{H,i}
+
\mathbf M_i
\mathbf G_i^{\mathrm{coll}}
}
$$

其中：

- $\mathbf U_{H,i}$：没有微泡碰撞时，背景血流使微泡产生的平移和旋转速度；
- $\mathbf M_i$：该位置的迁移率矩阵；
- $\mathbf G_i^{\mathrm{coll}}$：其他微泡对它产生的碰撞力和力矩。

对于圆形微泡，如果碰撞力沿两个圆心的连线作用，碰撞不产生额外力矩，因此：

$$
\mathbf G_i^{\mathrm{coll}}
=
\begin{pmatrix}
F_{x,i}^{\mathrm{coll}}\\
F_{z,i}^{\mathrm{coll}}\\
0
\end{pmatrix}.
$$

需要注意，血流“推力”不应再被单独任意设置成一个常数。背景流已经由现有流场

$$
\mathbf u(\mathbf x)
$$

给出。迁移率矩阵主要负责：

- 在主体流区，把碰撞力转化为相对血流的附加速度；
- 在近壁区，同时修正血流驱动、碰撞响应以及平移—旋转耦合。

---

## 三、需要预先计算的场

在流场求解完成后，预先保存以下网格场。

### 1. 速度场

$$
\mathbf u(x,z)
=
\begin{pmatrix}
u_x(x,z)\\
u_z(x,z)
\end{pmatrix}.
$$

### 2. 速度梯度场

$$
\nabla\mathbf u
=
\begin{pmatrix}
\dfrac{\partial u_x}{\partial x}
&
\dfrac{\partial u_x}{\partial z}
\\[3mm]
\dfrac{\partial u_z}{\partial x}
&
\dfrac{\partial u_z}{\partial z}
\end{pmatrix}.
$$

它用于计算：

- 局部流动旋转；
- 近壁剪切率；
- 分叉处流动方向变化。

### 3. 连续壁面距离场

设微泡球心到最近固体壁面的距离为：

$$
d_w(\mathbf x_i).
$$

微泡表面到壁面的真实间隙为：

$$
\boxed{
h_i
=
d_w(\mathbf x_i)-R_i
}
$$

其中 $R_i$ 为微泡半径。

### 4. 壁面法向场

$$
\boxed{
\mathbf n_i
=
\frac{
\nabla d_w(\mathbf x_i)
}{
\left\|\nabla d_w(\mathbf x_i)\right\|
}
}
$$

这里 $\mathbf n_i$ 指向血管内部。

你的程序目前已经能在亚网格位置插值得到连续壁面距离和法向，所以这一部分可以直接复用。

---

## 四、建立局部切向—法向坐标

迁移率矩阵在靠近壁面时，应在局部壁面坐标中计算。

若壁面内法向为：

$$
\mathbf n_i
=
\begin{pmatrix}
n_{x,i}\\
n_{z,i}
\end{pmatrix},
$$

定义局部切向：

$$
\boxed{
\mathbf t_i
=
\begin{pmatrix}
-n_{z,i}\\
n_{x,i}
\end{pmatrix}
}
$$

这样有：

$$
\mathbf t_i\cdot\mathbf n_i=0.
$$

局部坐标到全局坐标的变换矩阵为：

$$
\boxed{
\mathbf Q_i
=
\begin{pmatrix}
t_{x,i}&n_{x,i}&0\\
t_{z,i}&n_{z,i}&0\\
0&0&1
\end{pmatrix}
}
$$

局部速度和全局速度的关系为：

$$
\mathbf U_i^{\mathrm{global}}
=
\mathbf Q_i
\mathbf U_i^{\mathrm{local}}.
$$

局部力为：

$$
\mathbf G_i^{\mathrm{local}}
=
\mathbf Q_i^{T}
\mathbf G_i^{\mathrm{global}}.
$$

局部迁移率矩阵转换到全局坐标后为：

$$
\boxed{
\mathbf M_i^{\mathrm{global}}
=
\mathbf Q_i
\mathbf M_i^{\mathrm{local}}
\mathbf Q_i^{T}
}
$$

---

## 五、主体流区的迁移率矩阵

远离壁面时，微泡近似处于无限流体中。

球形微泡的自由空间平移迁移率为：

$$
\mu_T^\infty
=
\frac{1}{6\pi\mu R_i}.
$$

旋转迁移率为：

$$
\mu_R^\infty
=
\frac{1}{8\pi\mu R_i^3}.
$$

因此自由空间迁移率矩阵为：

$$
\boxed{
\mathbf M_i^\infty
=
\begin{pmatrix}
\dfrac{1}{6\pi\mu R_i}&0&0\\[3mm]
0&\dfrac{1}{6\pi\mu R_i}&0\\[3mm]
0&0&\dfrac{1}{8\pi\mu R_i^3}
\end{pmatrix}
}
$$

为了避免矩阵中同时出现 $R^{-1}$、$R^{-2}$ 和 $R^{-3}$ 导致数值尺度差别过大，建议在程序内部使用缩放变量：

$$
\widetilde{\mathbf U}_i
=
\begin{pmatrix}
V_{t,i}\\
V_{n,i}\\
R_i\Omega_i
\end{pmatrix},
$$

$$
\widetilde{\mathbf G}_i
=
\begin{pmatrix}
F_{t,i}\\
F_{n,i}\\
T_i/R_i
\end{pmatrix}.
$$

此时：

$$
\boxed{
\widetilde{\mathbf U}_i
=
\frac{1}{6\pi\mu R_i}
\widehat{\mathbf M}_i
\widetilde{\mathbf G}_i
}
$$

自由空间的无量纲矩阵变为：

$$
\boxed{
\widehat{\mathbf M}^{\infty}
=
\begin{pmatrix}
1&0&0\\
0&1&0\\
0&0&\dfrac34
\end{pmatrix}
}
$$

这一形式更适合代码实现和矩阵条件数检查。

---

## 六、近壁迁移率矩阵

定义无量纲间隙：

$$
\boxed{
\xi_i=\frac{h_i}{R_i}.
}
$$

目标论文的近壁迁移率矩阵描述了：

- 切向力产生切向移动；
- 法向力产生法向移动；
- 切向力可以引起旋转；
- 旋转力矩也可以引起切向移动；
- 微泡越靠近壁面，法向迁移率越小。

论文中的完整三维迁移率矩阵及其近壁修正来自靠近平面壁面的球体 Stokes 水动力结果。

在你的二维模型中，只保留：

$$
V_t,\qquad V_n,\qquad \Omega_y.
$$

按照目标论文补充材料中的记号，定义四个近壁系数：

$$
C_{Ft}
=
\frac{8}{15}\ln\xi-0.9588,
$$

$$
C_{Tt}
=
\frac{1}{10}\ln\xi-0.1895,
$$

$$
C_{Fr}
=
\frac{2}{15}\ln\xi-0.2526,
$$

$$
C_{Tr}
=
\frac{2}{5}\ln\xi-0.3187.
$$

再定义：

$$
\boxed{
\Delta_h
=
C_{Tt}C_{Fr}
-
C_{Ft}C_{Tr}.
}
$$

法向壁面修正因子采用论文中的近似：

$$
\boxed{
\Lambda
\approx
1+\frac{1}{\xi}
=
1+\frac{R_i}{h_i}.
}
$$

这些系数和 $\Lambda$ 的出处是目标论文补充材料中的近壁水动力模型。

使用前述缩放变量后，二维近壁无量纲迁移率矩阵可写为：

$$
\boxed{
\widehat{\mathbf M}^{\,w}(\xi)
=
\begin{pmatrix}
\dfrac{C_{Tr}}{\Delta_h}
&
0
&
\dfrac{3C_{Fr}}{4\Delta_h}
\\[4mm]
0
&
\dfrac{1}{\Lambda}
&
0
\\[4mm]
\dfrac{C_{Tt}}{\Delta_h}
&
0
&
\dfrac{3C_{Ft}}{4\Delta_h}
\end{pmatrix}
}
$$

于是：

$$
\boxed{
\begin{pmatrix}
V_t\\
V_n\\
R\Omega
\end{pmatrix}
=
\frac{1}{6\pi\mu R}
\widehat{\mathbf M}^{\,w}
\begin{pmatrix}
F_t\\
F_n\\
T/R
\end{pmatrix}
}
$$

这里有一个很重要的数值检查：

$$
C_{Tt}
\approx
\frac34 C_{Fr}.
$$

因此矩阵的两个平移—旋转交叉项应当非常接近。如果程序中两者差异很大，通常意味着公式、符号或单位实现有错误。

---

## 七、不要把近壁渐近公式直接用到血管中央

上述对数形式主要适用于：

$$
\xi=\frac hR\ll1.
$$

因此不能在：

$$
h\gg R
$$

时仍直接使用这些公式。推荐构造一个从近壁矩阵平滑过渡到自由空间矩阵的全域模型。

给出三个初始数值参数：

$$
\xi_{\min}=10^{-3},
$$

$$
\xi_{\mathrm{near}}=0.1,
$$

$$
\xi_{\mathrm{far}}=1.0.
$$

这些是建议的第一版工程参数，不是普适物理常数，后续必须做敏感性分析。

首先截断极小间隙：

$$
\boxed{
\xi_{\mathrm{eff}}
=
\max(\xi,\xi_{\min}).
}
$$

近壁公式只在：

$$
\xi\leq \xi_{\mathrm{near}}
$$

时按实际 $\xi$ 计算。在过渡区中，可固定使用：

$$
\xi^\star
=
\min(\xi_{\mathrm{eff}},\xi_{\mathrm{near}})
$$

计算近壁矩阵，避免把小间隙渐近公式外推到 $\xi\sim1$。

定义：

$$
s
=
\operatorname{clip}
\left[
\frac{
\xi-\xi_{\mathrm{near}}
}{
\xi_{\mathrm{far}}-\xi_{\mathrm{near}}
},
0,1
\right].
$$

平滑近壁权重为：

$$
\boxed{
w(\xi)
=
1-3s^2+2s^3.
}
$$

于是：

- 当 $\xi\leq0.1$ 时，$w=1$；
- 当 $\xi\geq1$ 时，$w=0$；
- 中间区域平滑变化。

全血管有效无量纲迁移率矩阵为：

$$
\boxed{
\widehat{\mathbf M}^{\mathrm{eff}}
=
\left(1-w\right)
\widehat{\mathbf M}^{\infty}
+
w
\widehat{\mathbf M}^{\,w}
}
$$

这就是你需要在整个血管中使用的位置相关迁移率矩阵。

---

## 八、背景血流产生的基础运动

### 1. 主体流区

微泡当前位置的背景平移速度为：

$$
\mathbf V_i^{\mathrm{bulk}}
=
\mathbf u(\mathbf x_i).
$$

背景流引起的角速度为流体涡量的一半：

$$
\boxed{
\Omega_i^{\mathrm{bulk}}
=
\frac12
\left(
\frac{\partial u_x}{\partial z}
-
\frac{\partial u_z}{\partial x}
\right)_{\mathbf x_i}
}
$$

因此：

$$
\mathbf U_i^{\mathrm{bulk}}
=
\begin{pmatrix}
u_x(\mathbf x_i)\\
u_z(\mathbf x_i)\\
\Omega_i^{\mathrm{bulk}}
\end{pmatrix}.
$$

第一版建议先不加入 Faxén 的有限尺寸二阶修正，因为它需要计算速度拉普拉斯，而当前网格上的二阶导数可能放大噪声。先验证基础迁移率后再加入。

### 2. 近壁区

将局部血流剪切率定义为：

$$
\boxed{
S_i
=
\mathbf t_i^{T}
\left(\nabla\mathbf u_i\right)
\mathbf n_i.
}
$$

它表示沿壁方向的流速，在法向上变化得有多快。

球心到壁面的距离为：

$$
H_i=R_i+h_i.
$$

按照目标论文，近壁剪切流产生的切向等效力为：

$$
\boxed{
F_{t,i}^{\mathrm{shear}}
=
F_s
6\pi\mu R_iH_iS_i
}
$$

旋转力矩为：

$$
\boxed{
T_i^{\mathrm{shear}}
=
T_s
4\pi\mu R_i^3S_i
}
$$

其中目标论文采用：

$$
F_s=1.7005,
\qquad
T_s=0.9440.
$$

这里采用 Goldman、Cox 与 Brenner（1967）Part II 在球体接触平面壁面极限下的原始结果。下游补充材料中的 $F_s=1.0075$ 视为对 Goldman 原值 $1.7005$ 的转录错误；$T_s=0.9440$ 保持不变。

这些公式直接来自论文中剪切流推动微泡和使其旋转的近壁模型。

近壁广义剪切载荷为：

$$
\widetilde{\mathbf G}_i^{\mathrm{shear}}
=
\begin{pmatrix}
F_{t,i}^{\mathrm{shear}}\\
0\\
T_i^{\mathrm{shear}}/R_i
\end{pmatrix}.
$$

近壁基础速度为：

$$
\boxed{
\widetilde{\mathbf U}_i^{\mathrm{wall}}
=
\frac{1}{6\pi\mu R_i}
\widehat{\mathbf M}_i^{\,w}
\widetilde{\mathbf G}_i^{\mathrm{shear}}
}
$$

---

## 九、将主体流动和近壁流动统一起来

先将主体流动转换到局部坐标：

$$
u_{t,i}
=
\mathbf t_i\cdot\mathbf u_i,
$$

$$
u_{n,i}
=
\mathbf n_i\cdot\mathbf u_i.
$$

对应的缩放主体速度为：

$$
\widetilde{\mathbf U}_i^{\mathrm{bulk}}
=
\begin{pmatrix}
u_{t,i}\\
u_{n,i}\\
R_i\Omega_i^{\mathrm{bulk}}
\end{pmatrix}.
$$

全域背景水动力速度定义为：

$$
\boxed{
\widetilde{\mathbf U}_{H,i}
=
\left(1-w_i\right)
\widetilde{\mathbf U}_i^{\mathrm{bulk}}
+
w_i
\widetilde{\mathbf U}_i^{\mathrm{wall}}
}
$$

因此：

- 远离壁面时：

  $$
  \widetilde{\mathbf U}_{H,i}
  =
  \widetilde{\mathbf U}_i^{\mathrm{bulk}};
  $$
- 极靠近壁面时：

  $$
  \widetilde{\mathbf U}_{H,i}
  =
  \widetilde{\mathbf U}_i^{\mathrm{wall}}.
  $$

这一步使整个血管中不会出现“突然从流场速度切换到近壁速度”的跳变。

---

## 十、微泡碰撞力如何进入迁移率矩阵

对于微泡 $i$ 和 $j$，圆心距离为：

$$
r_{ij}
=
\left\|
\mathbf x_i-\mathbf x_j
\right\|.
$$

定义一个很小的碰撞作用层 $\ell_c$，则压缩量为：

$$
\boxed{
\delta_{ij}
=
\max
\left[
0,
R_i+R_j+\ell_c-r_{ij}
\right].
}
$$

方向为：

$$
\boxed{
\mathbf n_{ij}
=
\frac{
\mathbf x_i-\mathbf x_j
}{
r_{ij}
}.
}
$$

碰撞力使用软排斥形式：

$$
\boxed{
\mathbf F_{ij}^{\mathrm{coll}}
=
k_{ij}\delta_{ij}\mathbf n_{ij}.
}
$$

并满足：

$$
\mathbf F_{ji}^{\mathrm{coll}}
=
-\mathbf F_{ij}^{\mathrm{coll}}.
$$

第 $i$ 个微泡的总碰撞力为：

$$
\boxed{
\mathbf F_i^{\mathrm{coll}}
=
\sum_{j\neq i}
\mathbf F_{ij}^{\mathrm{coll}}.
}
$$

碰撞刚度不建议任意给一个极大数，而可以根据迁移率自动设置。

定义两个微泡沿碰撞方向的相对迁移率：

$$
\boxed{
\mu_{ij}^{\mathrm{rel}}
=
\mathbf n_{ij}^{T}
\left(
\mathbf M_{VF,i}^{\mathrm{global}}
+
\mathbf M_{VF,j}^{\mathrm{global}}
\right)
\mathbf n_{ij}.
}
$$

希望重叠在时间尺度 $\tau_c$ 内消除，则设置：

$$
\boxed{
k_{ij}
=
\frac{1}
{\mu_{ij}^{\mathrm{rel}}\tau_c}.
}
$$

第一版可以使用：

$$
\tau_c=5\Delta t
\quad\text{至}\quad
10\Delta t.
$$

这样，碰撞力会根据微泡大小、黏度和近壁阻力自动调整，而不是所有位置都使用相同刚度。

对于圆形微泡的中心碰撞：

$$
T_i^{\mathrm{coll}}=0.
$$

将总碰撞力转到局部坐标：

$$
F_{t,i}^{\mathrm{coll}}
=
\mathbf t_i\cdot
\mathbf F_i^{\mathrm{coll}},
$$

$$
F_{n,i}^{\mathrm{coll}}
=
\mathbf n_i\cdot
\mathbf F_i^{\mathrm{coll}}.
$$

于是：

$$
\widetilde{\mathbf G}_i^{\mathrm{coll}}
=
\begin{pmatrix}
F_{t,i}^{\mathrm{coll}}\\
F_{n,i}^{\mathrm{coll}}\\
0
\end{pmatrix}.
$$

碰撞造成的附加速度为：

$$
\boxed{
\Delta
\widetilde{\mathbf U}_i^{\mathrm{coll}}
=
\frac{1}{6\pi\mu R_i}
\widehat{\mathbf M}_i^{\mathrm{eff}}
\widetilde{\mathbf G}_i^{\mathrm{coll}}
}
$$

目标论文主要处理单颗微泡，因此没有完整考虑多微泡碰撞和多颗粒水动力耦合；你的碰撞模块属于在原论文基础上的扩展。

---

## 十一、最终统一运动方程

局部缩放速度为：

$$
\boxed{
\widetilde{\mathbf U}_i
=
\widetilde{\mathbf U}_{H,i}
+
\frac{1}{6\pi\mu R_i}
\widehat{\mathbf M}_i^{\mathrm{eff}}
\widetilde{\mathbf G}_i^{\mathrm{coll}}
}
$$

展开为：

$$
\widetilde{\mathbf U}_i
=
\begin{pmatrix}
V_{t,i}\\
V_{n,i}\\
R_i\Omega_i
\end{pmatrix}.
$$

因此：

$$
\Omega_i
=
\frac{
\widetilde U_{i,3}
}{
R_i
}.
$$

全局平移速度为：

$$
\boxed{
\mathbf V_i
=
V_{t,i}\mathbf t_i
+
V_{n,i}\mathbf n_i.
}
$$

最终位置更新：

$$
\boxed{
\mathbf x_i^{\,n+1}
=
\mathbf x_i^{\,n}
+
\mathbf V_i^{\,n}\Delta t.
}
$$

旋转角度更新：

$$
\boxed{
\theta_i^{\,n+1}
=
\theta_i^{\,n}
+
\Omega_i^{\,n}\Delta t.
}
$$

建议采用二阶 Heun 方法，而不是只使用显式 Euler：

$$
\mathbf x_i^\ast
=
\mathbf x_i^n
+
\mathbf V_i^n\Delta t,
$$

在预测位置重新计算：

$$
\mathbf V_i^\ast,
$$

最后：

$$
\boxed{
\mathbf x_i^{n+1}
=
\mathbf x_i^n
+
\frac{\Delta t}{2}
\left(
\mathbf V_i^n+\mathbf V_i^\ast
\right).
}
$$

现有的连续路径检查、二分定位碰壁点和切向滑移仍然保留，作为防止数值穿壁的最后保护层。

---

## 十二、多个微泡对应的“整个系统迁移率矩阵”

如果当前有 $N$ 个微泡，那么每个微泡都有自己的 $3\times3$ 矩阵：

$$
\mathbf M_1,
\mathbf M_2,
\ldots,
\mathbf M_N.
$$

第一版暂时不计算微泡之间的远程水动力扰动，因此整个系统矩阵采用分块对角形式：

$$
\boxed{
\mathbb M
=
\begin{pmatrix}
\mathbf M_1&0&\cdots&0\\
0&\mathbf M_2&\cdots&0\\
\vdots&\vdots&\ddots&\vdots\\
0&0&\cdots&\mathbf M_N
\end{pmatrix}.
}
$$

整个系统的运动写为：

$$
\boxed{
\mathbb U
=
\mathbb U_H
+
\mathbb M
\mathbb G^{\mathrm{coll}}.
}
$$

虽然矩阵是分块对角的，但微泡并不是完全独立，因为碰撞力满足：

$$
\mathbf F_{ij}^{\mathrm{coll}}
=
-\mathbf F_{ji}^{\mathrm{coll}},
$$

所以一个微泡的位置会通过碰撞力影响另一个微泡。

以后若考虑多微泡水动力耦合，则需要加入非零的非对角块：

$$
\mathbf M_{ij}\neq0.
$$

但这不是当前第一版应实现的内容。

---

## 十三、分叉处如何处理

不需要为分叉额外设计随机选择规则。

微泡始终按照：

$$
\frac{\mathrm d\mathbf x_i}{\mathrm dt}
=
\mathbf V_i
$$

连续推进。

在主体流区且没有碰撞时：

$$
\mathbf V_i
\approx
\mathbf u(\mathbf x_i),
$$

因此轨迹与现有流场产生的流线接近。

发生碰撞或近壁修正后：

$$
\mathbf V_i
\neq
\mathbf u(\mathbf x_i),
$$

微泡可以偏离原流线，但仍由连续运动方程自然进入某一支路。

因此：

$$
\boxed{
\text{分支选择是迁移率方程积分的结果，不是额外设定。}
}
$$

你现有模型本来就是根据连续速度场推进，而不是按照血管中心线或血管编号直接指定路径，这一逻辑应保持不变。

---

## 十四、程序中建议增加的函数

建议将实现拆分为以下函数。

```text
interpolate_flow(position)
    返回 u_x, u_z, grad_u, viscosity

interpolate_wall_geometry(position)
    返回 wall_distance, wall_normal

build_local_frame(wall_normal)
    返回 tangent, normal, Q

mobility_free(radius, viscosity)
    返回自由空间迁移率矩阵

mobility_wall(gap, radius, viscosity)
    返回目标论文近壁迁移率矩阵

mobility_effective(gap, radius, viscosity)
    计算平滑权重并返回全域有效矩阵

background_hydrodynamic_velocity(...)
    计算主体流速度、近壁剪切速度及二者过渡

build_neighbor_list(...)
    查找可能碰撞的微泡

compute_collision_forces(...)
    计算所有成对碰撞力

advance_bubbles(...)
    根据总速度进行 Heun 时间推进

check_continuous_path(...)
    执行现有的连续路径和防穿壁检查
```

迁移率矩阵只依赖：

$$
\xi=\frac hR
$$

的无量纲部分可以预先建立查找表。例如在：

$$
10^{-3}
\leq\xi
\leq0.1
$$

范围内使用对数间隔采样，保存：

$$
\widehat{\mathbf M}^{\,w}(\xi).
$$

运行时只需要插值，然后乘上：

$$
\frac{1}{6\pi\mu R}.
$$

这样可以显著减少每个时间步反复计算对数和矩阵系数的成本。

---

## 十五、推荐的实际实现顺序

### 第一步：只实现自由空间迁移率

暂时关闭近壁修正和碰撞：

$$
w=0,
\qquad
\mathbf F^{\mathrm{coll}}=0.
$$

此时必须恢复：

$$
\boxed{
\mathbf V_i
=
\mathbf u(\mathbf x_i)
}
$$

因此新程序得到的轨迹应与旧程序基本一致。

这一测试证明新的软件结构没有改变原有流场推进逻辑。

### 第二步：加入近壁矩阵，但仍关闭碰撞

启用：

$$
\widehat{\mathbf M}^{\mathrm{eff}}.
$$

检查：

- 血管中央轨迹与旧模型接近；
- 微泡靠近壁面时切向速度平滑变化；
- 法向迁移率随 $h/R$ 减小而下降；
- 不出现速度突然跳变。

### 第三步：加入单对微泡碰撞

只放两个微泡，使用均匀直通道，检查：

- 两个微泡不会持续重叠；
- 碰撞力大小相等、方向相反；
- 碰撞后仍随背景流向前运动；
- 没有不合理的惯性反弹。

### 第四步：加入多微泡和分叉

恢复完整血管树和多个微泡，检查：

- 分支选择仍由连续轨迹产生；
- 低浓度无碰撞结果接近原模型；
- 提高浓度后，碰撞允许部分微泡跨越原流线；
- 不出现由局部坐标切换引起的轨迹折线。

---

## 十六、必须通过的验证标准

### 1. 自由空间极限

当：

$$
\frac hR\ge1
$$

时，应有：

$$
\widehat{\mathbf M}^{\mathrm{eff}}
\approx
\widehat{\mathbf M}^{\infty}.
$$

### 2. 近壁法向迁移率

当：

$$
h\rightarrow0
$$

时：

$$
\Lambda
=
1+\frac Rh
\rightarrow\infty,
$$

因此：

$$
\boxed{
M_{nn}\rightarrow0.
}
$$

### 3. 矩阵正定性

对任意测试载荷 $\widetilde{\mathbf G}$，应满足：

$$
\boxed{
\widetilde{\mathbf G}^{T}
\widehat{\mathbf M}^{\mathrm{eff}}
\widetilde{\mathbf G}
\ge0.
}
$$

如果出现负值，说明矩阵符号、坐标方向或系数实现有误。

### 4. 平移—旋转互易关系

应近似满足：

$$
\boxed{
\widehat M_{13}
\approx
\widehat M_{31}.
}
$$

### 5. 无碰撞直通道

关闭碰撞后，血管中央应满足：

$$
\mathbf V_i
\approx
\mathbf u(\mathbf x_i).
$$

### 6. 对称分叉

在完全对称分叉中，大量低浓度微泡的左右分支比例应接近：

$$
P_L\approx P_R\approx0.5.
$$

### 7. 时间步收敛

分别使用：

$$
\Delta t,
\qquad
\frac{\Delta t}{2},
\qquad
\frac{\Delta t}{4}
$$

比较：

- 轨迹终点；
- 分支比例；
- 最大微泡重叠；
- 近壁停留时间。

结果应逐渐收敛。

---

## 十七、当前方案的适用边界

这一版是：

$$
\boxed{
\text{全血管背景流场}
+
\text{单壁近场修正}
+
\text{分块对角迁移率}
+
\text{短程微泡碰撞}
}
$$

它暂时没有包含：

- 两个壁面对同一微泡的精确联合水动力作用；
- 微泡之间的远程流体扰动；
- 红细胞作用；
- 微泡形变；
- 黏附力；
- 声学力。

目标论文的矩阵主要针对一颗微泡靠近平面壁面的情况，而将其扩展为“主体流区—过渡区—近壁区”的全血管位置相关矩阵，是本方案对原论文的工程扩展。目标论文也明确把弯曲、分叉和多微泡水动力相互作用列为原模型未充分包含的因素。

如果某个微泡同时非常接近两侧壁面，例如第二近壁间隙也满足：

$$
\frac{h_2}{R}<1,
$$

第一版应记录：

```text
two_wall_warning = true
```

但暂时仍使用最近壁面模型。等单壁全域模型验证完成后，再增加双壁阻力叠加，而不应在第一版中同时引入。

---

## 最终应写入程序的核心更新式

整个实现最终可以浓缩为：

$$
\boxed{
\widetilde{\mathbf U}_i
=
(1-w_i)
\widetilde{\mathbf U}_i^{\mathrm{bulk}}
+
w_i
\widetilde{\mathbf U}_i^{\mathrm{wall}}
+
\frac{1}{6\pi\mu R_i}
\widehat{\mathbf M}_i^{\mathrm{eff}}
\widetilde{\mathbf G}_i^{\mathrm{coll}}
}
$$

然后：

$$
\boxed{
\mathbf V_i
=
V_{t,i}\mathbf t_i
+
V_{n,i}\mathbf n_i
}
$$

$$
\boxed{
\mathbf x_i^{n+1}
=
\mathbf x_i^n
+
\mathbf V_i\Delta t.
}
$$

这个版本能够保持原有二维连续流场和自然分支选择，同时使微泡大小、离壁距离、平移—旋转耦合以及微泡碰撞真正进入运动方程。

---

# chatgpt response

可以。**上述“全血管位置相关迁移率矩阵 + 微泡碰撞 + 流场插值 + 壁面约束”的实现仍然非常适合基于 Taichi 计算，而且没有必要放弃你现有的 Taichi 轨迹框架。**

原因是，这套算法的主要工作都是“对每个微泡重复执行相同的小规模计算”：

$$
\text{插值流场}
\rightarrow
\text{计算壁面间隙}
\rightarrow
\text{计算迁移率系数}
\rightarrow
\text{搜索邻居}
\rightarrow
\text{计算碰撞力}
\rightarrow
\text{更新位置}.
$$

Taichi 的 kernel 会把最外层循环自动并行到多核 CPU 或 GPU，因此可以让一个线程负责一颗微泡。([docs.taichi-lang.org](https://docs.taichi-lang.org/docs/master/kernel_function)) 你的现有程序已经采用二维连续流场插值和逐时间步拉格朗日推进，这部分可以直接保留。

当前真正需要改变的不是计算平台，而是粒子推进 kernel 内部的计算内容。

## 一、不要真的建立一个巨大的全局矩阵

假设有 $N$ 颗微泡，数学上可以把整个系统写成：

$$
\mathbb M
=
\begin{pmatrix}
\mathbf M_1&0&\cdots&0\\
0&\mathbf M_2&\cdots&0\\
\vdots&\vdots&\ddots&\vdots\\
0&0&\cdots&\mathbf M_N
\end{pmatrix},
$$

其中每个：

$$
\mathbf M_i\in\mathbb R^{3\times3}
$$

是第 $i$ 颗微泡的位置相关迁移率矩阵。

但是，程序中**不要分配一个 $3N\times3N$ 的大矩阵**。因为当前模型没有考虑微泡之间的远程水动力耦合，所以全局矩阵是分块对角的。Taichi 只需要并行执行：

```text
第0个线程 → 计算第0颗微泡的局部矩阵
第1个线程 → 计算第1颗微泡的局部矩阵
第2个线程 → 计算第2颗微泡的局部矩阵
……
```

每颗微泡独立计算：

$$
\boxed{
\widetilde{\mathbf U}_i
=
\widetilde{\mathbf U}_{H,i}
+
\frac{1}{6\pi\mu_iR_i}
\widehat{\mathbf M}_i^{\mathrm{eff}}
\widetilde{\mathbf G}_i^{\mathrm{coll}}
}
$$

即可。

Taichi 对 $2\times2$、$3\times3$、$4\times4$ 这类小矩阵比较合适，矩阵运算会在编译时展开；官方文档也明确把 $3\times3$ 列为适合使用的小矩阵规模。([docs.taichi-lang.org](https://docs.taichi-lang.org/docs/field))

不过，为了进一步提高速度，甚至不必在每个时间步创建一个完整的局部 `ti.Matrix`。可以直接用五个非零系数进行矩阵—向量乘法：

$$
\widehat{\mathbf M}^{\mathrm{eff}}
=
\begin{pmatrix}
m_{tt}&0&m_{tT}\\
0&m_{nn}&0\\
m_{\omega t}&0&m_{\omega T}
\end{pmatrix}.
$$

于是：

$$
\Delta V_t
=
\frac{1}{6\pi\mu R}
\left(
m_{tt}F_t
+
m_{tT}\frac{T}{R}
\right),
$$

$$
\Delta V_n
=
\frac{1}{6\pi\mu R}
m_{nn}F_n,
$$

$$
R\Delta\Omega
=
\frac{1}{6\pi\mu R}
\left(
m_{\omega t}F_t
+
m_{\omega T}\frac{T}{R}
\right).
$$

当前微泡中心碰撞不产生碰撞力矩，因此通常：

$$
T^{\mathrm{coll}}=0.
$$

此时碰撞速度修正进一步简化为：

$$
\Delta V_t
=
\frac{m_{tt}F_t^{\mathrm{coll}}}
{6\pi\mu R},
$$

$$
\Delta V_n
=
\frac{m_{nn}F_n^{\mathrm{coll}}}
{6\pi\mu R},
$$

$$
\Delta\Omega
=
\frac{m_{\omega t}F_t^{\mathrm{coll}}}
{6\pi\mu R^2}.
$$

这样比保存每颗微泡的九个矩阵元素更节省显存和内存带宽。

## 二、推荐的 Taichi 数据结构

永久保存在设备上的网格场可以写成：

```python
flow_velocity = ti.Vector.field(2, dtype=ti.f32, shape=(NX, NZ))
flow_gradient = ti.Matrix.field(2, 2, dtype=ti.f32, shape=(NX, NZ))

wall_distance = ti.field(dtype=ti.f32, shape=(NX, NZ))
wall_normal = ti.Vector.field(2, dtype=ti.f32, shape=(NX, NZ))
viscosity = ti.field(dtype=ti.f32, shape=(NX, NZ))
```

微泡状态可以写成：

```python
position = ti.Vector.field(2, dtype=ti.f32, shape=MAX_PARTICLES)
velocity = ti.Vector.field(2, dtype=ti.f32, shape=MAX_PARTICLES)
collision_force = ti.Vector.field(2, dtype=ti.f32, shape=MAX_PARTICLES)

radius = ti.field(dtype=ti.f32, shape=MAX_PARTICLES)
omega = ti.field(dtype=ti.f32, shape=MAX_PARTICLES)
active = ti.field(dtype=ti.i32, shape=MAX_PARTICLES)
bubble_uid = ti.field(dtype=ti.i64, shape=MAX_PARTICLES)
```

Taichi field 本身就是可由 Python 端和 Taichi kernel 共同访问的全局数据容器，并支持标量、向量、矩阵和结构体元素。([docs.taichi-lang.org](https://docs.taichi-lang.org/docs/field))

对于你的计算，建议优先采用这种“分开的字段”：

- 位置放在一个向量场；
- 速度放在一个向量场；
- 半径放在一个标量场；
- 有效状态放在一个整数场。

不建议把每颗微泡的全部信息都打包成一个很大的结构，因为某些 kernel 只需要位置和半径，而不需要读取其他状态。分开存储可以减少不必要的数据读取。

当前程序已经建立连续壁面距离场，并能在亚网格位置插值得到壁面间隙和壁面法向，这两项都可以原样保留。

## 三、迁移率矩阵计算本身不会成为主要性能瓶颈

对每颗微泡，迁移率计算通常只包括：

1. 双线性插值；
2. 计算
   $$
   h=d_w-R;
   $$
3. 计算
   $$
   \xi=\frac{h}{R};
   $$
4. 计算平滑权重 $w(\xi)$；
5. 读取或计算五个迁移率系数；
6. 完成一次很小的矩阵—向量乘法。

这些计算量很小。真正可能变慢的是：

$$
\boxed{\text{微泡—微泡邻居搜索和碰撞计算}}
$$

因为若直接让每颗微泡与所有其他微泡比较，则复杂度为：

$$
O(N^2).
$$

例如有 $10^5$ 颗微泡时，直接两两比较显然不可行。

因此，Taichi 高效实现的核心不是迁移率矩阵，而是建立二维空间分箱。

## 四、使用均匀网格建立微泡邻居表

设碰撞最大作用距离为：

$$
r_{\mathrm{cut}}
=
2R_{\max}+\ell_c,
$$

其中：

- $R_{\max}$ 是最大微泡半径；
- $\ell_c$ 是碰撞作用层厚度。

将碰撞网格单元尺寸设为：

$$
\boxed{
L_{\mathrm{cell}}
\geq
r_{\mathrm{cut}}.
}
$$

这样，一颗微泡只需要检查自己所在网格及周围八个网格：

$$
3\times3=9
$$

个单元，不需要检查整个血管中的所有微泡。

Taichi 字段可定义为：

```python
cell_count = ti.field(
    dtype=ti.i32,
    shape=(NCX, NCZ)
)

cell_particles = ti.field(
    dtype=ti.i32,
    shape=(NCX, NCZ, MAX_PER_CELL)
)

cell_overflow = ti.field(
    dtype=ti.i32,
    shape=()
)
```

每个时间步先执行：

```text
清空 cell_count
→ 把每颗微泡放入对应网格
→ 搜索周围九个网格
→ 计算碰撞力
```

对应的网格编号为：

$$
c_x
=
\left\lfloor
\frac{x-x_{\min}}{L_{\mathrm{cell}}}
\right\rfloor,
$$

$$
c_z
=
\left\lfloor
\frac{z-z_{\min}}{L_{\mathrm{cell}}}
\right\rfloor.
$$

把微泡加入网格时，可通过原子计数获得插入位置：

```python
slot = ti.atomic_add(cell_count[cx, cz], 1)
if slot < MAX_PER_CELL:
    cell_particles[cx, cz, slot] = i
else:
    ti.atomic_add(cell_overflow[None], 1)
```

必须保留 `cell_overflow`。如果它不为零，说明：

$$
\text{MAX\_PER\_CELL}
$$

设置过小，已有微泡没有进入邻居表，碰撞结果会漏算。

## 五、碰撞力建议采用“每颗微泡自己收集”的模式

有两种并行方式。

第一种方法是每对微泡只计算一次，然后同时写入：

$$
\mathbf F_i
\quad\text{和}\quad
\mathbf F_j.
$$

这需要大量原子加法。

第二种方法是每颗微泡单独遍历自己的邻居，只写自己的碰撞力：

```python
for i in range(MAX_PARTICLES):
    if active[i] == 1:
        force_i = ti.Vector([0.0, 0.0])

        for 周围九个网格:
            for 该网格中的微泡 j:
                if j != i:
                    force_i += pair_force(i, j)

        collision_force[i] = force_i
```

这样每一对微泡会被计算两次：

$$
i\rightarrow j,
\qquad
j\rightarrow i,
$$

但不同线程不会同时写同一个 `collision_force[i]`，从而减少原子竞争。

对 GPU 而言，**多算一次简单的两粒子距离，通常比大量争抢同一个原子写入位置更容易优化**。不过最终仍应通过 profiler 比较两种实现，不能仅凭经验决定。

## 六、迁移率系数最好使用查找表

近壁迁移率中包含：

$$
\ln\left(\frac hR\right)
$$

以及平滑过渡计算。如果每颗微泡、每个时间步都重新调用多个对数函数，仍可能产生一定开销。

可以预先建立一张一维查找表：

$$
\xi_k
=
10^{
\log_{10}\xi_{\min}
+
k\Delta
},
$$

并保存：

$$
m_{tt}(\xi_k),
\quad
m_{nn}(\xi_k),
\quad
m_{tT}(\xi_k),
\quad
m_{\omega t}(\xi_k),
\quad
m_{\omega T}(\xi_k).
$$

例如：

```python
mobility_lut = ti.Vector.field(
    5,
    dtype=ti.f32,
    shape=N_LUT
)
```

运行时根据：

$$
\xi=\frac hR
$$

找到相邻两个表项，再做线性插值：

$$
m(\xi)
\approx
(1-\alpha)m_k+\alpha m_{k+1}.
$$

这样，Taichi kernel 中不再需要反复计算多个对数表达式。

第一版也可以直接计算公式，因为矩阵计算本身通常不是主要瓶颈。建议顺序是：

1. 先直接计算公式，确保物理结果正确；
2. 用 profiler 测量；
3. 若迁移率 kernel 占时明显，再改为查找表。

## 七、推荐的 Taichi kernel 流程

一个时间步建议由四到五个 kernel 组成。

### 1. 清空碰撞网格和力

```python
@ti.kernel
def clear_step_data():
    for I in ti.grouped(cell_count):
        cell_count[I] = 0

    for i in range(MAX_PARTICLES):
        collision_force[i] = ti.Vector([0.0, 0.0])

    cell_overflow[None] = 0
```

### 2. 建立微泡网格表

```python
@ti.kernel
def build_cell_list():
    for i in range(MAX_PARTICLES):
        if active[i] == 1:
            cx, cz = position_to_cell(position[i])
            slot = ti.atomic_add(cell_count[cx, cz], 1)

            if slot < MAX_PER_CELL:
                cell_particles[cx, cz, slot] = i
            else:
                ti.atomic_add(cell_overflow[None], 1)
```

### 3. 计算碰撞力

```python
@ti.kernel
def compute_collision_forces():
    for i in range(MAX_PARTICLES):
        if active[i] == 1:
            f = ti.Vector([0.0, 0.0])
            cx, cz = position_to_cell(position[i])

            for ox, oz in ti.static(ti.ndrange((-1, 2), (-1, 2))):
                nx = cx + ox
                nz = cz + oz

                if valid_cell(nx, nz):
                    n_in_cell = ti.min(
                        cell_count[nx, nz],
                        MAX_PER_CELL
                    )

                    for k in range(n_in_cell):
                        j = cell_particles[nx, nz, k]

                        if j != i and active[j] == 1:
                            f += pair_collision_force(i, j)

            collision_force[i] = f
```

### 4. 计算迁移率并生成候选位置

```python
@ti.kernel
def predict_positions(dt: ti.f32):
    for i in range(MAX_PARTICLES):
        if active[i] == 1:
            u, grad_u = interpolate_flow(position[i])
            dw, normal = interpolate_wall(position[i])

            gap = dw - radius[i]

            vel, angular_vel = compute_mobility_velocity(
                position[i],
                radius[i],
                u,
                grad_u,
                gap,
                normal,
                collision_force[i]
            )

            velocity[i] = vel
            omega[i] = angular_vel
            candidate_position[i] = position[i] + dt * vel
```

### 5. 壁面检查并提交位置

```python
@ti.kernel
def enforce_wall_and_commit():
    for i in range(MAX_PARTICLES):
        if active[i] == 1:
            corrected = continuous_wall_check(
                position[i],
                candidate_position[i],
                radius[i]
            )
            position[i] = corrected
```

Taichi kernel 不能在另一个 kernel 内直接调用，因此这些 kernel 应由 Python 端按照固定顺序调度；已经编译的 kernel 会被缓存，后续调用不需要每次重新编译。([docs.taichi-lang.org](https://docs.taichi-lang.org/docs/master/kernel_function))

你的现有连续路径检查、碰壁点二分查找和切向滑移都可以继续使用，不需要因为引入迁移率矩阵而删除。

## 八、先继续用 Euler，再决定是否使用 Heun

你的现有轨迹更新是显式时间推进。 第一版迁移率实现建议继续使用：

$$
\mathbf x_i^{n+1}
=
\mathbf x_i^n
+
\mathbf V_i^n\Delta t.
$$

因为它只需要：

- 建立一次邻居表；
- 计算一次碰撞力；
- 计算一次迁移率速度。

如果使用严格的 Heun 方法：

$$
\mathbf x_i^\ast
=
\mathbf x_i^n+\mathbf V_i^n\Delta t,
$$

还要在预测位置上：

1. 重新构建粒子网格；
2. 重新计算碰撞力；
3. 重新插值流场；
4. 重新计算迁移率矩阵。

计算量接近翻倍。

因此推荐顺序为：

1. 使用 Euler 完成迁移率矩阵正确性验证；
2. 分别用

   $$
   \Delta t,\quad \frac{\Delta t}{2},\quad \frac{\Delta t}{4}
   $$

   做时间步收敛；
3. 如果 Euler 误差仍明显，再改为 Heun。

碰撞作用会产生较强的短时间尺度，因此应限制：

$$
\boxed{
\max_i
\left(
\|\mathbf V_i\|\Delta t
\right)
<
\eta
\min
\left(
\Delta x,
R_{\min},
\ell_c
\right)
}
$$

其中 $\eta$ 可以先取小于 $0.2$ 的安全系数，再通过收敛实验调整。

## 九、连续灌注也可以继续用 Taichi

以后加入连续灌注时，不建议不断扩大 Taichi 数组，而是预先设置：

$$
\text{MAX\_PARTICLES}
$$

个可复用的计算槽位。

需要区分：

$$
\boxed{
\text{计算槽位编号}
\neq
\text{永久微泡身份编号}
}
$$

例如：

- `slot = 37` 表示 GPU 上第 37 个存储位置；
- `bubble_uid = 105928` 表示这颗微泡的永久身份。

微泡离开以后，可以复用 `slot=37`，但新进入微泡获得新的 `bubble_uid`。这样不会因为连续灌注而无限增加设备数组。

你当前程序已经为微泡保存永久身份、出生帧和终止原因。 因此，只需要把“永久身份不复用”和“GPU 存储槽位不复用”分离开：

$$
\boxed{
\text{身份不可复用，存储槽位可以复用。}
}
$$

这会显著降低连续灌注时的设备内存压力。

## 十、不要每个时间步都把结果传回 NumPy

Taichi 场应一直保留在 GPU 或 CPU 设备上。不要在每一步执行：

```python
position_numpy = position.to_numpy()
```

因为这会产生频繁的设备—主机数据传输。

建议在设备上建立轨迹缓冲区：

```python
trajectory_buffer = ti.Vector.field(
    2,
    dtype=ti.f32,
    shape=(BUFFER_FRAMES, MAX_PARTICLES)
)
```

每隔若干模拟步保存一次：

$$
n\bmod s_{\mathrm{output}}=0.
$$

当缓冲区装满后，再一次性复制到 CPU 并写入文件。

Taichi 的 ndarray 可以与 NumPy 和 PyTorch 数据交换；如果外部数组与 Taichi 使用同一设备，可以避免额外的数据传输，而当 NumPy 在 CPU、Taichi 在 CUDA 时，则会发生设备间复制。官方文档也建议桌面 GPU 计算明确使用 32 位数据。([docs.taichi-lang.org](https://docs.taichi-lang.org/docs/master/ndarray))

## 十一、数值精度建议

生产计算可以优先采用：

```python
ti.init(
    arch=ti.cuda,
    default_fp=ti.f32,
    kernel_profiler=True
)
```

但不建议直接在程序中同时使用：

- 米；
- 牛顿；
- 帕秒；
- 微米级间隙。

因为在 `f32` 下，过大和过小数值混合会使矩阵计算条件较差。

建议统一使用一种缩放单位，例如：

$$
\text{长度：}\mu\mathrm m,
$$

$$
\text{时间：}\mathrm s,
$$

$$
\text{速度：}\mu\mathrm m/\mathrm s,
$$

$$
\text{力：}\mathrm{pN}.
$$

或者对整个系统进行无量纲化。

同时保留一个小规模 `f64` CPU 或 CUDA 验证版本。比较：

$$
\mathbf x^{f32}(T)
$$

和：

$$
\mathbf x^{f64}(T).
$$

若轨迹、分支比例和最小间隙差异可接受，再使用 `f32` 进行大规模计算。

## 十二、Taichi 中最值得优化的地方

迁移率 kernel 本身很轻。建议按以下优先级优化：

1. **邻居搜索是否仍有 $O(N^2)$ 循环；**
2. **粒子网格是否溢出；**
3. **轨迹数据是否频繁传回 CPU；**
4. **是否保存了不必要的每颗微泡 $3\times3$ 矩阵；**
5. **碰撞 kernel 是否存在大量原子写竞争；**
6. **壁面检查是否使用过多不定长循环；**
7. **迁移率对数公式是否值得改成查找表。**

Taichi 的外层循环会自动并行，必要时可以通过：

```python
ti.loop_config(block_dim=128)
```

或：

```python
ti.loop_config(block_dim=256)
```

调整 GPU block 大小，但应先使用默认设置，再根据实测调优。官方文档说明 `ti.loop_config` 可调整 CPU 并行线程数和 GPU block 大小。([docs.taichi-lang.org](https://docs.taichi-lang.org/docs/master/performance))

## 十三、必须使用 profiler，而不是猜测

可以这样初始化：

```python
ti.init(
    arch=ti.cuda,
    default_fp=ti.f32,
    kernel_profiler=True
)
```

运行若干预热步后：

```python
ti.sync()
ti.profiler.clear_kernel_profiler_info()

for _ in range(N_BENCHMARK_STEPS):
    simulation_step()

ti.sync()
ti.profiler.print_kernel_profiler_info("count")
```

Taichi 的 KernelProfiler 可以分别统计各个 kernel 的总时间、调用次数和平均时间；在 CUDA 后端还提供更深入的硬件性能分析。([docs.taichi-lang.org](https://docs.taichi-lang.org/docs/profiler))

重点观察：

```text
build_cell_list
compute_collision_forces
predict_positions
enforce_wall_and_commit
save_trajectory
```

各自占用多少时间。

通常预期是：

$$
\boxed{
\text{碰撞邻居搜索}
>
\text{壁面路径检查}
>
\text{迁移率矩阵计算}
}
$$

但这只是合理预期，最终必须以 profiler 结果为准。

## 十四、推荐采用的实现顺序

第一阶段只把原来的：

$$
\mathbf V_i=\mathbf u(\mathbf x_i)
$$

替换成新的 Taichi 函数结构，但暂时令：

$$
\mathbf F_i^{\mathrm{coll}}=0,
\qquad
\widehat{\mathbf M}^{\mathrm{eff}}
=
\widehat{\mathbf M}^{\infty}.
$$

此时新旧轨迹应基本一致。

第二阶段加入位置相关迁移率，但仍关闭碰撞，检查：

$$
h/R\gg1
$$

时能否恢复主体流场轨迹，以及：

$$
h/R\rightarrow0
$$

时法向迁移率是否下降。

第三阶段加入两颗微泡碰撞，验证：

$$
\mathbf F_{ij}
=
-\mathbf F_{ji},
$$

且不会持续重叠。

第四阶段加入均匀网格邻居搜索和大量微泡。

第五阶段再加入连续灌注和设备端轨迹缓冲。

这种顺序能保证每一步都可以与前一步比较，避免一次性修改太多模块。

## 最终结论

$$
\boxed{
\text{完全可以继续使用 Taichi，而且 Taichi 仍然是这套模型的合适实现工具。}
}
$$

具体来说：

- 流场、速度梯度、壁面距离和法向继续使用二维 Taichi field；
- 每颗微泡由一个并行线程计算局部迁移率；
- 不构造 $3N\times3N$ 全局矩阵；
- 不保存每颗微泡的完整 $3\times3$ 矩阵，而是直接计算五个非零系数；
- 微泡碰撞通过二维 cell list 搜索邻居；
- 轨迹推进、壁面检查和分支选择继续使用现有 Taichi kernel；
- 轨迹输出采用设备端缓冲，不在每一步传回 NumPy。

当前官方发布页将 Taichi v1.7.4 标记为最新版本，因此正式开发时建议固定具体版本并记录在环境文件中，避免不同版本导致 kernel 行为或性能变化。([GitHub](https://github.com/taichi-dev/taichi/releases))

这套实现的计算瓶颈不会是迁移率矩阵，而会是**碰撞邻居搜索、壁面路径检查和轨迹数据输出**。只要这三个部分采用并行网格搜索、固定容量设备字段和批量输出，迁移率模型可以继续在 Taichi 上进行高效的大规模计算。

---

Powered by [AI Exporter](https://saveai.net)
