
下面是一版可直接拆成 5 页 PPT 的“原理文字 + LaTeX 公式”。

## 第 1 页：瞬时运动趋势的核心思想

### PPT 正文

微泡不是简单地以局部血流速度随流运动。系统在每个时刻综合三类因素：

1. 血流和近壁剪切产生的背景平移与转动；
2. 微泡碰撞、分子黏附等载荷通过迁移率矩阵产生的附加运动；
3. 红细胞分布产生的确定性横向迁移。

因此，微泡的自由瞬时运动趋势可以概括为：

$$
\boxed{\mathbf U_{\mathrm{free}}=\mathbf U_{\mathrm{bg}}+\mathbf M\mathbf L_{\mathrm{ext}}+\mathbf U_{\mathrm{RBC}}}
$$

其中：

$$
\mathbf U=\begin{bmatrix}V_x\\V_z\\\Omega_y\end{bmatrix},\qquad \mathbf L_{\mathrm{ext}}=\begin{bmatrix}F_x\\F_z\\T_y\end{bmatrix}
$$

- $V_x,V_z$：微泡在二维血管平面内的平移速度；
- $\Omega_y$：微泡绕垂直于模拟平面的轴旋转的角速度；
- $\mathbf U_{\mathrm{bg}}$：背景血流和近壁剪切产生的运动；
- $\mathbf M$：受微泡尺寸、黏度和壁面距离控制的迁移率矩阵；
- $\mathbf L_{\mathrm{ext}}$：碰撞力、分子黏附力和黏附转矩；
- $\mathbf U_{\mathrm{RBC}}$：红细胞诱导的确定性漂移修正。

### 一句话讲解

> 系统先判断血流本身会怎样带动微泡，再用迁移率矩阵把各种力和转矩转换成附加速度，最终形成微泡当前时刻“想要怎样运动”的自由趋势。

---

## 第 2 页：为什么需要迁移率矩阵

### PPT 正文

在低雷诺数微尺度流动中，微泡的惯性可以忽略。微泡受到载荷后，不经历明显的“持续加速”过程，而是几乎立即形成与载荷对应的平移和转动速度。

因此，系统使用线性迁移率关系：

$$
\boxed{\text{速度响应}=\text{迁移率}\times\text{外部载荷}}
$$

以局部壁面坐标表示：

- $\mathbf n$：从血管壁指向管腔内部的单位法向；
- $\mathbf t=(-n_z,n_x)$：沿血管壁的单位切向；
- $R$：微泡半径；
- $g$：微泡表面到血管壁的间隙；
- $\xi=g/R$：无量纲壁面间隙。

迁移率的基础尺度为 Stokes 迁移率：

$$
\boxed{\mu_0=\frac{1}{6\pi\eta R}}
$$

其中，$\eta$ 为局部血液动力黏度。

### 局部迁移率矩阵

$$
\begin{bmatrix}\Delta V_t\\\Delta V_n\\R\Delta\Omega_y\end{bmatrix}=\mu_0\underbrace{\begin{bmatrix}m_{tt}&0&m_{tR}\\0&m_{nn}&0\\m_{Rt}&0&m_{RR}\end{bmatrix}}_{\widetilde{\mathbf M}(\xi)}\begin{bmatrix}F_t\\F_n\\T_y/R\end{bmatrix}
$$

### 各矩阵项的物理意义

| 矩阵项     | 物理含义                       |
| ---------- | ------------------------------ |
| $m_{tt}$ | 切向力引起沿壁平移             |
| $m_{nn}$ | 法向力引起远离或靠近壁面的运动 |
| $m_{tR}$ | 转矩同时引起沿壁平移           |
| $m_{Rt}$ | 切向力同时引起微泡转动         |
| $m_{RR}$ | 转矩引起微泡自转               |

### 一句话讲解

> 迁移率矩阵不仅描述“力产生速度”，还描述平移与转动之间的耦合：例如，沿壁拉动微泡时，微泡可能同时发生滚动。

---

## 第 3 页：壁面距离如何改变迁移率

### PPT 正文

微泡远离血管壁时，可近似采用自由空间迁移率：

$$
\widetilde{\mathbf M}_{\infty}=\begin{bmatrix}1&0&0\\0&1&0\\0&0&3/4\end{bmatrix}
$$

此时：

- 切向和法向平移能力基本相同；
- 平移和转动不发生明显耦合；
- 微泡主要跟随局部血流运动。

当微泡接近壁面时，微泡与壁面之间的液膜会产生额外阻力：

- 法向运动受到明显抑制；
- 沿壁滑动阻力增大；
- 切向平移与转动开始相互耦合；
- 同样大小的力不再产生同样大小的速度。

系统通过平滑权重在自由空间模型和近壁模型之间过渡：

$$
\boxed{\widetilde{\mathbf M}(\xi)=\left[1-w(\xi)\right]\widetilde{\mathbf M}_{\infty}+w(\xi)\widetilde{\mathbf M}_{\mathrm{wall}}(\hat{\xi})}
$$

其中：

$$
\hat{\xi}=\min\left[\max(\xi,\xi_{\min}),\xi_{\mathrm{near}}\right]
$$

$$
w(\xi)=\begin{cases}1,&\xi\leq\xi_{\mathrm{near}}\\1-3s^2+2s^3,&\xi_{\mathrm{near}}<\xi<\xi_{\mathrm{far}}\\0,&\xi\geq\xi_{\mathrm{far}}\end{cases}
$$

$$
s=\frac{\xi-\xi_{\mathrm{near}}}{\xi_{\mathrm{far}}-\xi_{\mathrm{near}}}
$$

### 产品层面的解释

- 距离血管壁较远：主要采用自由随流模型；
- 进入近壁影响区：逐渐增强壁面阻力和滚动耦合；
- 非常接近血管壁：采用近壁迁移率；
- 极小间隙会被正则化，避免迁移率公式出现数值奇异。

### 一句话讲解

> 迁移率不是固定常数，而是一套随壁面距离变化的“运动响应能力”：越接近壁面，微泡越难自由运动。

---

## 第 4 页：背景流、碰撞和分子黏附如何合成

### 4.1 背景流产生基础运动

远离血管壁时：

$$
\mathbf V_{\mathrm{bulk}}=\mathbf u(\mathbf x)
$$

$$
\Omega_{\mathrm{bulk}}=\frac12\left(\frac{\partial u_x}{\partial z}-\frac{\partial u_z}{\partial x}\right)
$$

其中，$\mathbf u(\mathbf x)$ 是微泡中心位置处的局部流体速度。

近壁时，系统计算局部剪切率：

$$
\dot{\gamma}=\mathbf t^{\mathsf T}\left(\nabla\mathbf u\right)\mathbf n
$$

剪切流在微泡上产生等效切向载荷和转矩：

$$
F_t^{\mathrm{shear}}=C_F\,6\pi\eta R\left[R+\max(g,0)\right]\dot{\gamma}
$$

$$
T_y^{\mathrm{shear}}=C_T\,4\pi\eta R^3\dot{\gamma}
$$

这些载荷通过近壁迁移率转换成沿壁滑动和旋转。背景运动最终写为：

$$
\mathbf U_{\mathrm{bg}}=\left[1-w(\xi)\right]\mathbf U_{\mathrm{bulk}}+w(\xi)\mathbf U_{\mathrm{wall}}
$$

### 4.2 外部载荷合成

系统将碰撞和分子黏附产生的载荷相加：

$$
\mathbf F_{\mathrm{ext}}=\mathbf F_{\mathrm{collision}}+\mathbf F_{\mathrm{bond}}
$$

$$
T_y^{\mathrm{ext}}=T_{\mathrm{bond}}
$$

其中：

- 碰撞力：在微泡过于接近时产生排斥，防止相互穿透；
- 黏附切向力：阻碍微泡沿壁滑动；
- 黏附法向力：改变微泡靠近或远离壁面的趋势；
- 黏附转矩：改变微泡滚动和旋转状态。

附加运动为：

$$
\boxed{\Delta\mathbf U_{\mathrm{ext}}=\mathbf M_{xz}\begin{bmatrix}F_x^{\mathrm{ext}}\\F_z^{\mathrm{ext}}\\T_y^{\mathrm{ext}}\end{bmatrix}}
$$

### 4.3 全局迁移率矩阵

局部迁移率会被旋转回全局 $X-Z$ 坐标：

$$
\boxed{\mathbf M_{xz}=\mu_0\begin{bmatrix}m_{tt}\mathbf t\mathbf t^{\mathsf T}+m_{nn}\mathbf n\mathbf n^{\mathsf T}&\dfrac{m_{tR}}{R}\mathbf t\\[8pt]\dfrac{m_{Rt}}{R}\mathbf t^{\mathsf T}&\dfrac{m_{RR}}{R^2}\end{bmatrix}}
$$

它建立了以下映射：

$$
\begin{bmatrix}F_x\\F_z\\T_y\end{bmatrix}\longrightarrow\begin{bmatrix}\Delta V_x\\\Delta V_z\\\Delta\Omega_y\end{bmatrix}
$$

### 一句话讲解

> 背景血流直接给出基础运动；碰撞和黏附先表现为力与转矩，再由同一个迁移率矩阵统一转换成附加平移和转动。

---

## 第 5 页：得到自由瞬时趋势，并进行壁面约束

### PPT 正文

综合所有确定性因素后：

$$
\boxed{\mathbf U_{\mathrm{free}}=\mathbf U_{\mathrm{bg}}+\mathbf M_{xz}\begin{bmatrix}\mathbf F_{\mathrm{collision}}+\mathbf F_{\mathrm{bond}}\\T_{\mathrm{bond}}\end{bmatrix}+\begin{bmatrix}\mathbf V_{\mathrm{RBC}}\\0\end{bmatrix}}
$$

这里的 $\mathbf U_{\mathrm{free}}$ 表示微泡在不考虑刚性壁面穿透限制时的瞬时运动趋势。

红细胞修正包含：

$$
\mathbf V_{\mathrm{RBC}}=\mathbf V_{\mathrm{drift}}+\mathbf V_{\mathrm{Fick}}
$$

其中：

- $\mathbf V_{\mathrm{drift}}$：红细胞诱导的确定性横向迁移；
- $\mathbf V_{\mathrm{Fick}}$：空间扩散系数变化产生的修正速度。

红细胞随机扩散不是迁移率矩阵中的力项，而是在确定性推进完成后，作为独立随机位移施加，并再次接受壁面、分叉和出口检查。

### 壁面约束

如果自由趋势会使微泡穿透血管壁，系统增加一个沿壁面法向的接触反力：

$$
\mathbf a=\begin{bmatrix}n_x\\n_z\\0\end{bmatrix}
$$

$$
\boxed{\mathbf U_{\mathrm{contact}}=\mathbf U_{\mathrm{free}}+\mathbf M_{xz}\mathbf a\,\lambda}
$$

接触反力大小为：

$$
\boxed{\lambda=\max\left[0,\,-\frac{g/\Delta t+\mathbf a^{\mathsf T}\mathbf U_{\mathrm{free}}}{\mathbf a^{\mathsf T}\mathbf M_{xz}\mathbf a}\right]}
$$

其约束关系可以表示为：

$$
0\leq\lambda\quad\perp\quad g+\Delta t\,\mathbf a^{\mathsf T}\mathbf U_{\mathrm{contact}}\geq0
$$

含义是：

- 如果微泡不会穿壁，则 $\lambda=0$，保留自由运动；
- 如果微泡存在穿壁趋势，则产生最小必要接触反力；
- 接触反力通过同一迁移率矩阵改变平移和转动；
- 修正后的微泡必须保持在有限尺寸可达的管腔内。

### 一句话讲解

> 瞬时运动趋势先由血流、迁移率载荷和红细胞漂移共同形成；如果该趋势导致穿壁，系统再施加最小必要接触反力，将自由趋势修正为物理可接受的运动。

---

## 总结页可用文字

### 标题

迁移率矩阵是“载荷”与“微泡运动”之间的核心转换器

### 正文

微泡的瞬时运动由“背景流动”和“载荷响应”共同决定。系统首先根据局部血流、速度梯度和壁面距离计算背景平移与转动；随后将碰撞力、分子黏附力和黏附转矩输入迁移率矩阵，获得相应的附加速度与角速度；再叠加红细胞诱导的确定性漂移，形成自由瞬时运动趋势。最后，系统通过不可穿透壁面约束对该趋势进行审核和修正，得到能够进入时间积分的实际运动状态。

### 最简公式链

$$
\boxed{\text{局部流场}\rightarrow\mathbf U_{\mathrm{bg}}}
$$

$$
\boxed{\text{碰撞与黏附载荷}\xrightarrow{\mathbf M}\Delta\mathbf U_{\mathrm{ext}}}
$$

$$
\boxed{\mathbf U_{\mathrm{free}}=\mathbf U_{\mathrm{bg}}+\Delta\mathbf U_{\mathrm{ext}}+\mathbf U_{\mathrm{RBC}}}
$$

$$
\boxed{\mathbf U_{\mathrm{free}}\xrightarrow{\text{壁面约束}}\mathbf U_{\mathrm{accepted}}}
$$

对应实现主要位于 [particle_mobility.py](E:/ULM/hatimb-particle_flow_simulator/ulm_microbubble_traj_gen/utils/particles/particle_mobility.py:151) 和 [particle_mobility_transport.py](E:/ULM/hatimb-particle_flow_simulator/ulm_microbubble_traj_gen/utils/particles/particle_mobility_transport.py:184)。
