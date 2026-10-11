# 

# ① 背景血流：先知道血是怎么流的

**目的**：碎块掉下来以后要靠血流推着走，所以先算好血管里每个位置的流速。

**公式**：不可压缩流体的稳态 Navier–Stokes 方程

$\rho(\mathbf u\cdot\nabla)\mathbf u=-\nabla p+\mu\nabla^2\mathbf u,\qquad \nabla\cdot\mathbf u=0$

左边是流体惯性；右边第一项是压力推动，第二项是黏性（流体内部的“摩擦”）。$\nabla\cdot\mathbf u=0$ 的意思是：流进来多少，就流出去多少。

⚠️ 这是**没有血栓**时的流场。血栓挡住血流后流场会变，这一点还没算。

入口给的是**抛物线流速**（Poiseuille 流），中间快、贴壁处为 0。管中心约 10 mm/s。管壁上流速为 0（无滑移），两个出口压力都设为 0。雷诺数 $Re=\rho\bar u D/\mu\approx13$，远小于 2000，是平稳的**层流**。**参数 **ρ = 1000 kg/m³，μ = 0\.001 Pa·s（这是**水**的参数，不是血）；网格有 26,712 个节点、146,230 个四面体。

$u(r)=2\bar u\Big(1-\frac{r^2}{R^2}\Big),\qquad \bar u=5\ \text{mm/s}$

> **出处**：
> 
> Navier–Stokes 方程和 Hagen–Poiseuille 流都是流体力学教科书内容（如 Batchelor 1967, *An Introduction to Fluid Dynamics*）；
> 
> 求解器是 SimVascular（Updegrove et al\. 2017, *Ann\. Biomed\. Eng\.* 45:525–541）。
> 
> 

---

# ② 血栓形状：在管壁上捏一个“三个鼓包”的血块

**目的**：做一个不规则、贴在母管内壁上的血块。它是人为设计的形状，不是真实病人的数据。

**公式**：在管壁上用两个坐标：沿血管方向的 $a$，绕管壁一圈方向的弧长 $b$。血块高度为

$h(a,b)=\frac{1}{N_h}\underbrace{\max(1-q^2,0)^{1.15}}_{\text{底面轮廓}}\Big[\underbrace{h_p(1-0.09\tanh 2u)}_{\text{平台}}+h_l\underbrace{\sum_{j=1}^{3}w_j\,e^{-\frac{(a-a_j)^2+(b-b_j)^2}{2\sigma_j^2}}}_{\text{三个鼓包}}\Big]$

底面轮廓是一个略微扭曲的椭圆，边缘处高度降到 0；上面先铺一层“平台”（上游稍厚），再叠三个钟形鼓包。这是项目自己设计的工程形状，没有文献依据。

---

# ③ 切成点、连上弹簧

**目的**：计算机算不了“一整块连续的肉”，只能算有限个点，所以要**把血块切成小点，再用弹簧连起来**。

**做法**：

1. **铺点阵**：在整个空间铺一层间距 Δx = 50 µm 的立方点阵，只保留落在血块内部的点。

2. **保守裁剪**：每个点代表一个 50 µm 的小方块。方块上 27 个采样点都必须在血管腔里，否则去掉，这一步去掉了 150 个贴壁的点。只有一个邻居的“毛刺点”也去掉（1 个）。

3. **锚点**：管壁外侧再铺 3 层固定不动的点（504 / 541 / 516 个，共 1,561 个），代表血管壁。

4. **连弹簧**：两点距离不超过 δ 就连一根弹簧。$\text{点 } i,j \text{ 相连}\iff\|\mathbf X_j-\mathbf X_i\|\le\delta,\qquad \delta=3.01\,\Delta x=150.5\ \mu\text{m}$

$V_i=(\Delta x)^3=1.25\times10^{-13}\ \text{m}^3,\qquad m_i=\rho V_i=1.25\times10^{-10}\ \text{kg}$

> **出处**
> 
> 用点阵离散近场动力学、δ 取约 3 倍点距，是 Silling \& Askari \(2005\) 的常用做法（*Computers \& Structures* 83:1526–1535）。
> 
> 

---

# ④ 弹簧多硬：从实验数据换算 K 和 G

**目的**：让整个点阵的软硬程度像真实血块。

**来源**：Liu et al\. \(2024\) 用 **Yeoh 超弹性模型**拟合了血块实验：

$W=C_1(\bar I_1-3)+C_2(\bar I_1-3)^2+C_3(\bar I_1-3)^3+\frac{1}{D}(J-1)^2$

其中 $C_1=470$ Pa、$C_2=170$ Pa、$C_3=-12$ Pa、$D=6\times10^{-5}\ \text{Pa}^{-1}$。

剪切变形 γ 时，$\bar I_1-3=\gamma^2$，所以 $W\approx C_1\gamma^2$。弹性理论里剪切能量是 $\tfrac12G\gamma^2$，两者相等，得到 $G=2C_1=940$ Pa。

纯体积变化时，$W=(J-1)^2/D$，而体积能量是 $\tfrac12K(J-1)^2$，所以 $K=2/D=33{,}333$ Pa。

这两个数交给 LPS 模型（之前讲过：K 管“胀”，G 管“歪”）。

**代价**：$C_2$、$C_3$ 描述的是“越拉越硬”，这一部分被丢掉了；Liu 模型里的损伤软化（Mullins 效应）也没有保留。

> **出处**：
> 
> Yeoh \(1993\), *Rubber Chem\. Technol\.* 66:754–771；
> 
> Liu et al\. \(2024\), *J\. Mech\. Phys\. Solids* 193:105858（表 S2）；
> 
> LPS：Silling et al\. \(2007\), *J\. Elasticity* 88:151–184。
> 
> 

---

# ⑤ 弹簧拉多长会断：用“二分法”调门槛

**断裂规则**：$s_{ij}=\frac{\|\mathbf x_j-\mathbf x_i\|-\|\mathbf X_j-\mathbf X_i\|}{\|\mathbf X_j-\mathbf X_i\|},\qquad s_{ij}>s_{00}\ \Rightarrow\ \text{永久断开}$

**门槛 s₀₀ 怎么定**：做一个小方块（8×8×8 个点）的数值拉伸实验，量出它最大能承受多大应力；调 s₀₀，直到这个强度落进 Liu 2024 的实验范围。

> **出处**：
> 
> 临界伸长断裂：Silling \& Askari \(2005\)；
> 
> 二分法：介值定理，数值分析教科书（如 Burden \& Faires, *Numerical Analysis*）。
> 
> 

文献里还有一种按**断裂能**定门槛的公式（Silling \& Askari 2005，键型模型）$s_0=\sqrt{5G_c/(9K\delta)}$。

反过来算，$s_{00}=0.15$ 对应的断裂能只有约 **0\.2 J/m²**，而 Liu 2024 测的血块是 **6–18 J/m²**，差了几十倍。也就是说，现在的模型**强度对了，但非常脆**：断开时吸收的能量远少于真实血块。这可能是它碎成小碎屑、而不是断成大块的原因之一。这个公式严格说是给键型模型的，用在 LPS 上只能算粗估。

---

# ⑥ 外加剪切：模拟微泡对血块表面的“刮”

**完整公式**（作用在表面点 $i$上）：$\mathbf F_i^{\text{load}}(t)=A_i\;\tau_0\;g(t)\;\exp\!\Big(-\frac{\|\mathbf x_i(t)-\mathbf c\|^2}{2\sigma^2}\Big)\;(\mathbf I-\mathbf n_i\mathbf n_i^{\top})\,\mathbf e_x$

⚠️此处面积 $A_i$ 和方向 $\mathbf n_i$ 都用的是**最初的形状**；血块裂开后新露出来的表面不受力。

$(\mathbf I-\mathbf n_i\mathbf n_i^{\top})\mathbf e_x$：把顺着血流的方向 $\mathbf e_x$ 压到表面上，变成贴着表面刮的方向。正对血流的那部分表面，分到的力就小。

钟形部分：中心 $\mathbf c$ 在血块上游一侧，σ = 0\.28 mm；按点**当前位置**计算。

$\tau_0=780$ Pa，由预算反推：让整个表面的**平均**剪切等于 200 Pa，最强的点约 771 Pa。

$\tau_0=\frac{200\ \text{Pa}\times\sum_iA_i}{\sum_iA_i\,\|\text{形状}_i\|}$

$g(t)$ 是一个平滑脉冲：前 1/4 时间慢慢升到 1，中间 1/2 保持，最后 1/4 慢慢降到 0。

$g(t)=\begin{cases}\tfrac12[1-\cos(4\pi t/T)] & 0\le t<T/4\\ 1 & T/4\le t<3T/4\\ \tfrac12[1+\cos(4\pi(t/T-\tfrac34))] & 3T/4\le t<T\\ 0 & t\ge T\end{cases}$

**脉冲时长 T = 1\.16 ms 的依据**：剪切波穿过血块高度所需的时间

$T_{\text{cross}}=\frac{H}{\sqrt{G/\rho}}=\frac{0.65\ \text{mm}}{0.97\ \text{m/s}}\approx0.67\ \text{ms},\qquad T\approx1.73\,T_{\text{cross}}$

意思是让力作用得够久，整个血块都能“感受到”它，而不是只有表面抖一下。

---

# ⑦ 逐步求解：每 0\.5 微秒算一次所有点怎么动

**总方程**（每个血块点）：

$m_i\ddot{\mathbf x}_i=\underbrace{\mathbf F_i^{\text{PD}}}_{\text{弹簧}}+\underbrace{\mathbf F_i^{\text{contact}}}_{\text{点碰点}}+\underbrace{\mathbf F_i^{\text{wall}}}_{\text{点碰管壁}}+\underbrace{\mathbf F_i^{\text{load}}}_{\text{外加剪切}}+\begin{cases}-500\,m_i\dot{\mathbf x}_i & \text{还附着}\\ \mathbf F_i^{\text{flow}} & \text{已脱落}\end{cases}$

锚点始终固定不动。

**\(a\) 点碰点**（弹簧断开后，两个点靠得太近就互相推开）

$F=15\,k\,\bar V\,\frac{d-r}{\delta}\quad(r<d),\qquad k=\frac{18K}{\pi\delta^4},\qquad d=1.0\,\Delta x$

- $k$ 是键型模型的标准“微模量”（Silling \& Askari 2005）；系数 15 是 LAMMPS 沿用的 EMU 手册约定。

- 项目**改了 LAMMPS 源码**，把 $d$ 从默认的 1\.35 Δx 改成 1\.0 Δx。原因前面讲过：斜对角邻居本来就只隔约 1\.41 Δx，1\.35 太接近，会引发连锁粉碎。



**\(b\) 点碰管壁**（弹簧加阻尼器，只在靠近管壁时起作用）

$\mathbf F^{\text{wall}}=\big[k_w\,(r_p-d)_++c\,(-\mathbf v\cdot\mathbf n)_+\big]\mathbf n,\qquad c=2\zeta\sqrt{k_w\,m}$

$k_w=50$ N/m，$r_p=20\ \mu$m，$\zeta=0.5$。$(\cdot)_+$ 表示负数取 0：离壁超过 20 µm 不受力；离开管壁时不加阻尼，只在撞向管壁时减速。



**\(c\) 时间推进：速度 Verlet 法**（LAMMPS 的 `fix nve`）

$\mathbf v^{n+\frac12}=\mathbf v^n+\frac{\Delta t}{2m}\mathbf F^n,\quad \mathbf x^{n+1}=\mathbf x^n+\Delta t\,\mathbf v^{n+\frac12},\quad \mathbf v^{n+1}=\mathbf v^{n+\frac12}+\frac{\Delta t}{2m}\mathbf F^{n+1}$

$\Delta t=0.5\ \mu$s，共 4,660 步，即 2\.33 ms。对弹簧振子，Verlet 法稳定的条件是 $\omega\Delta t<2$。管壁弹簧 $\omega=\sqrt{k_w/m}=6.3\times10^5\ \text{s}^{-1}$，$\omega\Delta t=0.32$，满足。材料里纵波速度约 5\.9 m/s，穿过一个点距要 8\.5 µs，也远大于 0\.5 µs。

> **出处**：
> 
> Verlet \(1967\), *Phys\. Rev\.* 159:98；
> 
> Swope et al\. \(1982\), *J\. Chem\. Phys\.* 76:637；
> 
> 弹簧–阻尼接触：Cundall \& Strack \(1979\), *Géotechnique* 29:47–65；
> 
> 稳定性：Hairer, Lubich \& Wanner \(2006\), *Geometric Numerical Integration*。
> 
> 

---

# ⑧ 本次新增：每一步检查谁掉了，掉了立刻受血流拖

**\(a\) 怎么判断“掉了”**

只看还没断的血块内部弹簧，把点分成一团一团**互相连着的团**。一个团里只要还有一根界面弹簧连着管壁，就算附着；一根都没有，就算脱落。

$\text{团 }\mathcal C\text{ 附着}\iff\mathcal C\text{ 里至少有一个点，还有一根没断的界面弹簧}$

程序每一步都看 LAMMPS 里弹簧的总数，一变就重新分团，所以**掉下来的那一步就知道**。最早的一次脱落发生在 0\.2965 ms（节点 267）。

> 出处：
> 
> 连通分量算法（Hopcroft \& Tarjan 1973, *Commun\. ACM* 16:372–378），代码用 SciPy 实现。
> 
> 



**\(b\) 血流拖力：按“团”整体算**

先求团的质心位置 $\mathbf x_c$ 和质心速度 $\mathbf v_c$，再在质心处查血流速度 $\mathbf u$：

$\mathbf F=\frac{M(\mathbf u-\mathbf v_c)}{\tau},\qquad \tau=\begin{cases}\dfrac{\rho\,d_{\text{eq}}^2}{18\,\mu\,f(Re)} & \text{团里} \ge 3 \text{ 个点}\\[6pt] 100\ \mu\text{s} & 1\text{–}2\text{ 个点（人为设定）}\end{cases}$$f(Re)=1+0.15\,Re^{0.687},\qquad Re=\frac{\rho\,d_{\text{eq}}\|\mathbf u-\mathbf v_c\|}{\mu}$

- **τ 的含义**：碎块追上血流速度大约需要的时间。例如一个 4 个点的碎块，$d_{\text{eq}}\approx98\ \mu$m，τ ≈ 0\.5 ms。

- **1–2 个点的碎屑**分辨率太低，不当成真实颗粒，直接规定 τ = 100 µs。



**\(c\) 查血流速度：线性插值**

质心落在流场网格的某个四面体里，就用四个顶点的速度按距离加权平均：如果质心落在网格外面，不外推，这一步不给力。



**\(d\) 把团的总力分给每个点**

$\mathbf F_i=\frac{m_i}{M}\mathbf F$

这样分有两个好处。

1. 总力不变：$\sum_i\mathbf F_i=\mathbf F$。

2. 不会凭空让碎块转起来。绕质心的力矩为$\sum_i(\mathbf x_i-\mathbf x_c)\times\mathbf F_i=\frac1M\Big[\sum_im_i(\mathbf x_i-\mathbf x_c)\Big]\times\mathbf F=\mathbf 0$    



**\(e\) 阻尼只给还附着的点**：−500 m v 是让计算稳定的“数值阻尼”，不是真实物理。脱落的点不再加它，改受血流拖力。旧版到 2\.33 ms 才打开血流拖力，之前碎屑以中位约 0\.86 m/s 的速度乱飞（而当地血流只有约 5 mm/s），然后突然被拉慢，所以视频看起来“先停后动”，还有 79 个点撞上了管壁。

新版碎块的速度和血流基本一致，代价是计算时间变成约 3 倍，因为每一步都要调用一次 Python。

> **出处**：
> 
> Stokes \(1851\) 阻力；
> 
> Schiller \& Naumann \(1933\), *Z\. Ver\. Dtsch\. Ing\.* 77:318–320；
> 
> 质心定理是经典力学内容（如 Goldstein, *Classical Mechanics*）；
> 
> 线性插值是有限元标准内容。
> 
> 

---

# ⑨ 交接：2\.33 ms 之后，跟踪碎块走 5 秒

- **目的**：LAMMPS 每步只走 0\.5 µs，要算 5 秒得一千万步，太贵。所以 2\.33 ms 之后换一个便宜的方法：每个碎块当成一个不变形的整体，只算它的质心。

- **公式**：用同一个拖力公式。在一小步里假设血流速度 $\mathbf u$ 不变，方程 $\dot{\mathbf v}=(\mathbf u-\mathbf v)/\tau$ 有**精确解**：

$\mathbf v(t+\Delta t)=\mathbf u+(\mathbf v-\mathbf u)\,e^{-\Delta t/\tau}$$\mathbf x(t+\Delta t)=\mathbf x+\mathbf u\,\Delta t+(\mathbf v-\mathbf u)\,\tau\big(1-e^{-\Delta t/\tau}\big)$

- **定理**：这是线性常微分方程的精确解，所以就算步长比 τ 大很多也不会算炸。代码规定每步不超过 5 ms，且一步位移不超过 20 µm。

- **管壁处理**：如果碎块里某个点要钻进管壁 5 µm 以内，就沿法向把它推回来，并去掉朝向管壁的速度（不弹、不摩擦）。

---

# 碎块统计和本次结果

碎块大小用等效球直径 $d_{\text{eq}}=(6N_k\Delta x^3/\pi)^{1/3}$ 衡量。

- 第一次内部弹簧断开：0\.1205 ms；第一次有碎块脱落：0\.2965 ms。

- 最后共有 217 个点脱落（约 11\.3%），**约 88\.7% 仍附着**。

- 208 个碎块里：201 个是单点，6 个是两点，1 个是 4 点。

- 整个过程中最大的碎块是 13 个点，没有达到 200 µm 所需的 34 个点。

所以结论是 **附着主体 \+ 表面碎屑脱落**，还没有出现大块断裂。



