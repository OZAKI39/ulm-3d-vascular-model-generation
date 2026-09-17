> From: https://chatgpt.com/c/6a54af07-a098-83ee-a4ff-ccb4ef2b52ad

# chatgpt response

## 审查结论

上一版方案保留了三个正确方向：利用红细胞碰撞造成的横向输运、允许不同分支具有不同血细胞比容，以及在微泡进入无细胞层后停止红细胞平均壁向作用。你提供的文档将它们分别表述为剪切诱导扩散、局部血细胞比容依赖和动态无细胞层截止。

但上一版仍有四处需要修正。

第一，红细胞诱导的横向扩散本质上是随机碰撞造成的群体色散，不能严格等同于作用在单颗微泡上的确定性力。Vahidkhah 和 Bagchi 发现，红细胞碰撞使微粒的壁面法向扩散比布朗扩散高数个数量级；同时，靠近无细胞层边缘的颗粒还可能发生快速、非连续的“瀑布式”迁移。  因此，正式模型可以用扩散系数确定**平均迁移时间尺度**，但不能把它描述成真实的单次红细胞碰撞，也不应加入未经标定的随机白噪声。

第二，上一版使用的分叉相分离系数不是 Pries–Secomb 2005 修正模型的完整形式。若采用该经典网络模型，应使用母血管当前的排出血细胞比容和 2005 年修正后的 $A$、$B$、$X_0$，而不能始终使用根血管血细胞比容，也不能混用旧版系数。Pries 模型的目的正是同时描述 Fåhræus 效应和分叉处红细胞—血浆不同比例分流。([PubMed](https://pubmed.ncbi.nlm.nih.gov/2208609/?utm_source=chatgpt.com))

第三，不能由

$$
H_t/H_D
$$

直接假设红细胞核心浓度等于 $H_D$，再反推出无细胞层厚度。$H_D$ 是红细胞通量意义上的排出血细胞比容，并不等于红细胞核心区的局部体积分数。该推导虽然量纲正确，却包含未经验证的核心浓度假设。更可靠的做法是使用小鼠皮层血管的体内 THG/3PEF 测量结果。论文 [Label-free assessment of hemodynamics in individual cortical brain vessels using third harmonic generation microscopy](<../references/Label-free assessment of hemodynamics in individual cortical brain vessels using third harmonic generation microscopy.pdf>) 将单侧无细胞层厚度定义为完整管腔直径与红细胞占据区直径之差的一半。结果显示，CFL 在毛细血管中随直径快速增加，在约 $8\,\mu\mathrm m$ 直径处达到约 $1\,\mu\mathrm m$，并在更大的小动脉和小静脉中保持在约 $1\,\mu\mathrm m$。

第四，将横向扩散系数进一步乘以未经验证的微泡碰撞截面积、壳层刚度或经验常数，会重新引入不可辨识参数。当前最稳妥的闭合只在与参考研究相近的中等血细胞比容、微米级近球形载体和多列红细胞流动条件下使用，不额外拟合尺寸指数。

经过上述修正，最终方案只新增一个生理输入：

$$
\boxed{H_{D,0}}
$$

即根血管的排出血细胞比容。其余量均由血管拓扑、流量、连续几何、现有流场和固定文献关系确定。

---

## 局部血细胞比容

设母血管为 $p$，两条子血管为 $1$ 和 $2$，流量满足

$$
Q_p=Q_1+Q_2.
$$

定义流入第一条子血管的血液流量比例

$$
\boxed{
F_Q=\frac{Q_1}{Q_p}.
}
$$

设母血管排出血细胞比容为 $H_{D,p}$，母、子血管直径分别为 $D_p,D_1,D_2$，单位均为 $\mu\mathrm m$。Pries–Secomb 相分离模型定义

$$
\boxed{
A
=
-13.29
\frac{1-H_{D,p}}{D_p}
\frac{
(D_1/D_2)^2-1
}{
(D_1/D_2)^2+1
},
}
$$

$$
\boxed{
B
=
1+
6.98
\frac{1-H_{D,p}}{D_p},
}
$$

以及

$$
\boxed{
X_0
=
0.964
\frac{1-H_{D,p}}{D_p}.
}
$$

在

$$
X_0<F_Q<1-X_0
$$

时，定义

$$
z
=
\frac{F_Q-X_0}{1-2X_0},
$$

红细胞通量进入第一条子血管的比例为

$$
\boxed{
F_E
=
\frac{
1
}{
1+
\exp
\left[
-A-B\ln\left(\dfrac{z}{1-z}\right)
\right]
}.
}
$$

在极低流量分支中采用模型本身的截断：

$$
\boxed{
F_E=
\begin{cases}
0,
&
F_Q\leq X_0,\\[1mm]
1,
&
F_Q\geq1-X_0.
\end{cases}
}
$$

这些表达式来自 Pries–Secomb 微血管相分离模型的修正版；该模型说明，小于约 $30\,\mu\mathrm m$ 的血管中，红细胞通量分配可能明显偏离血液体积流量分配。([imagwiki.nibib.nih.gov](https://www.imagwiki.nibib.nih.gov/physiome/jsim/models/webmodel/NSR/phaseseparation))

两条子血管的排出血细胞比容为

$$
\boxed{
H_{D,1}
=
H_{D,p}
\frac{F_E}{F_Q},
}
$$

$$
\boxed{
H_{D,2}
=
H_{D,p}
\frac{1-F_E}{1-F_Q}.
}
$$

由此严格满足红细胞通量守恒：

$$
\boxed{
Q_pH_{D,p}
=
Q_1H_{D,1}
+
Q_2H_{D,2}.
}
$$

这里必须使用每一级母血管自身的

$$
H_{D,p},
$$

不能在整个血管树中反复使用同一个根血管 $H_{D,0}$。

微泡实际经历的是血管段内部的红细胞体积分数，因此还需要把排出血细胞比容 $H_D$ 转换为管内血细胞比容 $H_t$。采用 Pries 的 Fåhræus 关系：

$$
\boxed{
\frac{H_t}{H_D}
=
H_D+
(1-H_D)
\left[
1+
1.7e^{-0.415D}
-
0.6e^{-0.011D}
\right],
}
$$

即

$$
\boxed{
H_t
=
H_D
\left\{
H_D+
(1-H_D)
\left[
1+
1.7e^{-0.415D}
-
0.6e^{-0.011D}
\right]
\right\}.
}
$$

该关系描述红细胞在小血管中平均运动快于血浆，使管内血细胞比容通常低于排出血细胞比容。([imagwiki.nibib.nih.gov](https://www.imagwiki.nibib.nih.gov/physiome/jsim/models/webmodel/NSR/fahraeuseffect?utm_source=chatgpt.com))

如果相分离方程产生

$$
H_D<0
\quad\text{或}\quad
H_D>1,
$$

不应通过强行截断来维持运行，而应判定该分叉的流量、直径或模型适用范围不一致。

---

## 无细胞层（CFL）与微泡的近壁位置

### 先理解 CFL

红细胞在细小血管中流动时，大多集中在血管中间。靠近血管壁的位置会形成一圈红细胞较少的区域，这就是无细胞层：

$$
\mathrm{CFL}=\text{cell-free layer}.
$$

这里的“无细胞”表示细胞明显较少，并不表示任何时候都完全没有细胞。

可以把血管想成一根透明吸管。吸管的内径是完整血管直径，中间移动的一条红色带是红细胞集中经过的区域，吸管内壁与红色带之间的空隙就是 CFL。

### 论文怎样测量 CFL

[THG/3PEF 论文](<../references/Label-free assessment of hemodynamics in individual cortical brain vessels using third harmonic generation microscopy.pdf>) 使用两种图像观察同一条血管：

- 3PEF-FITC 显示血浆填满的范围，用来测量完整血管腔直径 $D_{\mathrm{3PEF-FITC}}$。
- THG 突出显示红细胞等血细胞经过的范围，由此得到 $D_{\mathrm{THG}}$。

完整血管直径减去中间血细胞区域的直径，得到左右两侧 CFL 的总厚度。因此，平均单侧 CFL 厚度为

$$
\boxed{
\delta_{\mathrm{CFL}}
=
\frac{
D_{\mathrm{3PEF-FITC}}
-
D_{\mathrm{THG}}
}{2}
}
$$ 

例如，

$$
D_{\mathrm{3PEF-FITC}}=10\,\mu\mathrm m,
\qquad
D_{\mathrm{THG}}=8\,\mu\mathrm m,
$$

则一侧 CFL 的厚度为

$$
\delta_{\mathrm{CFL}}
=
\frac{10-8}{2}\,\mu\mathrm m
=
1\,\mu\mathrm m.
$$

这两个直径都是图像信号的半高全宽，也就是在信号最高值一半的位置测得的宽度。

论文测量了 4 只小鼠中的 52 条血管。Figure 1E 显示：在较细的毛细血管中，CFL 随血管直径增加；结果部分写的是血管直径接近 $8\,\mu\mathrm m$ 时，CFL 达到约 $1\,\mu\mathrm m$；讨论部分则把进入稳定阶段概括为直径大于约 $10\,\mu\mathrm m$。这表示 CFL 大约在 $8$ 至 $10\,\mu\mathrm m$ 附近逐渐进入 $1\,\mu\mathrm m$ 左右的平台，而不是在某个精确直径突然改变。

图中的血管直径大约覆盖 $3$ 至 $50\,\mu\mathrm m$。其中两个较大小动脉出现了负的测量值。这不是说 CFL 真有负厚度，而是因为血管壁本身产生的 THG 信号干扰了边界判断。

实验使用的大分子荧光物质可能不能完全进入血管壁表面的糖萼层，因此测得的血管腔和 CFL 都可能略小。论文认为该测量可能更接近红细胞与糖萼起始位置之间的距离，而不一定包括完整糖萼。

### 程序怎样由血管直径估算 CFL

论文通过两种图像的直径差直接测量 CFL，没有给出一个只依赖血管直径 $D$ 的拟合公式。当前程序只有 $D$，没有每条血管的 $D_{\mathrm{THG}}$，所以选择结果部分明确提到的 $8\,\mu\mathrm m$ 作为平台起点，并采用简单近似：

$$
\boxed{
\delta_{\mathrm{CFL}}(D)
=
\begin{cases}
\displaystyle
1\,\mu\mathrm m
\frac{D}{8\,\mu\mathrm m},
&0<D<8\,\mu\mathrm m,\\[8pt]
1\,\mu\mathrm m,
&D\geq8\,\mu\mathrm m.
\end{cases}
}
$$

例如，

$$
\delta_{\mathrm{CFL}}(4\,\mu\mathrm m)
=0.5\,\mu\mathrm m,
$$

$$
\delta_{\mathrm{CFL}}(8\,\mu\mathrm m)
=
\delta_{\mathrm{CFL}}(50\,\mu\mathrm m)
=1\,\mu\mathrm m.
$$

这条分段关系是本项目根据论文趋势构造的计算近似，不是论文作者发表的公式。超出约 $3$ 至 $50\,\mu\mathrm m$ 的实验覆盖范围时，程序仍会计算，但会将该血管标记为超出测量范围。

CFL 还可能受到红细胞比例、流量、红细胞聚集和变形能力影响。因此，这个近似只表示程序使用的平均 CFL 厚度，不能代替每条真实血管的直接测量。

### CFL 怎样决定微泡能否完全进入近壁区域

设微泡半径为 $R_i$，那么微泡直径就是

$$
2R_i.
$$

再设 $g_i$ 为微泡表面到血管壁的距离。微泡远离血管壁的一侧到壁面的距离就是

$$
g_i+2R_i.
$$

当这一距离刚好等于 CFL 厚度时，微泡远离壁面的一侧正好接触红细胞集中区域的边缘：

$$
g_i+2R_i=\delta_{\mathrm{CFL}}.
$$

因此，程序使用的目标间隙为

$$
\boxed{
g_{\ast,i}
=
\left[
\delta_{\mathrm{CFL}}(D_i)-2R_i
\right]_+
},
$$

其中

$$
[x]_+=\max(x,0).
$$

这可以分成两种情况理解：

- 如果 $\delta_{\mathrm{CFL}}>2R_i$，CFL 比微泡直径更厚，微泡能够完全放入 CFL。
- 如果 $\delta_{\mathrm{CFL}}\leq2R_i$，CFL 不足以完整容纳微泡，此时 $g_{\ast,i}=0$。

$g_{\ast,i}=0$ 只表示微泡需要靠到壁面附近才能最大程度进入 CFL，并不允许微泡穿过血管壁；是否接触血管壁仍由已有的接触处理保证。

---

## 红细胞增强横向输运

Vahidkhah 和 Bagchi 在

$$
H_{t,\mathrm{ref}}=0.24,
\qquad
\dot\gamma_{\mathrm{ref}}
=
1000\,\mathrm{s^{-1}}
$$

条件下得到约

$$
40\!-\!50\,\mu\mathrm m^2/\mathrm s
$$

的红细胞增强壁面法向扩散系数。取其中值：

$$
\boxed{
D_{\perp,\mathrm{ref}}
=
45\,\mu\mathrm m^2/\mathrm s.
}
$$

该横向输运源于红细胞—颗粒反复碰撞，比布朗扩散强两至三个数量级；论文同时表明，碰撞频率和单次碰撞引起的横向位移共同决定边集效果。

围绕这一参考状态，采用最低阶局部缩放：

$$
\boxed{
D_{\perp,i}
=
45
\left(
\frac{H_{t,i}}{0.24}
\right)
\left(
\frac{|\dot\gamma_i|}{1000\,\mathrm{s^{-1}}}
\right)
\frac{\mu\mathrm m^2}{\mathrm s}.
}
$$

这里没有新的待拟合系数。该关系不是论文直接给出的普适经验公式，而是以碰撞频率随红细胞浓度和剪切运动增加为依据，对参考状态作线性局部延拓。因此只在

$$
0.15\lesssim H_t\lesssim0.30
$$

以及与参考剪切率同数量级的多列红细胞流中解释为定量闭合；超出该范围时不能通过调节系数强行外推。

不再额外乘以

$$
(a_R+R_i)^2
$$

或其他微泡尺寸函数。现有文献证明载体尺寸和形状会影响红细胞碰撞，但没有提供适用于含气球形微泡的统一尺寸指数。Vahidkhah 和 Bagchi 的结果还显示，在红细胞富集核心区中，不同颗粒形状的扩散系数差异相对较小，而形状差异主要影响接近无细胞层后的快速迁移和壁面接触。 因此，当前微泡半径只通过有限尺寸几何和 $g_\ast$ 进入模型，不引入未经验证的扩散尺寸缩放。

设微泡从血管中心附近移动到目标间隙 $g_{\ast,i}$ 需要跨越的代表性距离为

$$
\boxed{
\ell_i
=
\frac{D_i}{2}
-
R_i
-
g_{\ast,i}.
}
$$

该量必须满足

$$
\ell_i>0.
$$

根据

$$
\langle y^2\rangle
=
2D_\perp t,
$$

定义代表性的红细胞边集时间尺度：

$$
\boxed{
\tau_{\mathrm m,i}
=
\frac{
\ell_i^2
}{
2D_{\perp,i}
}.
}
$$

这里的 $\tau_{\mathrm m}$ 只用于把群体横向色散转换成**条件平均位置的演化时间尺度**。它不表示每颗微泡真实经历了连续、平滑的红细胞作用力，也不试图复现单次瀑布事件。

---

## 血管尺度切换

取红细胞有效直径

$$
\boxed{
d_R=8\,\mu\mathrm m,
}
$$

并定义

$$
\chi_i
=
\frac{D_i}{d_R}.
$$

Takeishi 和 Imai 发现，当

$$
\chi\leq1.25
$$

时，红细胞主要形成单列团块流，微米颗粒被相邻红细胞之间的涡旋捕获，而不是发生普通边集；当

$$
\chi\gtrsim1.5
$$

时，红细胞逐渐形成多列流，微粒才明显进入近壁贫细胞层。

因此采用没有自由陡峭度参数的平滑激活函数：

$$
\boxed{
\eta(\chi)
=
\begin{cases}
0,
&
\chi\leq1.25,\\[2mm]
3s^2-2s^3,
&
1.25<\chi<1.5,\\[2mm]
1,
&
\chi\geq1.5,
\end{cases}
}
$$

其中

$$
s
=
\frac{\chi-1.25}{0.25}.
$$

在狭窄毛细血管中，

$$
\eta=0.
$$

这不表示红细胞没有影响，而是明确表示：该区域的主要机制是团块捕获和单列红细胞输运，当前确定性边集闭合不对其进行虚假近似。Takeishi 和 Imai 还发现，团块捕获随血细胞比容降低而增强、随剪切率升高而降低，进一步说明其不能被表示成一个简单的壁向速度。

---

## 最终的确定性红细胞修正

设 $\mathbf n_i$ 为由最近连续壁面指向管腔内部的单位法向。红细胞造成的条件平均边集速度定义为

$$
\boxed{
\mathbf V_i^{\mathrm{RBC}}
=
-
\eta(\chi_i)
\frac{
[g_i-g_{\ast,i}]_+
}{
\tau_{\mathrm m,i}
}
\mathbf n_i.
}
$$

其中

$$
[x]_+
=
\max(x,0).
$$

该式具有以下确定性质：

当

$$
g_i>g_{\ast,i}
$$

时，微泡的条件平均位置逐渐向无细胞层移动。

当

$$
g_i\leq g_{\ast,i}
$$

时，

$$
\boxed{
\mathbf V_i^{\mathrm{RBC}}=\mathbf0.
}
$$

因此红细胞模型不会在微泡已经进入无细胞层后继续把它压向壁面。

当

$$
D_i/d_R\leq1.25
$$

时，

$$
\boxed{
\mathbf V_i^{\mathrm{RBC}}=\mathbf0,
}
$$

从而避免把毛细血管团块流误表示成普通边集。

在局部 $D,H_t,\dot\gamma$ 近似不变时，间隙满足

$$
\boxed{
g_i(t)-g_{\ast,i}
=
\left[
g_i(0)-g_{\ast,i}
\right]
\exp
\left[
-\eta(\chi_i)
\frac{t}{\tau_{\mathrm m,i}}
\right].
}
$$

这条解析关系只描述平均迁移趋势。真实红细胞悬浮液中的单颗微泡仍会表现出碰撞涨落和瀑布式跳跃；确定性模型不声称重现这些瞬时轨迹。

若微泡处于血管中心线等最近壁面不唯一的位置，则根据对称性规定

$$
\boxed{
\mathbf V_i^{\mathrm{RBC}}=\mathbf0.
}
$$

模型不能在完全对称状态下人为选择左壁或右壁。微泡只有在入口位置、原有流线、微泡碰撞或血管几何已经打破对称性后，才使用唯一最近壁面的法向。

包含红细胞影响后的自由广义速度为

$$
\boxed{
\mathcal U_i^{\mathrm{free}}
=
\mathcal U_i^0
+
\begin{pmatrix}
\mathbf V_i^{\mathrm{RBC}}\\
0
\end{pmatrix}.
}
$$

红细胞项只修改条件平均平移速度：

- 不修改迁移率矩阵；
- 不产生额外转矩；
- 不修改分子成键与解离参数；
- 不代替预测性不可穿透约束；
- 不直接增加分子键力。

边集、初始壁面接触和牢固黏附仍应分别处理。Vahidkhah 和 Bagchi 明确指出，这三个阶段由不同机制决定，高近壁积累并不必然意味着高稳定黏附。

---

## 模型输入、边界与验收

该模型唯一新增的生理输入为

$$
\boxed{
H_{D,0}.
}
$$

其余数值均固定为：

$$
d_R=8\,\mu\mathrm m,
$$

$$
H_{t,\mathrm{ref}}=0.24,
$$

$$
|\dot\gamma|_{\mathrm{ref}}
=
1000\,\mathrm{s^{-1}},
$$

$$
D_{\perp,\mathrm{ref}}
=
45\,\mu\mathrm m^2/\mathrm s.
$$

不再引入：

$$
K_c,\quad
\Gamma_{\mathrm m},\quad
\text{随机种子},\quad
\text{瀑布事件概率},\quad
\text{Sigmoid 陡峭度},
$$

也不根据靶区是否产生成键反向调节任何红细胞参数。

该方案是一种**单向耦合的平均边集闭合**。它要求现有 CFD 速度场已经被解释为宏观平均血流。若流场仍采用纯血浆黏度，则不能把最终结果描述为完全自洽的全血流动模拟；局部血细胞比容在本方案中只用于修正微泡横向输运，而不重新求解血液黏度和流量。

正式结果必须满足：

1. 每个分叉处红细胞通量守恒：

   $$
   Q_pH_{D,p}
   =
   Q_1H_{D,1}
   +
   Q_2H_{D,2}.
   $$
2. 在

   $$
   D/d_R\leq1.25
   $$

   时，红细胞边集修正严格为零。
3. 在

   $$
   g\leq g_\ast
   $$

   时，红细胞修正严格为零。
4. 在恒定直血管中，数值间隙应收敛到

   $$
   g(t)-g_\ast
   =
   [g(0)-g_\ast]
   e^{-\eta t/\tau_m}.
   $$
5. 当

   $$
   H_{D,0}=0
   $$

   时，应严格恢复当前无红细胞轨迹。
6. 引入红细胞影响后，首先比较靶区几何暴露量

   $$
   \boxed{
   \mathcal E_T
   =
   \sum_i
   \int_0^T
   A_{T,i}(t)\,\mathrm dt.
   }
   $$

   只有当

   $$
   \mathcal E_T>0
   $$

   后，才启动分子成键参数比较。当前基线中最小靶区间隙仍大于最大捕获距离，因此红细胞模型首先需要证明它改变了近壁暴露，而不是直接证明它提高了成键率。

最终模型可压缩为：

$$
\boxed{
\mathbf V_i^{\mathrm{RBC}}
=
-
\eta\!\left(\frac{D_i}{8\,\mu\mathrm m}\right)
\frac{
[g_i-g_{\ast,i}]_+
}{
\ell_i^2/
\left(2D_{\perp,i}\right)
}
\mathbf n_i,
}
$$

其中

$$
\boxed{
g_{\ast,i}
=
\left[
\delta_{\mathrm{CFL}}(D_i)-2R_i
\right]_+,
}
$$

$$
\boxed{
D_{\perp,i}
=
45
\left(
\frac{H_{t,i}}{0.24}
\right)
\left(
\frac{|\dot\gamma_i|}{1000\,\mathrm{s^{-1}}}
\right)
\frac{\mu\mathrm m^2}{\mathrm s}.
}
$$

该方案保留了分叉血细胞比容异质性、红细胞增强横向输运和无细胞层截止三个必要机制，同时删除了随机噪声、任意拟合系数和过多状态变量。其输出应被解释为**红细胞对球形微泡条件平均边集的确定性预测**，而不是单次红细胞碰撞轨迹的直接重建。

---

Powered by [AI Exporter](https://saveai.net)
