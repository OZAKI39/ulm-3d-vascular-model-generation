"""Write the Chinese reading guide and English gallery from verified records."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def main():
    c=json.loads((HERE/'COMPUTE_VALIDATION.json').read_text())
    m=json.loads((HERE/'MEDIA_VALIDATION.json').read_text())
    q=json.loads((HERE/'LOCAL_VERIFICATION.json').read_text())
    assert c['all_pass'] and m['all_pass'] and q['all_pass']
    def stat(region,quantity):return next(r for r in c['statistics'] if r['region']==region and r['quantity']==quantity)
    w=stat('J1_r5um','WSS_raw');nt=stat('J1_r5um','NormalViscousTraction_abs');angle=stat('J1_r5um','TractionAngleFromTangent')
    table='\n'.join('| '+str(r['region'])+' | '+r['quantity']+' | '+r['unit']+' | '+' | '.join(f"{r[k]:.6g}" for k in ['min','mean','P5','P50','P95','max'])+' |' for r in c['statistics'])
    stage=[('01_wall_normals','1. Wall and outward normals',f"{c['node_glyph_count']} gold arrows show outward unit node normals. All arrows have the same length; original wall triangles are visible."),
        ('02_velocity_gradient','2. Near-wall velocity gradient','Color shows |G n|, the magnitude of the normal derivative of velocity, in s⁻¹. This is not a single tensor component or WSS.'),
        ('03_viscous_traction','3. Viscous traction','Gold arrows show mu (G + G^T) n in Pa, including the normal viscous component. No pressure term is included.'),
        ('04_tangential_wss','4. Tangential WSS vectors','Cyan arrows show the tangential projection, at the same facet locations and linear scale as stage 3. This is wall-on-fluid traction.'),
        ('05_wss_magnitude','5. WSS magnitude and vectors','The existing continuous WSS color display (0–55 Pa) is combined with black raw facet vectors. Nodal color averaging does not change the saved raw values.')]
    links='\n'.join(f'| {title} | [MP4](animations/J1_{key}.mp4) | [4K PNG](figures/{key}_overview_4k.png) · [PDF](figures/{key}_overview_4k.pdf) |' for key,title,_ in stage)
    text=rf'''# Junction 1：WSS 五步计算过程动画说明

本目录基于**当前 `{c['case']}` 新流场**的真实速度与壁面网格生成。所有计算是读取现有解后的轻量后处理；没有启动 CFD，没有修改求解器、边界条件、网格或原始结果。

[打开五段动画预览](OPEN_RESULTS.html)

| 步骤 | 动画 | 静态图 |
|---|---|---|
{links}

每段 1920×1080、24 fps、432 帧、18 秒。18 秒是稳态场的相机旋转时间，不是血流物理时间。五段保持同一局部视野、相机路径、固定 Z 轴和 6° 仰角，复用现有黑色背景、turbo 配色、字体、坐标网格、旁置箭头及周期渐变标签。服务器使用 RTX 4090/EGL 渲染，CPU libx264 编码。界面文字全部英文。

## 位置、区域及数据身份

- J1 中心：**(92, 49, 111) μm**。
- 显示范围：原始壁面三角形中心距 J1 小于 **11 μm**，共 **{c['local_facets']} 面、{c['local_nodes']} 点**。不切割或重建三角形，不添加端盖。窗口边缘的开口是显示裁剪边界，不是新设置的血管出口。
- 另列 **5 μm** 球区统计，不能将它与 11 μm 显示区混用。
- 仅选择已确认 WALL tag=1；隐藏内部流体单元，入口、出口人工截面不参与 WSS。
- Toward inlet / Toward O2 / Toward J2 的引线落在对应支路的真实壁面；J1 引线指向分叉中心。分支方向来自已核实 ROI 图路径，不按边界标签大小猜测。`data/region_and_annotations.json` 保存路径采样点、引线落点、当前四端口坐标及几何核对误差。历史资料仅提供几何定位，不复用其中旧流场或 WSS。
- 输入 `../input_data/frozen_flow/steady_flow_mean_2p0_mmps.vtu`，其 SHA256 为 `{c['source_field_sha256']}`；输入别名对应实际当前算例，以源计算记录和哈希为准。
- `audit/source_lock.json` 锁定输入与父目录渲染代码。源 WSS 文件为 `../input_data/field_diagnostics/data/wall_wss_si.vtp`。

## 五步分别表示什么

### 1．壁面及向外单位法向

先从体网格找到每个边界三角形唯一的相邻四面体。三角形边向量叉积归一化后，用“相邻四面体中心 → 面中心”检验符号，保证法向指向流体域外部。

主动画按要求显示**节点法向**：将完整壁面上该节点相邻三角形的面积加权法向求和后归一化。先在完整壁面生成，再截取 J1，避免显示区边缘产生人为法向偏差。均匀空间抽样 {c['node_glyph_count']} 个节点，金色箭头长度均为 {c['normal_glyph_length_um']} μm；长度不代表应力。平滑曲面节点通常没有唯一的平面法向，因此它代表相邻面片的平均方向。

**后续原始 WSS 使用每个平面三角形自己的单位法向**，不是这些节点平均法向。这与项目实际实现一致。两个数组分别为 `NodeOutwardNormal`（point data）和 `Outward_normal`（cell data）。

### 2．完整速度梯度与壁面映射

用相邻四面体四个节点的三维速度对 P1 有限元形函数求导，直接得到 9 个分量：

$$G_{{ij}}=\frac{{\partial u_i}}{{\partial x_j}},\qquad \mathbf G=\nabla\mathbf u.$$

复用生产函数 `audit/wss_core.py::p1_gradients`，不对速度模长求梯度、不在无滑移壁面的零速度之间做表面差分。P1 单元内梯度恒定，所以每个壁面三角形直接取其唯一相邻四面体的梯度，不需要跨体单元任意插值。

图中选取的可读标量为：

$$q_n=\|\mathbf G\mathbf n\|=\left\|\frac{{\partial\mathbf u}}{{\partial n}}\right\|\quad[\mathrm{{s}}^{{-1}}].$$

它是**沿壁面法向的速度变化率向量的模长**，不是矩阵的主特征值，不是单一坐标分量，也不是上一轮动画的等效剪切率 $\sqrt{{2\mathbf D:\mathbf D}}$。9 个梯度分量全部保留在 VTP 和 CSV；`VelocityGradient_s_inv` 按 xx,xy,xz,yx,yy,yz,zx,zy,zz 排列。

色标 0–16000 s⁻¹。云图沿用完整壁面原始标量的面积加权节点显示，数据仍保存原始 cell 值；不通过显示裁剪改变数值。

### 3．黏性牵引力

实际动力黏度 **μ={c['mu_Pa_s']} Pa·s**，沿用当前恒黏度牛顿流体模型：

$$\mathbf t_\mu=\mu(\mathbf G+\mathbf G^T)\mathbf n\quad[\mathrm{{Pa}}].$$

金色箭头只表示**黏性牵引**。式中没有压力，不能把它的法向分量称为“压力干扰”。总牵引是 $-p\mathbf n+\mathbf t_\mu$，但压力部分在切向投影中消失。

当前 J1 半径 5 μm 内，$|\mathbf t_\mu\cdot\mathbf n|$ 的面积加权平均为 **{nt['mean']:.6f} Pa**，最大 **{nt['max']:.6f} Pa**；牵引偏离局部切平面的角度 P95={angle['P95']:.6f}°，最大 {angle['max']:.6f}°。因此不能预设该分量处处“微小”。

本次还验证了当前 P1 无滑移壁面上的代数关系 $\mathbf n\cdot\mathbf t_\mu=2\mu\nabla\cdot\mathbf u$，最大差 {c['normal_identity_max_error_Pa']:.3g} Pa。这说明该法向分量与当前离散速度的局部散度一致，不能用它证明存在压力混入。理想光滑、不可压缩、静止无滑移壁面下该分量应为零；本轮没有重新评估速度离散误差或网格收敛。

### 4．切向投影后的 WSS 矢量

$$\boldsymbol\tau_w=\mathbf t_\mu-(\mathbf t_\mu\cdot\mathbf n)\mathbf n.$$

青色箭头为投影后矢量。全壁面最大 $|\boldsymbol\tau_w\cdot\mathbf n|$ 为 **{c['tangency_max_error_Pa']:.3g} Pa**，接近浮点精度。

步骤 3、4、5 使用同一组 **{c['face_glyph_count']} 个面片**，最小空间间距 {c['glyph_spacing_um']} μm；抽样不依赖应力大小。三段使用相同线性比例 **{c['traction_glyph_scale_um_per_Pa']:.8f} μm/Pa**，不逐箭头归一化、不放大低值，也不隐藏反向矢量。为避免与壁面深度重叠，箭头起点沿局部法向外移 **{c['glyph_offset_um']} μm**；只影响图形位置，不修改矢量。

箭头方向与**起点所在三角形平面**相切，不代表有限长度的直箭头能沿弯曲曲面贴合整段路径；遮挡由真实 3D 深度关系决定。

**方向约定必须注意：**本项目用流体向外法向及 $+\mu(\mathbf G+\mathbf G^T)\mathbf n$，这里显示的是壁面对流体的切向牵引。流体对壁面的作用是其反向 $-\boldsymbol\tau_w$，在 VTP 中另存 `FluidOnWallShear_Pa`。本次没有为“看起来顺流”而翻转项目原始矢量。解析简单剪切检验及当前无滑移 P1 邻壁单元核对均表明，此处原始矢量与邻壁切向速度方向相反。因此不能把这些箭头直接当作流线，也不能仅凭其指向宣布存在回流、分离或真实分流比例。

### 5．WSS 模长云图与黑色矢量叠加

$$\mathrm{{WSS}}=\|\boldsymbol\tau_w\|\quad[\mathrm{{Pa}}].$$

主动画直接复用现有 `WSS_display_Pa` 节点字段，保持连续云图的现有显示方式；它是相邻原始面片 **WSS 模长的面积加权平均**，不是平均矢量的模长。黑色箭头仍来自原始面片的 `WSSVector_Pa`。色标保持 **0–55 Pa**、线性。

另存 `figures/05_raw_facet_comparison_000_4k.png` 和 `_180_4k.png`，在同视角、同色标下直接显示 `WSS_raw_Pa` 面片原始值，不进行节点平均。二者区别属于显示方式，不能把连续云图的平滑程度当成数值精度证据。

当前 5 μm J1 区域原始 WSS **{w['min']:.6f}–{w['max']:.6f} Pa**，面积加权均值 **{w['mean']:.6f} Pa**。没有为符合“分叉脊必红、外侧必蓝”的预期选点或改色；蓝色表示在统一色标下较低，不等于数值为零。本轮不据此诊断压力边界条件、回流或 WSS 精度。

## 原始数值摘要

以下均值及分位数均按原始三角形面积加权。原始面片数值没有截断、平滑或归一化。

| 区域 | 量 | 单位 | 最小 | 均值 | P5 | P50 | P95 | 最大 |
|---|---|---|---:|---:|---:|---:|---:|---:|
{table}

## 已执行核验

- 生产梯度/牵引链路重建 WSS，与现有原始 WSS 最大差 **{c['WSS_reproduction_max_error_Pa']:.3g} Pa**。
- 独立使用重心形函数的叉积公式求梯度，与生产线性方程求解方法核对；局部全部单元梯度最大差 **{q['independent_gradient_max_abs_error_s_inv']:.3g} s⁻¹**，WSS 最大差 **{q['independent_WSS_max_abs_error_Pa']:.3g} Pa**。
- 单位节点法向、面片向外法向、唯一体单元对应、仅 WALL 标签、节点显示值继承、数值单位、投影切向性及简单剪切符号检验均通过。
- 五段各 432 帧在服务器与本地完整解码；相机范围、Z 轴、标签净空、视频尺寸/帧率和 4K 图片尺寸通过检查。布局初次失败的记录保留在 `logs/attempt1_layout.err`；修正局部坐标标题候选位置后重新执行，未放宽净空检查。
- 这验证的是**后处理与可视化一致性**，不等于新的 CFD 验证或网格无关性证明。

## 文件与复现

- `data/J1_wss_pipeline_si.vtp`：真实局部壁面；坐标 m；完整梯度、节点/面片法向、牵引、投影矢量、原始与显示 WSS、父四面体 ID。
- `data/J1_facet_values.csv`：每个显示面片的坐标（μm）、面积（μm²）、法向、9 个梯度分量、牵引及 WSS。
- `data/J1_summary.csv`：5 μm 与 11 μm 两个区域的面积加权统计。
- `data/glyph_samples.npz`：三个应力阶段共享的面片抽样 ID，以及节点法向抽样 ID。
- `J1_camera.json`、`J1_annotations.json`：全部 432 帧的相机和标签。
- `COMPUTE_VALIDATION.json`、`MEDIA_VALIDATION.json`、`LOCAL_VERIFICATION.json`：数据与媒体核验。
- `prepare_pipeline.py`、`render_pipeline.py`、`verify_pipeline.py`：数据派生、服务器绘图、本地独立核验。

本地核验命令：
```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -B {HERE}/verify_pipeline.py
```

服务器目录：`/workspace/roi_only_flow_visualization_20260928/j1_wss_pipeline/`。渲染使用该目录的私有 `supervisord.conf`，仅包含 `j1_pipeline`，autostart=false、autorestart=false。复现应使用新的输出目录：
```bash
/root/particle8_2_runs/env/bin/python -B render_pipeline.py --output /workspace/roi_only_flow_visualization_20260928/j1_wss_pipeline/rerender
```
脚本拒绝覆盖已完成 MP4；服务器长期渲染通过 supervisor 管理。`prepare_pipeline.py` 同样拒绝覆盖已完成数据。

关于牵引的法向/切向分解可对照 [COMSOL 的应力与运动方程说明](https://www.comsol.com/multiphysics/stress-and-equations-of-motion)；不同软件的法向符号约定不能混用，例如 [OpenFOAM wallShearStress 文档](https://doc.openfoam.com/2312/tools/post-processing/function-objects/field/wallShearStress/)明确说明其法向约定。本项目实际依据的是上面记录的现有 P1 FEM 代码，并未使用 OpenFOAM 求解。
'''
    (HERE/'README_ZH.md').write_text(text)
    sections='\n'.join(f'''<section><h2>{title}</h2><p>{desc}</p>
<video controls loop preload="metadata" poster="figures/{key}_overview_4k.png" src="animations/J1_{key}.mp4"></video>
<p><a href="animations/J1_{key}.mp4">MP4</a> · <a href="figures/{key}_overview_4k.png">4K PNG</a> · <a href="figures/{key}_overview_4k.pdf">PDF</a> · <a href="figures/{key}_rotation_views.png">Four views</a></p></section>''' for key,title,desc in stage)
    (HERE/'OPEN_RESULTS.html').write_text('''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Junction 1 | WSS calculation stages</title><style>body{max-width:1280px;margin:32px auto;padding:0 20px;background:#0b0d11;color:#edf2f7;font:17px/1.6 system-ui}video,img{display:block;width:100%;background:black}a{color:#90caff}section{margin:38px 0}h1{font-size:28px}h2{font-size:23px}</style>
<h1>Junction 1 | WSS calculation stages</h1><p>Current ROI-only balanced flow · Real wall surface · Same camera in all five stages · 1080p / 24 fps / 18 s each</p>
<p>J1 = (92, 49, 111) µm. Display: wall facet centroids within 11 µm. Crop openings are display boundaries, not vessel outlets. The field is steady; only the camera rotates.</p>
<p>Traction arrows use the outward-fluid normal and the existing wall-on-fluid sign convention. They are not flow streamlines or proof of recirculation.</p>
<p><a href="README_ZH.md">Chinese reading guide and formulas</a> · <a href="data/J1_wss_pipeline_si.vtp">Scientific VTP data</a> · <a href="data/J1_summary.csv">Region statistics CSV</a></p>
'''+sections+'''<section><h2>Raw facet comparison</h2><p>Same viewpoint and 0–55 Pa scale; no nodal averaging in these two images.</p>
<img src="figures/05_raw_facet_comparison_000_4k.png" alt="Raw facet WSS at the initial angle"><p><a href="figures/05_raw_facet_comparison_180_4k.png">Opposite view</a></p></section></html>''')
    print('READING_GUIDE_AND_GALLERY_READY')

if __name__=='__main__':main()
