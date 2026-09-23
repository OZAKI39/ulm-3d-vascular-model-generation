# FEM 流线、速度矢量、压力与 WSS 旋转可视化

本目录包含独立可运行的完整渲染代码、所需输入数据副本、依赖版本和最终结果。无需引用项目其他目录中的 Python 模块，也不需要重新运行 FEM 求解或流线积分。

## 直接使用

```bash
cd /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/rotate_visualization

# 先生成两张 4K 图片及四个角度的预览
./run.sh --stills-only --output preview

# 生成完整动画、静态图、数据记录与浏览器预览页
./run.sh

# 生成放大的全血管压力与 WSS 动画、4K 图及记录
./run.sh --surface-fields

# 独立验证压力与 WSS 的数据、相机、文字及视频
./run.sh --validate-surface-fields
```

`run.sh` 会检查依赖，优先使用本目录的 `.venv`，其次使用本次实际渲染环境，再尝试上一级项目环境和 `python3`。也可以显式指定环境：

```bash
FEM_ROTATE_PYTHON=/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python ./run.sh
```

迁移到新机器后，可在本目录创建环境：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

输出默认写入 `results/`。打开 `results/OPEN_RESULTS.html` 可预览视频、下载 MP4 与 4K 图片。若需要保留已有结果，使用 `--output 新目录`。

压力与 WSS 默认写入 `results/surface_fields/`，预览页为该目录中的 `OPEN_RESULTS.html`。本目录顶层的 `OPEN_RESULTS.html` 汇总四种视图。可使用 `./run.sh --surface-fields --stills-only --output preview_surface` 先生成静态预览。

使用其他数据路径时，保持 `input_data/` 所示的输入文件结构：

```bash
./run.sh --case /path/to/mean-2p0-mmps --output /path/to/new_results
```

## 内容

- `render_visualization.py`：完整渲染、局部裁剪、速度箭头、固定 Z 轴相机、刻度与网格、文字排版、视频编码及解码验证。
- `render_surface_fields.py`：压力和 WSS 的完整渲染入口，复用本目录的固定 Z 轴、网格及渐变文字代码，不依赖其他项目模块。
- `validate_surface_fields.py`：独立核对压力原值、WSS 面积加权节点值、全表面逐帧投影边界、端口箭头及文字连续性。
- `validate_results.py`：独立验证物理坐标、速度矢量插值、固定旋转轴、缩放比例、视频文件和输入校验值。
- `inspect_transitions.py`：从最终 MP4 提取文字换位过程，生成渐变审核拼图；使用安装了本目录依赖的 Python 运行。
- `run.sh`：运行入口。
- `requirements.txt`：本次实际运行环境的 Python 依赖版本；Linux/WSL 环境使用系统 DejaVu Sans 字体。
- `input_data/`：重现本次图像所需的原始输入副本，包括完整 FEM 场、已有 96 条流线、壁面、入口出口及来源记录。
- `INPUT_MANIFEST.json`：输入来源和 SHA256。
- `results/`：最终两段动画、4K 图、相机轨迹、局部矢量数据与验证记录。
- `PACKAGING_PROVENANCE.json`：独立版本与原渲染脚本的对应关系。

## 最终显示约定

全血管保持 96 条已计算流线，并标注 `Inlet`、`Outlet 01`、`Outlet 02`、`Outlet 03`。顶部为 `Full vessel | FEM streamlines`，底部仅为 `96 streamlines`。

端口文字移至血管旁边，用细引线和三角箭头指向真实边界中心。规划时检查整个旋转周期，端口文字矩形与血管表面投影保持至少 22 像素间距（以 1920 × 1080 为基准）；淡入、淡出过程也参与检查。独立验证还检查文字矩形不与任何投影表面三角形的包围框重叠。

所有文字均直接绘制为透明叠加层，不绘制黑色背景矩形。端口标签、坐标标题在需要换位时采用 21 帧（0.875 秒）的升余弦交叉淡化；刻度文字在拥挤或接近视野边缘时平滑淡入、淡出。文字位置按完整周期预先规划，末帧与首帧也连续衔接；固定的标题、计数和颜色条文字保持不变。渐变只控制显示，不改变科学数据。

4K 静态图只显示当前主要标签位置，保持文字清晰，便于直接放入 PPT；动画及从 MP4 提取的审核帧保留真实渐变状态。刻度文字淡出时，坐标轴上的刻度线仍保留。

局部显示 68 个三维速度箭头，顶部为 `Local detail | FEM velocity vectors`，底部仅为 `68 vectors`。每个箭头由实际四面体重心处的 P1 速度矢量确定；箭头方向表达流向，颜色表达速度大小，箭头长度统一用于清晰显示。输入速度为 m/s，颜色条为 mm/s。

两幅图均保持原始物理坐标中的 Z 方向竖直，并绕通过各自视图中心、平行于 Z 的固定轴旋转。相机仰角为 6°，旋转期间中心、距离、仰角和缩放比例不变。视频为 1920 × 1080、24 fps、18 秒，转满一周；18 秒为展示时长。

刻度、三维参考网格及坐标标题使用原始物理坐标，单位为 µm。参考轴经过相机记录中的参考点，不将该点错误标为原始坐标零点；刻度为绝对坐标值。全血管刻度间隔为 20 µm，局部为 5 µm。颜色条仅包含 `Speed (mm/s)` 与数值刻度。

局部窗口以原速度局部图的分叉焦点为中心，使用 X/Y/Z 半轴分别为 18/18/11.5 µm 的椭球裁剪：横向取景宽度为 36 µm，聚焦分叉并覆盖两侧分支。

## 压力与 WSS

二者采用与最新全血管流线相同的固定相机取景比例，相比旧压力/WSS 动画，单位物理长度在画面中的显示比例增大约 33.4%。Z 轴保持竖直，入口出口标签旁置并用箭头指示；坐标轴、µm 单位、刻度、网格、透明文字及 0.875 秒渐变全部沿用。

压力顶部保持 `Full vessel | Pressure field`，底部仅为 `Steady FEM pressure on the vessel surface`。颜色条仅含 `Pressure (Pa)` 和数字；范围保持 −5～2500 Pa，数字刻度保持 0、500、1000、1500、2000、2500。

WSS 顶部保持 `Full vessel | Wall shear stress`，底部仅为 `Derived from FEM velocity gradient`。颜色条仅含 `WSS (Pa)` 和数字；范围保持 0～55 Pa，数字刻度保持 0、10、20、30、40、50、55。

使用已有 `Pressure_Pa` 和 `WSS_display_Pa`，不重算求解。压力绘制于包含端口的完整表面；WSS 仅绘制于物理血管壁，沿用由 FEM 速度梯度得到的壁面切向黏性应力模长，再按相邻三角形面积加权得到节点显示值。没有调整数值范围、裁掉异常值或平滑几何；曲面法向光照只影响显示。每段视频仍为 1080p、24 fps、18 秒，另外提供 4K 静态图。

## 参数位置

脚本顶部的 `SIZE`、`FPS`、`FRAMES`、`VIEWPORT`、`DETAIL_RADII_UM`、`ELEVATION_DEG` 控制画面、旋转与局部区域。`velocity_glyphs()` 控制速度箭头长度和空间采样；`CoordinateGrid` 控制网格；`annotate()` 控制固定标题与颜色条；`AnnotationPlan` 控制透明文字、端口指示箭头、避让和周期渐变。`FADE_RADIUS` 为渐变窗口的半宽，当前为 10 帧。

`fitted_z_center()` 和 `Turntable` 根据完整旋转周期自动选择固定取景范围，并检查全部几何在视野内。直接强行缩小相机尺度可能裁掉完整血管，请保留视野检查。

## 验证

```bash
./run.sh --validate
```

验证结果写入 `results/INDEPENDENT_VALIDATION.json`。图像生成流程也会逐帧检查相机、视野、视频帧数和不重复帧，并检查源文件未被修改。两个 `*_annotations.json` 记录每帧文字中心、外框、透明度和端口箭头目标；独立验证覆盖末帧至首帧的透明度与位置变化、端口文字遮挡和箭头投影准确性。`BUNDLE_VALIDATION.json` 记录当前独立版本的运行验证。
