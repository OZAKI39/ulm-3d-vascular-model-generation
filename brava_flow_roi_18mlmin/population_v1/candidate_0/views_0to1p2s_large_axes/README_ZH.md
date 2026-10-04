# 四视角动画：放大坐标轴文字

本次仅调整坐标轴文字的显示尺寸，并重新渲染原有 4 个动画：

| 内容 | 原尺寸 | 新尺寸 |
|---|---:|---:|
| 坐标轴名称及单位 `X/Y/Z (mm)` | 24 px | **40 px** |
| 坐标刻度数字 | 18 px | **30 px** |

均放大约 **66.7%**。字体、文字颜色、网格颜色及线宽保持原样。刻度数字到其真实刻度位置的屏幕间距从 23 px 增至 38 px，以容纳放大的文字；坐标轴标题按新文字尺寸重新避让，拥挤或被血管遮挡的刻度文字仍采用原有避让规则。刻度的物理位置、数值及单位不变。

4 个方位角仍为相对基准的 240°、330°、60°、150°，仰角均为 25°。所有相机位置、主图大小、右上角角度标注、粒子颜色和大小、透明度、物理时间 0–1.2 s、逐帧低速筛选和慢放速度均与上一版一致。

每个动画为 1920×1080、24 fps、262 帧，时长 10.916667 s。未启动 CFD 或轨迹求解。

## 检查和附件

- `OPEN_RESULTS.html`：四视角播放和下载。
- `CONTACT_SHEET.png`：四视角 1.2 s 末帧总览。
- 各角度目录的 `animations/`：新版 MP4；`figures/`：原始渲染及编码后抽帧。
- 各角度目录的 `data/layout.json`：字号、刻度间距和相机；`coordinate_annotations.json`：真实文字边界框。
- `data/QC.json`：4 个视频完整解码通过；坐标文字无相互重叠、无裁切；与旧版的逐帧位置、速度、可见性和相机逐项一致。
- `data/provenance/`：修改前的绘图源码副本。

原始轨迹、流场和已有动画保留。流场共享绘图代码只新增了可选字号参数，默认值仍为原来的 24/18 px，不改变其他调用的默认样式。

## 复现

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 env/bin/python -u -B scripts/render_population_four_views.py --output-dir /workspace/brava_flow_roi_18mlmin_20260928/population_v1/candidate_0/views_0to1p2s_large_axes --axis-title-size 40 --axis-tick-size 30 --axis-tick-offset 38
OPENBLAS_NUM_THREADS=2 env/bin/python -u -B scripts/check_population_four_views.py --output-dir /workspace/brava_flow_roi_18mlmin_20260928/population_v1/candidate_0/views_0to1p2s_large_axes --compare-with /workspace/brava_flow_roi_18mlmin_20260928/population_v1/candidate_0/views_0to1p2s
```

本地：`/home/lzy/projects/brava_flow_roi_18mlmin/population_v1/candidate_0/views_0to1p2s_large_axes/`。
服务器：`/workspace/brava_flow_roi_18mlmin_20260928/population_v1/candidate_0/views_0to1p2s_large_axes/`。
