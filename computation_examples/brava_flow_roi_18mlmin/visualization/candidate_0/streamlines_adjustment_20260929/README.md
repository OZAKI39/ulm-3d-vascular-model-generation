# candidate 0 流线动画显示调整（2026-09-29）

按用户最新要求，采用整周不超出画面的最大固定倍率；不再强制放大 30%。
保持原有 1920×1080、24 fps、固定 Z 轴、固定旋转中心、透明管壁、线宽、配色、0–600 mm/s 色标、英文文字及旁置渐隐箭头标签。
整周旋转时间由 18 秒改为 36 秒（20°/s → 10°/s）。

相机取景同时考虑真实壁面和现有流线的全部顶点，以整个 360° 连续旋转的解析投影包络计算固定倍率；在 1080p 画面边缘留 8 px，4K 静态图相应为 16 px。
不随角度自动缩放，不移动物理几何或端口，不裁切血管、不重新计算流线。
最终倍率和逐帧投影边界见 `CAMERA_FIT.json`；视频解码、数据哈希和标签检查见 `RENDER_VALIDATION.json`。

## 复现

在项目根目录执行（Python 需含现有 PyVista、VTK、NumPy、SciPy、Pillow、imageio-ffmpeg 环境）：

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow \
  /root/particle8_2_runs/env/bin/python -B visualization/render_streamlines_fitted.py --duration-s 36 --margin-px 8
/root/particle8_2_runs/env/bin/python -B visualization/render_streamlines_fitted.py --publish
```

上述命令只处理 candidate 0 的流线动画及对应 PNG/PDF，先写入本目录的 `rendered/`，校验后再替换页面所引用的文件。
其他流场动画、candidate 15、原始流场及微泡文件保持不变。未启动 CFD 或微泡计算。
`previous/` 保留调整前的相关文件，`PROTECTED_BEFORE.json` 记录受保护文件的原始哈希。

`MEDIA_VALIDATION.json` 的各个 `media` 条目分别记录帧数与时长：流线为 864 帧，其他动画仍为 432 帧。
总览页面仅对流线视频增加版本参数以刷新缓存，页面样式不变。

## 本次输出

固定倍率为原来的 1.016144349 倍（约放大 1.614%）。原始视图已接近完整一周旋转的横向取景上限，因此保持完整血管入镜时，可增加的倍率有限。

已检查九个方位的预览，标签避让检查通过；完整视频解码为 864 帧、24 fps、36 秒。
