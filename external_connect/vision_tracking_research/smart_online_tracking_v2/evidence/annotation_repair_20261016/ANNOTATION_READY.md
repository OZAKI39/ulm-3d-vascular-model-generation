# 人工标注工作台已准备

使用官方 **micro-sam 1.8.14 Annotator2D + napari 0.9.2**，加载原版 SAM ViT-B（不是 SAM 2）。本工程只补充类别分配、版本化保存和换帧按钮。窗口已经在 WSLg 实际打开；保存和重新读取一致，详见 `evidence/annotation_gui_smoke.json` 和 `evidence/annotation_window.png`。

**当前正式专家标签为 0 张：EXPERT_GROUND_TRUTH_PENDING。** 软件窗口检查保存于独立测试目录，不计入人工标签。模型生成的任何候选都不是真值。

## 启动

在 WSL 终端执行：

```bash
bash /home/lzy/projects/external_connect/vision_tracking_research/smart_online_tracking_v2/start_annotation.sh
```

从第 11 张继续：在命令末尾添加 `--index 10`。只用 napari 手工描画：添加 `--manual-only`。不需要连接相机，当前候选帧已经保存在本地。首次打开一帧需要计算图像嵌入，后续可从缓存读取。

## 逐张操作

1. 窗口标题显示当前序号、图像 ID；左下方显示训练/验证/测试组。先浏览整幅图，滚轮缩放，按住空格拖动画面。
2. **确认一个真实单泡**：右侧 SAM 工具已装载原版 `Natural Images (SAM)` 自定义权重。选中 `point_prompts` 层，在泡内点击一个正点，按 `S` 分割；必要时用负点修正。查看 `current_object` 是否准确贴合可见边界。
3. 左下 `Object class` 选择 `SINGLE_BUBBLE`，点 **Assign SAM current mask to class**。该区域成为待审核草稿。在 `point_prompts` 层清除上一个目标的点，再标下一个对象。
4. **确认微泡团**：能从实际观察确认多个泡构成一个可见团时，圈出其二维外轮廓，分配到 `AGGREGATE`。仅因像素粘连或两泡投影重叠时不能确认，使用 `UNCERTAIN`。
5. 背景纹理/杂质选 `BACKGROUND_ARTIFACT`；无法判断选 `UNCERTAIN`；明显失焦且边界不可靠选 `OUT_OF_FOCUS`。可直接选相应 napari Labels 层，用多边形（快捷键 `3`）或画笔（`2`）修改，用 `M` 换到新实例编号，橡皮擦删除错误区域。
6. 随时点 **Save draft** 或 `Ctrl+S`。该操作不表示完成专家审核。每次保存保留独立修订文件，不覆盖前一版。
7. 检查整个原始画面是否漏标。填写真实审核人名字，勾选 **I reviewed the full frame, including missed objects**，再点 **Confirm full-frame expert review**。类别像素冲突会拒绝确认。空画面也只有在确实完整核查后才能确认。
8. 点 **Next frame** / `Alt+→` 继续，或 **Previous frame** / `Alt+←` 返回。换帧会保存有修改的草稿；重开已确认而未改动的帧不会自动取消确认。关闭窗口前请先保存。

如果整幅图过于密集，先保存草稿和备注；不要勾选“整幅检查完成”来代替局部检查。本轮正式全图指标只使用完整人工核查的图像。先完成少量高质量帧比匆忙确认所有帧更有用。

## 短片段怎么标

候选中有连续 8 帧片段，分别为原视频 120–127、450–457、780–787 帧。只跟踪能确实识别的少量目标即可。

对同一片段中能可靠对应的对象，在左下填写 `Optional short track ID`（如 `A`），选择该类别和该帧 `Instance ID for short track`，点击 **Set short track ID**，然后保存。不同帧的实例编号可以不同，短轨迹 ID 必须人工确认；看不清就留空，不凭运动方向猜身份。

## 数据保存与继续开发

标注位于 `data/annotations/<image_id>/`：每次保存生成 NPZ 掩码、JSON 修订记录和 `latest.json`。五类掩码保存原始像素坐标，没有缩放。重新打开同一帧会自动读取最新修订。不要把 `evidence/gui_smoke_sandbox/` 当成正式标注。

初始提供 72 张视频候选，训练/验证/测试各 24 张。后续实时候选单独隔离，不自动混入训练集。现有素材尚未被专家确认覆盖各种焦平面、突然转向和聚集/分离条件。

完成后在聊天中告诉 Codex“已完成哪些图像/多少张标注，可以继续核查与训练”。Codex 将检查真实保存文件、确认状态、类别冲突、数据划分和标签 SHA，再冻结评估；不会仅凭一句“开始”就假定已有专家真值。未经审核的预测不得用来填写 Precision/Recall 或团分割准确率。
