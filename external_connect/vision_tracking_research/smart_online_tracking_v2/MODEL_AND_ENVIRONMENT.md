# 模型、来源与环境

> **最新续研：2026-10-16，已完成离线人工ROI评价。** 以下较早记录按历史保留；当前状态以文末续研和[完整报告](results/supervised_20261016_v1/REPORT.md)为准。

- Vision：Python3.12独立venv安装层，只读引用上级研究venv；OpenCV-headless4.10.0.84、Norfair2.3.0、Cellpose3.1.1.3、Torch2.11.0+cu128。有效锁定见evidence/requirements-vision-effective.lock。不修改Windows Python3.14环境。
- Annotation：独立Python3.12 venv，micro-sam1.8.14、napari0.9.2、PyQt6，NumPy2.2.6。独立锁定见requirements-annotation.lock，避免与Norfair的NumPy1.x要求冲突。两环境pip check通过。
- 原SAM ViT-B仅辅助标注。官方权重URL：https://dl.fbaipublicfiles.com/segment_anything/sam_vit_b_01ec64.pth ，SHA-256 `ec2df62732614e57411cdcf32a23ffdf28910380d03139ee0f4fcbe91eb8c912`。原SAM官方代码与权重Apache2.0；SA-1B数据有单独研究条款，本轮未下载训练数据。实际采用micro-sam加载，不安装/运行SAM2。
- Cellpose唯一模型cyto3，官方下载 https://www.cellpose.org/models/cyto3 ，SHA-256 `2dc3087a8abd7da46d1ab0ddd5824639933cc3ff63b382af3fa1939a392db93c`。3.1.1.3代码LICENSE为BSD3条件；该版官方README明确cyto3等模型训练涉及CC-BY-NC数据。不能从代码BSD许可推断所有训练数据或商业用途无额外条件。未下载Cellpose训练数据，未将用户视频上传第三方。许可原文快照和URL/SHA保留在evidence/licenses。
- micro-sam MIT，napari BSD3，Norfair BSD3，OpenCV Apache2.0；具体依赖还有各自许可。官方项目：https://github.com/computational-cell-analytics/micro-sam ，https://github.com/napari/napari ，https://github.com/tryolabs/norfair ，https://github.com/opencv/opencv 。micro-sam源克隆只用于核查，实际运行来自锁定的PyPI包，不能混淆两者版本。
- cyto3选择依据：[Cellpose3官方模型文档](https://cellpose.readthedocs.io/en/v3.1.1.1/models.html)。仅运行一个轻量替代模型，没有竞赛式扫描。

权重、环境和外部仓库未加入Git。中文界面使用用户Windows系统已有微软雅黑字体，通过用户字体目录符号链接供WSLg显示，没有改动设备工程。


## 2026-10-16 已有标注的离线续研（supervised_20261016_v1）

本轮复用既有cyto3 3.1.1.3及环境，未安装新模型家族。一次100 epoch SGD微调耗时15.613 s，验证失败，仍选预训练。拒绝检查点SHA `b1f80ccb122778b131ce0c99ed50fb77757d0fc0a08afc190402743c469a3c3c`；权重只在本地results内，未入Git。后台原SAM提示掩码与自动检测分开评价。

完整参数、失败证据、28项回答和输出路径见[本轮报告](results/supervised_20261016_v1/REPORT.md)。
