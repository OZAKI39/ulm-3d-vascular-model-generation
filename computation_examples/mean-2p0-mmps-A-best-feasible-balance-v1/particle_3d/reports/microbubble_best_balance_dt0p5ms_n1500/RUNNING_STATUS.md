# 本批计算与交付已完成

1500 条，名义 dt=0.5 ms；O1=0、O2=133、O3=1182，接触支持静止 185 条，执行失败 0。保持原有粒径分布，没有补充小粒径覆盖样本。

服务器交付标记：`DELIVERY_COMPLETE.json`；本地独立核验：`LOCAL_VERIFICATION.json`（PASS）。已核验 9000 个轨迹结果文件、114 个原有源码文件及 288 帧视频解码。

结果入口：`OPEN_RESULTS.html`；中文说明：`MICROBUBBLE_RESULTS_ZH.md`；复现与路径：`REPRODUCE.md`。

8 CPU 工作进程完成运动积分；RTX 4090 完成 FP64 数据复核及 NVIDIA EGL 渲染。视频使用 CPU libx264 编码。未启动 CFD。
