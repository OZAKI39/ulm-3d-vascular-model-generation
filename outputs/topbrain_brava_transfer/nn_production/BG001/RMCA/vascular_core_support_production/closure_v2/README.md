# s1-11 本轮交付入口

**状态：NEEDS_ADJUSTMENT。生产候选为空。这里的 STL / 3MF 是诊断文件，不是已验收的打印方案。**

首先阅读 [中文结果报告](vascular_core_support_safe_removal_report.md)。本轮完成六条既有姿态与间隙路线、32 次实际修补切片；没有通过局部接受条件的修补，原基线保留。诊断展示为 O000_R20_R3_C3_LAYER_I05：中心线 98.74371859%，端口 99.17355372%，逐层 99.36224476%，支撑自身仍有 284 个未承托区域。没有启用覆盖闭合后的拆除几何优化。

- [selected/](selected/)：独立血管、独立支撑、拼接审阅 STL，以及两个独立 normal part 的 X1C_project.3mf。Auto Support OFF。
- [QC_SAFE/](QC_SAFE/)：17 张检查图，包含实际失败点、真实组件、几何载荷路径及 Bambu 原生预览。
- [拆除限制指南](SAFE_SUPPORT_REMOVAL_GUIDE.md)：按真实组件编号列出阻碍；没有可直接照做的生产剪切顺序。
- [reports/](reports/)：任务原文、原基线、真实运行日志、重放日志及测试日志。
- [validation/](validation/)：199 项通过的测试记录、源文件前后哈希、交付校验值、两线宽和三线宽连接片的真实切片试验。
- [support_closure_iterations.csv](support_closure_iterations.csv)：每次实际指标、回退及停止证据。
- [candidates/](candidates/)：逐轮支撑、Bambu 3MF/G-code、逐点与逐层失败明细；失败记录保留。

命令：`./run_python.sh s1-11_vascular_core_support.py --stage closure`。仅重放已有结果：`--stage review`。退出码 2 表示制造检查未通过。详细接口见项目 `docs/S1_11覆盖闭合与安全拆除.md`。

源保护：485 个文件，0 个改变；未向打印机上传或启动任务。软件检查通过不等于物理拆除已验证。
