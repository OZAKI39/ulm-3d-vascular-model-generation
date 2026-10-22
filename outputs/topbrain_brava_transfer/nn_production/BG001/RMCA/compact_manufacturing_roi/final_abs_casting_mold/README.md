# BG001 RMCA 一体式 ABS 倒模件

当前状态：READY_FOR_HUMAN_REVIEW。四面墙与底板封闭，顶部开放，血管及四个端口与盒体为同一个实体。

请先阅读 [中文评估报告](abs_casting_mold_report.md)，再查看 QC/04、05、06、09、12、14、15、16 的连接、间距、支撑通路与切片图。

在 Bambu Studio 中打开 final/BG001_RMCA_BALANCED_ABS_casting_mold_with_vascular_core.3mf 可检查本机实际生成的切片。它使用 P1S 0.4 mm、Bambu ABS、0.20 mm Standard 工艺，并已选择侧放姿态。未发送给打印机。最终打印模型为同目录的 *_oriented.stl；不带 oriented 后缀的 STL 是开口向上的倒模坐标版本。reference/ 中是单独参考零件，diagnostic/future_PDMS_volume.vtp 是未来 PDMS 体积预览，均不要与最终整体重复叠加打印。

打印后需将整体转回开口向上再灌注。盒内支撑约 16.24 g，顶部可达性仅经过几何代理检查，尚未证明实物支撑可拆除。O3 附近最薄 PDMS 约 1.688 mm。

复现时，在项目根目录使用原虚拟环境并指定一个新的输出目录（已有目录会拒绝覆盖）：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python s1-6_abs_casting_mold.py --output-root outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold_review_02
```

本次完成 134 项测试；XML 与完整日志位于 test_results.xml 和 test_results.log。真实 Bambu 命令、原始 3MF、G-code 和支持量在 bambu/ 与 tables/ 中。最初未应用喷头坐标偏移的误报已修复，初次检查证据保留在 diagnostic/initial_coordinate_audit/，复核未改动任何原始切片。1701 个受保护源文件哈希全部不变。
