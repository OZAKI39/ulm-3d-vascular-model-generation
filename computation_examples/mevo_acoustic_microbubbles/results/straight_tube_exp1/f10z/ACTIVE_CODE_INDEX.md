# F10Z active code

- `common.py`：路径、执行身份、历史源码和输入SHA保护。
- `prepare.py`：确定种子初态、动作和预算预登记；仅新空目录运行。
- `run.py`：继承F10S积分器，A对照、B收敛/动作、C有条件时长扩展。
- `diagnostics.py`：原场/响应/F2的CPU与GPU核查、冻结阻力分解和一致梯度检查。
- `state.py`：真实粒子ID、时间、SI分布指标、XY投影与HDF5导出，不补全失败轨迹。
- `remote.py`：Git执行归档、输入哈希、独立Supervisor标签；不会停止其他任务。
- `archive.py`：回传及逐文件SHA核验。
- `analyze.py`：运动、网络事件、形状、数值噪声、局部灵敏度和分组预测基准。

测试：`tests/f10z/test_interfaces.py`。物理核心仍为F10S_ACTIVE_CODE_MAP所列文件，不新增第二套粒子引擎。

- `enrich.py`：复制原始列定义及压力证据；采用原平移黏度重建阻力并与原残差核对。
- `quality.py`：实际力分解、形态/主轴、收敛与资源汇总。
- `figures.py`：8 组 PNG/PDF，只有真实保存数据。
- `videos.py`：原始样本帧映射、慢放、H264 全解码和 PTS 检查。
- `dependencies.py`：实际环境分发版本、许可证与包 RECORD 身份。
- `deliver.py`：35 问中文报告和分层 Gate，不制造条件性模型结果。
- `finalize.py`：旧源码/数据保护、逐状态来源、PDF 检查、安全进度更新和 SHA 清单。
