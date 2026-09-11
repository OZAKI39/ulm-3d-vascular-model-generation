# 授权后的修复与复验

**尚未完成同等质量比较，qualified_speedup=null。** 授权已使用；本次停止原因是物理质量关口失败，不是等待授权或显存不足。

原生库已隔离编译，未替换原安装。新增 local→halo 顺序、空 halo 早退和原生出错状态记录；另外发现并修正项目分段观测：Mirheo 每次 run() 都重新分类并清除对象力，因此改为每阶段一次连续调用、由原生插件保存中间状态。旧代码与错误的 classifier_corrections=0 原始字段保留，解释见 [run_boundary_audit.json](run_boundary_audit.json)。一个 CPU 行为回归先失败、修复后通过。



实际控制、首次复发、半步结果与质量失败详见 [report_zh.md](report_zh.md)。原生补丁为 [local_before_halo_v2.patch](native/local_before_halo_v2.patch)，未改容量或反弹公式。新 Python 连续观测代码位于 py_scripts/single_rbc_repair/continuous_worker.py 与 continuous_protocol.py，各次实际源码、spec 和执行记录保存在对应 runs 目录。

[授权前修复记录原文](fix_log_before_authorized_controls.md) 保留，不将短测通过称为完整修复。
