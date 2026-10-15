# P0代码安全审计与上下文精简

审计实际候选 28 个文件；**安全删除0个**。已有正式代码、离线数学模型、复现脚本、冻结source_snapshot及失败implementation_history保留。相同SHA的开发副本与冻结副本具有不同可追溯角色，不属于无价值重复。CMake编译探针不能仅凭临时名称证明已经无用，保留UNCERTAIN。

完成了逐文件SHA、Git状态、受限Python import/C++ include和运行入口检查；详细路径/依赖见CODE_CLEANUP_INVENTORY.json。删除白名单为空，没有删除操作或需要恢复的文件。历史科学数据仅流式SHA，不将大型日志/NPZ加载到上下文。

创建根目录ACTIVE_DEVELOPMENT_ENTRYPOINTS.md，共27行；只需优先读取当前模块及其直接依赖。旧PROJECT_ACTIVE_CONTEXT的历史叙述未改。保留全部现有未提交修改；Git tracked diff前后完全一致。当前源码、构建配置、冻结model导入与参考态F/U零值检查通过。旧生产与两套已接受诊断库SHA均核对一致。

本轮P0严格先于任何新C++动力学功能。零删除是审计结论，不将保留历史证据称为清理失败。
