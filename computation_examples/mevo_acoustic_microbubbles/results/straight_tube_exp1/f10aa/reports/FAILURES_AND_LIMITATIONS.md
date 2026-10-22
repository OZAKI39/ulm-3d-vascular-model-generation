# 失败、负结果与适用范围

1. 新隔离环境首次 profiling 在依赖安装未完成时启动，两个 attempt 仅 import 失败，无物理输出。随后使用已验证官方已装包独立复制，恢复记录保存。旧项目环境未改动。
2. 运行时切换到四线程的 OpenBLAS Cholesky 发生 SIGSEGV；随后在进程启动前固定线程数，独立测试N32/N64/N128及原B0通过。失败不是接受物理步，保留 faulthandler 和输入 SHA，不降容差绕过。具体结果：[(32, 'B1', 'PASS', 0), (64, 'B1', 'PASS', 0), (128, 'B1', 'PASS', 0), (128, 'B0', 'PASS', 0)]。
3. 初版 CuPy errstate 名称错误、Warp profiling callback indent 参数错误，分别修复；保存旧日志与不可变提交，已完成固定测试未选择性重跑取快值。
4. CuPy 在当前规模比优化CPU慢；GPU几何与因子转移未形成生产后端。Packed读回数量减少但字节略增，128泡没有可论证的额外整体收益。
5. 单条正式后端只有一次完整重复、每N一套主要分散初态；接触制造状态不是新物理证据。性能改善不证明 PDMS/自由液面/声流/气泡接触/电声映射已校准。
6. 256 仅合法短步预检，长轨迹 NOT_RUN。没有秒级、AI/world-model训练和真实PZT闭环。
7. N128延长原CPU参考在1.659ms耗尽240s预算，不能称完整2ms已通过参考；优化版dt/cache完整通过也不补齐缺失参考。
8. 图/视频显示的是2ms物理过程放慢播放，未通过插值、聚类简化、球径放大或篡改轨迹展示“聚成大团”。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
