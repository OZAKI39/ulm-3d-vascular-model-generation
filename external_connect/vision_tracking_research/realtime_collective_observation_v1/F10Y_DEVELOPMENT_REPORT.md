# F10Y-Vision 开发与验证报告

已交付可操作的Qt群体观测系统、两种统一输入、有限回放与故障验证、科研图和视频。真实相机与真实声学聚集均未验收；合成图像测试为19/20，保留额外拒判失败。

在家启动：

```bash
bash /home/lzy/projects/external_connect/vision_tracking_research/realtime_collective_observation_v1/start_replay.sh
```

科学执行快照：`40ce97aab4a8416f8319aecb89e674c84304ab12`。最终UI执行快照：`2b3707bfa0c7aa31e27aaacdb1d0f9eb708d26e3`。分支 `dev/f10y-realtime-observation-20261017`，交付提交由DELIVERY_IDENTITY.json记录。没有自动推送GitHub。

1. 是。Windows 的 LatestFrameStream 保持原代码，由原 windows_frame_server.py 调用；本轮只运行合成替身验证，没有访问相机。

2. 是。ConfiguredWindowsBridge 直接继承 WindowsLatestBridge._loop/read_exact/snapshot，扩展可配置路径、有限会话轮换与回收。

3. 没有修改原相机工程、配置或服务。两处原相机代码 SHA 前后相同；只在新 F10Y 模块开发。

4. 是。LIVE_CAMERA 与 REALTIME_REPLAY 共用 FramePacket、Pipeline、Observer、Qt 界面。失败不静默切换来源。

5. 是。PyAV读取每帧真实PTS，独立线程按PTS/speed调度；0.5/1/2×均实测，暂停恢复会平移调度基准。

6. 是。原片937帧，1476×1120；PTS0–31.1688s、间隔33.3ms，完整科研视频逐帧保留原PTS。

7. 是。源/分析之间单最新槽，记录队列32；正式1×接收937，处理693、跳过244，记录丢弃0，分别统计。

8. 是。无新帧、暂停、EOF、断流、会话切换、错误和TTL过期都会使运动/聚集输出失效，ZMQ下游收到null而非补零。

9. 每个源独立几何；录像mask绑定SHA和分辨率。LIVE初始无验证mask，需原图绘制/加载并设置固定壁面参考；不自动缩放旧mask。

10. 是，全部937帧完成独立完整分析；实时模式只对实际进入分析的帧做全通道扫描，不用小ROI冒充全通道。每帧有效PIV平均431.9窗，检测与有效速度覆盖分开。

11. 是，复用F10X冻结OpenCV候选检测器。质量分数是启发式，弱对比和密集重叠漏检限制保留，未审核区域没有被当新负例训练。

12. 可以在PIV支持足够的区域输出图像平面平均方向、速度及一致性；没有支持时null，不建立长期单泡ID。

13. 可描述二维候选分布收缩与局部汇聚证据；真实录像不足以证明物理聚集。

14. 是。按有效观测面积统计密度，并输出最近邻/kNN、Rg、协方差、长短比和方向。

15. 可以显示DBSCAN候选密集区域及空间持续性，但不是物理泡团身份或融合。密集后发生检测合并时可能失去可观测性。

16. 部分通过。冻结20例：数学20/20，实际图块+检测+CuPy PIV图像链路19/20；误报0，共5例拒判，其中1例为额外漏检。FINAL_GATE保守设SYNTHETIC_AGGREGATION_TEST_PASS=false。

17. 可能误判的风险包括恒定数量的候选替换、光学叠加、密集检测合并、视野进出、反射与非仿射PIV伪向量。本批密集团seed95001候选60→40，触发数量变化拒判，保留为失败。

18. 是，在合成纯平移、旋转对照中没有误报收缩；先移除平均平移，再检查径向与散度，并要求NN/Rg变化一致。

19. 允许描述可见数量变化，但变化超过门槛时聚集结论拒判；这是避免消失/漏检产生假聚集，不是强行保持固定泡数。

20. 是。保留SNR、纹理、空间验证、mask/注册、时间连续性和TTL；没有把无效向量填零。

21. 是。WSLg实际可见窗口和Qt按钮/关闭测试已运行。start_replay.sh一条命令启动；最新原图持续刷新，淡色较早图层注明分析帧和年龄，Measured frame可查看严格对应帧。

22. 正式1×输入约30.030Hz（937帧/原始时间轴）；分析22.147Hz，新鲜有效运动22.019Hz。2×约60.06Hz仅为软件压力输入，不是实机帧率。

23. 逐阶段p50/p95详见REALTIME_PERFORMANCE.md。解码、全通道扫描、PIV、聚集、发布、异步序列化分别测量；Windows→WSL合成BGR管道另测，不将其充当相机延迟。

24. 正式headless decode+analysis p50/p95=48.16/56.05ms；release→完成软件年龄=61.98/78.63ms；新鲜运动22.02Hz。最终GPU实际UI+录屏有效14.25Hz、UI刷新19.83Hz；Qt画布paint p95=10.22ms。未达到完整33.3ms p95目标。

25. 正式1×headless跳过244帧；最终UI+录屏跳过459帧。异步录屏/GUI负载单独列出，不能混用。

26. 正式headless CPU平均125.8%（100%=一逻辑核），GPU设备平均4.40%。本机Rust+可见GUI实测有效4.24Hz，明显慢于服务器；本机RTX4060/CuPy缓存后实测有效10.35Hz，一键脚本已优先使用本机GPU；不能把服务器22Hz宣传为本机速度。

27. 是。最终核心+Qt循环605.22s，20个会话，处理11425帧，记录丢弃0，错误None，线程正常回收。RSS 768.1–812.0MiB，GPU内存762–762MiB。最后仅显示/paint统计改动另做完整33s录屏和Qt回归；长期科学核心与正式快照一致。

28. 是。事件按钮只记录JSON/曲线，代码未连接AFG或PZT输出。

29. 是。UI演示定时事件为SIMULATED_EVENT；实际Windows传输测试是合成图像。前后窗口比较始终SOFTWARE_SEGMENTATION_ONLY，不给声学因果结论。

30. 是。至少6组PNG+PDF、3段完整解码MP4、真实Qt界面截图；见figures/、videos/、evidence/video_validation.json。完整真实视频937帧，合成视频610帧，错误案例可见。

31. 使用start_live.sh，先按原流程确认相机，然后关闭占UDP口的旧预览，点击LIVE Start；核查现场变化、分辨率、mask及固定参考。详细15步见LAB_CAMERA_INTEGRATION_GUIDE.md。

32. 统一数据协议、实际PTS释放、有限队列、跳帧计数、过期/会话、源切换、坏包、慢处理、Qt操作/关闭、全通道分析、图像合成趋势、事件、记录和长期回放已验证。

33. 真实相机现场新鲜度、曝光时钟/端到端光学延迟、实际60Hz接收、LIVE遮罩/注册、长时间实机重连、真实二维聚集精度及真实PZT事件同步仍待实验室验收。

34. 没有。完整真实录像得到{'NOT_OBSERVABLE': 609, 'NO_CONSISTENT_CONVERGENCE': 328}，没有足够一致汇聚证据，也没有可信PZT日志；REAL_PZT_ON_OFF_RESPONSE_VERIFIED=NOT_VERIFIED。

35. 需要同光路尺度标定、现场接入验收、设备执行日志/同步误差、预登记OFF/ON对照、独立人工二维事件真值和密集失焦/合并失败对照。之后再考虑经带宽/压缩验证的远程GPU接入；不直接开展闭环或世界模型训练。

完整限制和Gate见FINAL_GATE.json；所有通过状态都有适用范围。
