# M4 状态兼容

新增M4-channel-3.0，包含窗口并行数组、mask哈希/有效性、safe/active/valid/fresh/stale/age、实际PTS、ROI、backend/device、每组局部six及reason。原v2/v1消费者不被直接替换。

兼容适配必须显式选择当前局部空间群；不选择则输出EXPLICIT_LOCAL_GROUP_REQUIRED和null，禁止自动把全通道拟合为一个群。在线演示明确记录选择的第一个有效空间群key，交由原ZMQ MockConsumer。session、frame_id、PTS、calibration_status和模型身份保留。

100ms TTL后参数null、valid=false、stale=true；帧不连续清空历史，不保留假新鲜速度。tests/test_contracts.py实际验证原JSON schema、显式局部适配和expiry；回放保存consumer_rows/watchdog。consumer接收成功只表明接口工程可用，不表示科学控制准入。
