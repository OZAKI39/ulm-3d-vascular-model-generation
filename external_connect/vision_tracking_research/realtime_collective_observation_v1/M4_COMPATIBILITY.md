# M4 观测状态

F10Y_M4_1.0 是新观测 schema，不冒充旧有带身份泡团的 v1/v2。保留旧数学与 TTL 思想，但目前没有可信材料 ID，因此 material_id=null、six_parameters=null。输出群体 mean velocity 仅在有效 PIV 与足够候选支持时有效。

字段分别包含 valid（运动）、detection_valid、aggregation_valid、state_age_ms、ttl_ms、session_id、frame_id、published_monotonic、clock_scope、reason、velocity、aggregation、training_ready=false。publisher 绑定随机 localhost 端口，记录在 session/publisher.json，使用原成熟 pyzmq PUB/SUB。

消费者必须自行按接收后的本地 elapsed 加消息 state_age_ms 继续检查 TTL。断流/EOF/暂停/停止/会话变化/处理错误/过期时数值不能续用。跨主机不能拿 published_monotonic 直接相减；本版没有部署跨公网 M4。

测试用实际 SUB socket 收到了无效化消息，并检查速度变为 null。数据源仍无法提供相机曝光时刻；这不是世界模型训练或闭环控制准入。
