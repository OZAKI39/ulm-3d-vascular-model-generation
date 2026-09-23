# Particle-8.2 计算来源审计

本地：localhost，可用逻辑 CPU=24，RAM=7991440 kB。职责仅开发、轻量测试、至多 10 条 smoke、调度和审核。
远程：f7c62a262077，SSH root@50.115.148.16:4159，可用逻辑 CPU=16，RAM=64963600 kB。正式大计算必须在该主机执行；容器内存/CPU 上限另见 JSON。

**PARTICLE8_1_HEAVY_COMPUTE_HOST = NOT_PROVEN**。P8.1 环境文件记录了 WSL 平台，但没有逐轨迹/分片的主机身份与 SSH 计算回执。这不足以证明正式重积分在远程执行，也不能仅据此断言全部重积分发生在 WSL。
已检查 2235 个历史元数据、日志和清单，保存其 SHA256；历史文件不改写。

完整 hostname、uname、nproc、lscpu、free、df、Python、git、ulimit 和 GPU 查询原始输出见 COMPUTE_PROVENANCE.json。
后续正式分片须携带主机、PID、worker、源码提交、配置与 Frozen 输入散列；任何本地正式分片都会使验证失败。
