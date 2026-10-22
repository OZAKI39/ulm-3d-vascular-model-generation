# Phase 8：独立 realization GPU batch

已实际运行 `test_phase8_batch_gpu.py`，batch=2/4 的来源 ledger、每组时间分段一致，末状态仅浮点舍入级差异。GPU kernel 叠批 field/wall/resistance，而非把所有组塞进同一时间步或相互碰撞。每组 RNG/address/clock/particle IDs 独立。

性能另在 phase9 中测试 16×1 泡、各 20 步；不是生产额外注射。详见 `../server/phase8_batch_gpu.json` 和 `../server/phase9_sparse_batch_benchmark.json`。无真实多组统计置信区间研究。
