# BraVa CFD 独立环境

已完成迁移与 CPU/MPI/CUDA 小型验证。当前项目位于 `/workspace/brava_cfd/project`，Python 环境位于 `/workspace/brava_cfd/envs/cfd`。

```bash
source /workspace/brava_cfd/activate.sh
cd /workspace/brava_cfd/project
python -c "import sys; print(sys.executable)"
brava-mpi -n 2 /bin/hostname
```

当前求解器入口是 `/workspace/brava_cfd/bin/brava-solver`，MPI 入口是 `/workspace/brava_cfd/bin/brava-mpi`。旧 `project/gpu_solver_fix/svmultiphysics` 保留为 CUDA 13 历史证据，未被新构建覆盖。后续计算使用新的输出目录；不在冻结结果目录中直接运行。

[迁移与验证报告](MIGRATION_AND_VALIDATION_ZH.md) · [当前机器状态](MIGRATION_STATUS.json) · [独立环境版本](reports/requirements.lock.txt) · [完整性检查](reports/PAYLOAD_INTEGRITY.json) · [环境验证](reports/ENVIRONMENT_VALIDATION.json)。

开发前先读 `project/AGENTS.md` 和 `project/ACTIVE_DEVELOPMENT_SCOPE.md`。原科学门槛仍有效，迁移通过不代表局部场精度通过，不自动恢复旧队列。当前原生 CFD/Python 字段处理环境已验证；旧 SimVascular GUI/嵌入式 Python、微泡批次环境及其它工程未迁入。
