# 复现与恢复

本轮服务器根 `/workspace/y_clot/project`，执行目录 `analysis/physics_constrained_load_inverse_design_v1`，输出 `verification/physics_constrained_load_inverse_design_v1/20261016T144922Z`。本地审查副本保存全部本轮轻量结果、SQLite、原生状态及图表；既有原生库留在服务器，不重复迁移无关大文件。

## 读取与重建图表

```bash
source /workspace/y_clot/env.sh
cd /workspace/y_clot/project/analysis/physics_constrained_load_inverse_design_v1
python analyze.py
```

不会重新运行力学；会从已保存原生状态和历史生成图表/报告。

## 恢复同一研究

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/workspace/y_clot/envs/inverse_opt_v1/bin/python /workspace/y_clot/project/analysis/physics_constrained_load_inverse_design_v1/run_inverse_design.py --limit 4 --next
```

已完成四条有效样本时不再调用新原生求解。SQLite与SAMPLER_STATE.pkl需要一起保留；pickle仅加载本项目自己生成的文件。完整原生作业由supervisor管理，配置/脚本保存于本轮目录。单个原生作业若中断，研究恢复会用原参数从参考态重算；不完整历史先移入interrupted_attempt目录，不声称半途机械checkpoint恢复。已有SUMMARY直接读回，无效数值作业为FAIL而非零目标观测。

首次运行顺序为 `run_inverse_design.py --prepare`（既有clot环境）→ optimizer `--limit 2` →退出/恢复→ `--limit 4 --next`。若需要全新独立科研复现实验，先复制此隔离工作流并修改common.py的RUN到新目录，再冻结设计；不要删除或覆盖本轮DB。

## 实际MEVO资产检查

```bash
source /workspace/mevo_acoustic_microbubbles/activate.sh
python /workspace/y_clot/project/analysis/physics_constrained_load_inverse_design_v1/mevo_feasibility.py
```

只读当前实际APECSS/群集数据及历史k-Wave基底。源路径、原始SHA、执行commit见MEVO_ASSET_AUDIT.json，保存的紧凑实测数据见MEVO_ACTUAL_EVIDENCE.npz。

## 环境和范围

优化器独立venv `/workspace/y_clot/envs/inverse_opt_v1`，只读复用 `/venv/main/lib/python3.12/site-packages` 中已安装torch，新增Optuna/BoTorch等安装在独立venv；没有改写/升级原torch环境。`optimizer_requirements.txt`包含继承可见包，ENVIRONMENT.json列本轮实际使用版本。原生子进程仍使用 `/workspace/y_clot/envs/clot/bin/python`，通过LD_LIBRARY_PATH首先加载冻结V3库，运行时逐路径核对SHA。MEVO使用自己的venv。

`DESIGN_FREEZE.json`先于任何力学试验；`PORT_ADAPTER_FIX.json`保留初始化名称修正及当时源码SHA。trial0运行0步被拒绝，trial1是完全相同优化器候选的有效重试。最后增加的中断归档分支只处理无SUMMARY的不完整历史，不改变已执行计算；已执行mechanical_evaluator版本另存provenance。

所有图PNG/SVG及OPEN_RESULTS.html采用实际保存数据。不存在新完整V3或新的MEVO声场扫描。本轮完成后停止。
