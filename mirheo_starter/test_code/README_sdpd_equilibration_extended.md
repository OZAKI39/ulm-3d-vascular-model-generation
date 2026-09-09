# 固定参数 SDPD 后段观察与恢复核查

此扩展复用现有运行入口、连续推进、完整块统计、共享 GPU 账本、授权登记和离线 Plotly 浏览器检查。原生 Mirheo/CUDA/MPI、单位、第一阶段血管与历史原始数据未修改。

当前构建的 SDPD 相互作用没有 checkpoint/restart override，继承 MirObject 空实现，无法证明其内部 RNG 被保存/恢复。内核存在序列化函数不代表相互作用实际调用它们。源码与本地 `.so` 符号证据写入 `restart_contract.json`。因此正式恢复链被阻止；文件归档完整或 CPU 测试通过均不能替代真实恢复验证。

旧运行没有原生 checkpoint。新设计为独立冷启动，到 `t*=0.80` 停止；固定正式窗口 `(0.60,0.80]`，趋势窗口 `(0.40,0.60]`，200 步采样。保留原始门槛，生产 `selection=null`。

## 已实现

- 历史时序、固定窗口描述量、HAC 趋势敏感性、逐观测量状态、压力成本条件估算。
- 六个同相位快照的保守/耗散 virial CPU 重建及随机残差；与原 CPU 审核交叉核对。
- 退出后 checkpoint 完整版本归档、内部引用和哈希、粒子 ID 对齐、时间验证、绝对样本链、正式链停止保护。
- 原生非零 dt 三进程对照：A 4000 步；B 保存 2000 步状态后退出；新进程恢复再推进 2000 步。首步和短段的容差在 YAML 中事前固定，真实 GPU 测试尚未运行。
- 逐段调度与共享累计预算。当前原生 RNG 契约使正式调度在第一道检查处停止；并未伪造可恢复能力。
- 每次 `u.run` 重建 checkpoint 调度器，200 步块起点都会保存。原生 scratch 用 PingPong，退出后归档成完整独立版本；保存触发额外一步明确计费、排除在承诺轨迹之外。该 I/O 成本尚未实测。

## 命令

```bash
cd /home/lzy/projects/mirheo_starter
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --preflight-only
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --analyze-only
.venv/bin/python -B -m test_code.review_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --browser-test
.venv/bin/python -B -m test_code.review_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --open
.venv/bin/python -B -m unittest discover -s test_code -p 'test_*.py'
```

最终包通过明确的配置、delivery_record、manifest 和内容哈希定位，不按 mtime 选结果。新增报告属于历史补充证据；不替换旧 `(0.30,0.40]` 结果。

## 授权与限制

新增授权当前是 0。完整条件请求上限 1450 秒：恢复对照 3×30 秒，条件长实验 4×340 秒；原余额 71.715047528953 秒，需要追加整数 **1379 秒**。旧 493 秒仅限旧实验且已使用完毕。并发 1、单任务 ≤600 秒。

此请求不会自动消耗全部额度：只有原生 RNG 持久化已经证明且真实恢复 PASS 时，才允许四个正式段。当前构建缺失该能力，所以**批准预算本身不能解除长链阻塞**。失败仍计费；不自动重试、加长观察或改成超过 600 秒的单进程。

真实用户批准后，将真实消息保存为文本，并制作明确绑定 `budget_request.json` SHA-256 与其中 `authorization_scope` 的确认 JSON，字段沿用现有授权登记接口：`authority=explicit_user_message`、`user_message_file`、`additional_seconds`、`budget_request_sha256`、`scope`。配置和脚本不会自行创建这份确认。

```bash
# 仅在真实批准并形成证据文件后
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --register-authorization /absolute/path/to/real_confirmation.json
.venv/bin/python -B -m py_scripts.run_sdpd_equilibration --config py_scripts/sdpd_equilibration_extended.yaml --execute
```

恢复诊断用 `diagnostic_restore` 明确允许试读缺失 RNG 的粒子/时间归档，以测量不连续性；该例外只适用于诊断任务，验证结果仍由 RNG 契约强制判为 `RESTART_NOT_VALIDATED`，不能进入正式链。

旧 CPU 验证记录只证明旧代码。当前只读重分析遇到旧证书时不再报错或谎称获得当前验证，而记录无当前 CPU 证书；执行入口仍必须有与全部当前代码匹配的 PASS 记录。该兼容修复有修改前的失败回归记录。

浏览器自动核查与人工验收分开，人工始终 `PENDING`。本轮不提交或推送 Git。
