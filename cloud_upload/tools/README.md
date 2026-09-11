# WSL → Vast.ai 冻结文件传输

只使用 Python 标准库、已配置的 SSH 别名 `vast-mirheo` 和 rsync。原项目只读；不编译、不安装、不运行项目入口。

在本地 WSL 执行：

```bash
python3 -B /home/lzy/projects/cloud_upload/tools/transfer.py prepare
python3 -B /home/lzy/projects/cloud_upload/tools/transfer.py upload TRANSFER_ID
python3 -B /home/lzy/projects/cloud_upload/tools/transfer.py verify TRANSFER_ID
python3 -B /home/lzy/projects/cloud_upload/tools/transfer.py status TRANSFER_ID
```

将 `TRANSFER_ID` 替换为实际编号。也可先 `source /home/lzy/projects/cloud_upload/current_transfer.env`，再使用 `"$TRANSFER_ID"`。

- `prepare` 只建立新编号；存在未完成传输时拒绝另建编号。源码变化或前置检查失败会停止，保留证据。
- `upload` 复核已有冻结副本后，先 dry-run，再执行一次 rsync。中断后重复同一命令续传，不重新准备副本。已验证编号仅执行只读校验。
- `verify` 重新核对全部 source SHA-256、数量、大小、执行位、链接、包内文件集合及控制文件哈希。重新校验回执存放于工作区 runtime_logs，避免改写完成副本。
- 不使用 `--delete`、`-L`，不修改 SSH 配置或 known_hosts。大小上限 3 GiB；保守预计上传后剩余至少 5 GiB。
- `selection.json` 固定本轮案例来源及少量 runs/build 白名单。排除编译目录，但白名单保留 CMakeCache 作为来源证据，以及原生第三方项目中两个名为 build 的目录内的文本包装源码。
- `source` 中保留原文及 WSL 绝对路径。`metadata/path_mapping.json` 只记录映射。已存在的云端 `.venv` 保留且做前后指纹核查。
- `TRANSFER_VERIFIED_INPUTS_PARTIAL` 表示传输清单已通过校验，仍有云端库构建、路径和运行依赖等未验证项。详见 `metadata/input_dependency_check.json`；历史预算不是云端执行授权。

来源的未提交文件和未跟踪新文件纳入实际文件副本；不复制 Git 元数据。仅传输已落盘内容，无法包含编辑器中未保存的缓冲区。

上传运行日志和命令真实退出码保存在 `cloud_upload/runtime_logs/`；每次传输内另含 `metadata/transfer_log.txt` 和 `metadata/remote_verification.json`。最终只读校验的回执为 `runtime_logs/<TRANSFER_ID>-final-delivery.json`，其中包括最终控制记录已同步的核查。
