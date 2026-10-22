# 小鼠血管工程规则

- 唯一日常开发和 Git 根目录：`/home/lzy/projects/new_mouse_vessel_flow`。先读 README 与 `docs/CLEANUP_REPORT_ZH.md`。
- 只修改 `mouse_vessel_flow/` 中的权威实现；不同算例用 `cases/*/case.json` 区分，不复制求解程序。
- 直接在此目录检查 diff、测试、按文件暂存、commit 和 push；长期分支为 `dev/mouse-vessel-flow`，不得推送 `main` 或强制改写历史。
- 禁止为普通 GitHub 同步建立完整 clone、worktree、镜像、快照或打包工程。一次性小文件应在任务结束时清理。
- 冻结 CFD、正式轨迹、原始输入、checkpoint 和验收证据不得覆盖。大数据留在本地，身份清单进入 Git；不切片绕过大小限制、不擅改 LFS。
- 不自动修改 BraVa、打印、血栓、声学或其他独立项目；共享环境、共享 Git 库不属于小鼠可随意删除的缓存。
- 删除前核对真实路径、Git 归属、未提交/未推送修改、唯一性、依赖和保留副本。原始工作区或用途不明内容列清单等待用户确认。
- 未明确授权不得启动新 CFD、正式微泡批次或大型实验。已取消细档血管/O2 敏感性保持 `SKIPPED_BY_USER`。
- 小检查：`.venv/bin/python -m pytest tests -q`；只读全核验：`.venv/bin/mouse-flow validate --case all --full-integrity --output validation/新的目录名`。
- 修改受身份约束的源码后，记录原因和差异，再更新 `configs/runtime_source_identity.json`；不改原始科学身份。
- 每次交付说明修改文件、实际测试、临时文件是否清理、commit/push 状态及未解决事项。未执行项目写 `NOT_RUN`。
