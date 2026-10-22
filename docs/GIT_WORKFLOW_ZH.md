# Human 工程的直接 Git 开发方式

唯一工作目录为 `/home/lzy/projects/new_human_vessel_flow`。本目录的 `.git` 独立保存对象与索引，不依赖旧共享 Git 库或外部 worktree。

长期分支：`dev/human-vessel-flow`；远端：`https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`；上游：`origin/dev/human-vessel-flow`。

旧 `sync/brava-flow-project-20261022` 是历史发布快照，继续保留在远端。新工程目录结构已经重组，因此新开发分支从本目录的独立根提交开始；旧源码来源、commit 与 SHA 由 provenance 记录对应，不伪造为同一目录历史。没有修改 main、Mouse 或其他分支。

正常开发直接修改 `src/`，通过普通 import 调用。查看真实差异后，按文件提交：

```bash
cd /home/lzy/projects/new_human_vessel_flow
git status
git diff
bash scripts/run_small_tests.sh validation/一个新的测试目录
git add src/具体修改文件.py tests/相关测试.py
git diff --cached --stat
git diff --cached --check
git commit -m "说明具体修改"
git push
```

只运行与修改相关的检查；上面的完整小测试不表示每次文档修改都必须重跑。冻结数据只读检查使用 `human-flow verify --output validation/新目录`。不要覆盖已有验证目录。

Git 保存当前源码、配置、必要原生源码和补丁、轻量输入、报告与文件身份。大型网格、场、微泡状态、动画与缓存留在本地；`configs/active_data_manifest.json` 给出保留资料的身份。clone 开发分支不会自动下载全部本地科研数据。

不要创建日期命名的同步分支、完整 clone、worktree、镜像或工程压缩包来完成普通推送。旧共享同步入口已经退出并指向当前目录，不再自动复制或发布。完整 CFD、GPU 群体与 P2 科学验证不是 Git 同步步骤。

初始提交 `4a2d4e7` 保存已整理工程以及本次移除的一次性迁移工具。历史脚本仅用于追溯，不能据此恢复旧队列。需要追溯单文件时可用 `git show 提交:相对路径`，无需另建运行工程。

2026-10-23 经用户确认删除原 BraVa 工程前，原仓库的 4 个分支、全部 37 个历史/reflog 提交及未挂接对象已保存在本目录 `.git` 的 `refs/archive/brava-original/` 与 `refs/archive/brava-unreferenced-objects/`。这些是本地历史引用，不是额外工作树，也不参与普通 `git push`。完整历史约占 1 GB；远端原发布分支不改写。详情见 [历史保留凭据](cleanup_confirmed_20261023/history_preservation_receipt.json)。

例如可直接查看原文件：`git show a4b238b:ACTIVE_DEVELOPMENT_SCOPE.md`。重新 clone 当前开发分支只会取得已发布的当前工程，不包含上述仅本地保留的原仓库历史或大型科研数据。
