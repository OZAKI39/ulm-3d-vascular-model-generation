# 直接在当前项目中开发

唯一开发和 Git 根目录为 `/home/lzy/projects/new_mouse_vessel_flow`。当前分支为 `dev/mouse-vessel-flow`，远端为 `https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git`。

本分支以已验证的新工程建立独立提交历史。首个提交 `e68fe08` 保存清理前的新工程代码，包括后来退休的一次性迁移脚本。旧同步分支 `sync/mouse-vascular-flow-project-20261022` 保持原样，已核对远端提交为 `987d0638944f09e1acc1f415b26bf5b0479c3f1d`。两者布局不同，没有把旧快照检出覆盖新文件，也没有重写任何旧分支。

本地 `.git` 是此目录内的独立对象库，没有外部 `commondir` 或 `alternates`；只有本项目一个工作树。为核对历史仅取回旧同步分支的浅层提交/树元数据，没有另建 clone、镜像或数据副本。旧源码来源对应表见 `audit/source_map.json`；必要时可从原远端历史读取，而不是再复制工程。

日常流程：

```bash
cd /home/lzy/projects/new_mouse_vessel_flow
git status
git diff
.venv/bin/python -m pytest tests -q
git add -- 本次修改的文件
git diff --cached --stat
git diff --cached --check
git commit -m "说明具体修改"
git push
```

Git 管理源码、配置、测试、依赖说明、必要的小输入、结果身份清单和精选验收报告。冻结流场、原始轨迹、checkpoint、完整数据和旧逐文件审计仍保留在本地；不因 Git 忽略它们而删除。`configs/active_data_manifest.json` 记录这些数据的 SHA256 和来源。**新电脑仅 clone 代码不能自动获得所有科研数据**，应按清单另行迁移并校验数据，再运行带真实数据的检查。

当前分支首次整理后的最新提交、远端一致性与 upstream 见 `cleanup_20261023/GIT_STATE.json`。该记录描述已验证的提交，不递归嵌入包含自身的提交号。共享旧 Git 库仍服务其他项目，未清空对象、删除其他分支或回收历史。

不再使用 GitHub 同步工程：本项目四个旧迁移/整理脚本已从工作树移除，可在 `e68fe08` 查阅；共享历史导出/发布工具中只禁用了小鼠路径，BraVa 路径未改动。没有发现小鼠同步定时任务。根目录 `AGENTS.md` 明确约束以后的 Codex 操作。
