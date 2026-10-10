# MEVO 项目同步快照 · 2026-10-05

> 本文件是 2026-10-05 的历史快照。当前 F5 基线、新服务器环境及
> 2026-10-10 本地备份删除状态见
> [当前项目入口](docs/CURRENT_PROJECT_STATE_20261010.md)。以下旧运行状态不代表当前状态，也不构成恢复旧任务的指令。

本目录从正在使用的 `mevo_acoustic_microbubbles` 工程同步，源提交为
`a6c3153f484113d37ffb5e9b4b57f0c6594ef396`（`stage4-bifurcation-control`）。
原有 278 个 Git 跟踪文件按原字节保留；没有修改模型、配置、轨迹或原始结果身份。
这次提交不运行模拟，也不恢复已暂停的后续开发控制器。

## 从这里阅读

- [原工程说明与历史阶段](README.md)、[Stage 4 运行说明](docs/STAGE4_QUICKSTART.md)
- [当前运行状态快照](sync/RUNTIME_SNAPSHOT.json)：控制器暂停，既有 J 子进程仍运行
- [Stage 4I.5 协议](docs/STAGE4I5_PROTOCOL.md)、[后端决策](results/stage4/performance_refactor/backend_decision/STATE.json)
- [新版 B2/B3 动画播放页](results/stage4/visualization/B2_B3_clean_enlarged_v7/OPEN_RESULTS.html)
- [3000 ms 后群内重排核查](results/stage4/visualization/B2_B3_clean_enlarged_v7/shape_audit/SHAPE_AUDIT.html)
- [导出内容及排除规则](sync/EXPORT_POLICY.json)、[逐文件哈希](sync/FILE_MANIFEST.json)、[验证报告](sync/VERIFICATION.json)

下载或克隆后在浏览器打开 HTML。GitHub 文件预览不会直接执行本地 HTML 播放器。
B2/B3 原始保存状态、声场、血管表面和最新动画均保留；球群内部显示放大 60 倍，
动画主场景比 v6 再放大约 12%，1080p、60 fps、全程 0.30×。没有平滑或重写科学轨迹。

## 结果身份与未完成部分

Stage 3 状态是 `MINIMAL_SIMULATOR_VALIDATED_WITH_SCOPE`。
正式 Stage 4I B2/B3 都来自 `results/stage4/primary_secondary/20261003T132343245938Z/`，
每例 256 个初始微泡、6 s 观察、301 个保存状态。
B2 为 A=76、B=0、观察截止未分配=180；B3 为 A=0、B=88、观察截止未分配=168。
未出流的泡不能计为成功。近壁、近接触与有限雷诺数限制仍然成立。

Stage 4I.5 的已完成比较推荐 `P1_ordered` 用于 J/K。Stage J 尚未完成，
其活跃运行 `20261005T032931_171586Z_B2_seed20261002` 的 HDF5 和接触事件流没有上传。
保留配置、源码快照、身份与进度；这些记录只是文件级观察，不代表阶段验收通过。
原始历史状态文件可能保留 `RUNNING`，应结合 `sync/RUNTIME_SNAPSHOT.json` 的进程观察解读。

## 保留和排除

保留源代码、配置、测试、依赖锁定、许可证、Stage 1–4 的轻量日志/报告，
正式统计实验轨迹、响应表、关键声场/CFD 输入、校验数据、源码快照和最新动画。
`sync/server_logs/` 额外保存本项目的 65 份 Supervisor 日志快照。
`sync/SOURCE_HISTORY.bundle` 保存原独立工程全部分支/标签的完整 Git 历史，约 0.9 MB。

排除数 GB 的逐接触事件流、可重建的响应节点缓存、重复性能回放、
旧动画、逐帧 VTK 导出、虚拟环境、下载包/编译产物和活跃 J 二进制文件。
所有排除文件的相对位置、观测大小和原因在 `sync/EXCLUDED_FILES.csv.gz`。
本机和服务器原始数据没有因导出被删除。

**这里是精选数据快照，不是完整 gold 归档。** 原 gold manifest 保持原哈希，
但其两份大型 `events.jsonl` 未上传，因此不要对这个精选副本执行完整 gold 封存验收，
也不要重新封存、覆盖或把它称作完整的 16 GB gold 备份。
参考风格视频属于另一个 clot-clearance 项目，Arial 为系统字体，二者不随本快照上传；原路径与哈希见 `sync/EXTERNAL_REFERENCES.json`，不影响直接播放交付 MP4。
历史页面可能链接到排除项；以本同步入口、文件清单和最新 v7 播放页为准。

## 校验、恢复源码和使用环境

在本目录运行以下命令，只读取已有数据，不启动任何模拟：

```bash
python3 sync/verify_snapshot.py
# 同时检查正式 B2/B3 HDF5 的结果身份，需要 numpy/h5py：
python3 sync/verify_snapshot.py --scientific
```

若要恢复独立项目完整源码历史，在新目录执行：

```bash
git clone -b stage4-bifurcation-control sync/SOURCE_HISTORY.bundle /path/to/new-mevo-source
```

本目录保留原始 `pyproject.toml` 和 `requirements.lock.txt`。
一般 Python 环境可使用 `python -m pip install -e '.[dev,visualization]'`；
CUDA、k-Wave、APECSS 与 SimVascular 的环境和范围见原 README、第三方许可证及各阶段文档。
虚拟环境、可执行文件和第三方源码需按原文档安装，不包含在仓库中。

原始配置、来源清单、渲染字体及部署脚本有 `/workspace/...`、`/home/lzy/...`、
`/mnt/c/Windows/Fonts/arial.ttf` 等历史路径，保留这些内容是为了不改变身份哈希。
查看结果与上述哈希校验不依赖这些位置；重新运行部署或渲染前需在副本中配置本机路径和字体。
旧命令中的 `submit`、`run_post_candidates` 会启动计算，不能把它们当作只读校验命令。
本同步不改动 BraVa 算例、旧兼容链接或另一个 clot-clearance 项目。
