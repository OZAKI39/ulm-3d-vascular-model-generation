# Human / BraVa 工程规则

- 唯一日常开发与 Git 根目录：`/home/lzy/projects/new_human_vessel_flow`。先读 README、`docs/CLEANUP_REPORT_ZH.md` 和当前任务。
- 当前自有 Python 实现只在 `src/` 维护，通过普通 import 使用；同一功能不复制到各个算例。
- 长期分支为 `dev/human-vessel-flow`。直接在本目录检查 diff、按文件暂存、commit、push；禁止强制推送、改写历史或覆盖其他项目分支。
- 不为 GitHub 同步创建 clone、worktree、镜像、快照、工程压缩包或拼接运行工作区。一次性小型测试文件应清理，不删除必要验证证据。
- 冻结 CFD、原始网格、checkpoint、原 1565 个微泡的完整状态和账本不得覆盖。大数据留在本地，源码、轻量材料及 SHA256 身份进入 Git。
- 旧群体有 141 条观察截止、2 条数值失败；不能算成功出流，也不能标成新固定灌注协议的结果。固定协议为 1.240815660 s 灌注 + 10 s 观察，dt=0.5 ms、浓度 4e9 m⁻³，尚无正式新大批量结果。
- 新 BraVa 边界设计用中心线阻力估算出口压力；不得改写历史校准场。candidate 0/15 是同一解的展示姿态。
- Gmsh V2/P1 科学精度 FAIL、P2/P1 资源停止、V3 仅提案的状态必须保留。工程测试、守恒、CPU/GPU 一致和动画流畅不能替代物理精度验证。
- 未获明确授权不得启动新 CFD、P2/P1、网格研究、轨迹收敛或正式群体，不恢复已停止队列。诊断目录内原脚本是历史证据，不是自动运行入口。
- 保护 Mouse、Clot、打印、声学、RBC 等独立项目与共享环境。旧原始工程删除需用户确认；缓存/同步副本删除前核对身份、Git、依赖、唯一性和使用进程。
- 未经明确授权不登录远端修改、安装或计算。第三方本地依赖在 `new_vessel_flow_solvers`，不能误删其他项目使用的软件。
- 小测试：`bash scripts/run_small_tests.sh validation/新目录`。独立性：`.venv/bin/python -B -I scripts/check_independence.py validation/新目录`。只读数据：`.venv/bin/human-flow verify --output validation/新目录`。
- 每次交付说明修改、实际测试、临时文件、commit/push 和未解决问题；未执行项写 NOT_RUN。
