# 新聊天入口：microbubble 与墙模型进度

先读 [当前状态](CURRENT_PROJECT_STATUS.md) 和 [机器可读状态](CURRENT_STATE.json)，再读 [最新固定小球墙报告](fixed_multiblob/FIXED_MULTIBLOB_WALL_FEASIBILITY_AUDIT_REPORT.md)。本目录在此前 HemoCell/Stage4/RBC 交接之后新增，保留旧分支内容。

研究目标仍是 microbubble 的轨迹、速度、近壁停留及未来黏附。RBC 是后续背景流体环境；本批 microbubble 算例中 RBC、adhesion、buoyancy、lift 和 two-way coupling 均未启用。不要自动进入下一阶段。

已完成 SonoVue number-weighted 连续分布采样、最小 LAMMPS 技术引擎、Palabos 冻结速度场单向耦合和 overdamped 被动输运。冻结场来自新介质 5000 步 smoke，不能写成新介质已充分收敛的正式长期解。稳定 LAMMPS 为 22Jul2025 Update 6，提交 9c5ab448c78a14fd534619622162ba418d6a1fb1；2Sep2026 候选源代码门槛失败，未成为当前安装。

当前 pair 近场模型是 rigid no-slip sphere 技术假设；normal、tangential、transverse rotation、TR 和硬非重叠已验证，twist 尚缺，不能宣称完整 rotational lubrication。它尚不包含变形、壳层或可变滑移 SonoVue 物理。

墙模型现有两条证据：

- RMBW lookup/local-plane 的代码与合成算例通过且参考限制保留，但真实 vascular STL 的局部平面有效性仅 4/10000，均为平端盖。真实 8 泡 Case J 未接受任何时间步，物理时间为 0；问题不能归结成 lookup table 实现失败。
- 最新 Pecnut fixed hard spheres 平墙审查得到 FAIL_WALL_REPRESENTATION。最优已测诊断 β=.125、HEX、2 层的核心最坏矩阵作用误差仍为 88.678%，门槛 5%。细化改善了相位敏感性，但未得到合格连续墙近似。结论限定在名义平面、已测分辨率和资源上限；真正 regularized multiblob 后端尚未运行，不能一起否定。

最新审查保留 5266 个有效 6×6 查询、3 个失败查询及数值性质核查。2,000 墙球装配完成，查询失败和进程异常均保留，根因未定；β=.0625/.03125 未通过资源预检。曲壁和真实 STL 条件阶段未运行。人工视觉审阅仍 PENDING。

本次 Git 同步只归档，没有运行求解器、改物理合同或替换正式源代码。下一阶段仅建议 BEM_CURVED_WALL_FEASIBILITY_AUDIT；需用户另行授权。历史阶段 README 中的 NEXT_STAGE 不代表当前执行指令。

数据范围与复现限制见 [README](README.md)。最新版原始矩阵在 fixed_multiblob/remote_raw/；约 543 MiB 的派生汇总 HDF5 没有进入 Git，但原始矩阵、冻结输入、finalizer 和图表已包含。全部遗漏文件的身份和本地位置见 OMITTED_ARTIFACTS.tsv。
