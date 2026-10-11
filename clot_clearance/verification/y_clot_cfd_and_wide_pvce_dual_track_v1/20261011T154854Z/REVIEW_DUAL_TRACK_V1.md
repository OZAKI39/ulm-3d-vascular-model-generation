# Y_CLOT_CFD_AND_WIDE_PVCE_DUAL_TRACK_V1 共同审查

run_id：`20261011T154854Z`。A、B 独立验收，没有单一总科学 PASS。A 的 12 项回答见相邻 A 报告，B 的 12 项回答见相邻 B 报告；下列覆盖共同问题 13–19。

13. **正式算例和资源。** A 共 5/8：低载三态、旧保存态恢复一次、宽薄高载候选一次；B 共 2/2：含栓基线和一次局部体网格加密。没有新动态断裂、无栓重算或第三个网格。A 已记录正式墙钟 89.431 s、CPU 89.038 s，另有恢复墙钟上界 1.564 s；恢复 CPU/RSS/迭代遥测因报告序列化异常不可恢复，不能以零代替。B 两次原生墙钟合计 2081.426 s（34.690 min）。已记录正式墙钟相加 2170.857 s，另加上述恢复上界；不把并行任务墙钟相加冒充端到端时长。完整时长、CPU 下界、内存/GPU 样本与当前报告时刻见 `RESOURCE_AND_COST_SUMMARY.json`。A 所记录峰值 RSS 749.63 MiB；B 峰值采样 RSS 分别 6643.37/6788.71 MiB。15.36 CPU 配额下采用 A 单进程、B 8 MPI，优先 B。4090 沿用原生 PETSc GPU 路径；无压力测试或依赖升级。

14. **实际审查图。** 共 14 组 PNG＋SVG。A1 点阵/控制自由度/尺寸比，A2 周期影像邻居和原生一致性，A3 三个参考切口，A4 低载收敛/能差/G，A5 旧态恢复及宽态停止，A6 未接受高载候选真实键，A7 三个实际低载平衡。B1 原生 Y 与 V3 材料点，B2 真实体网格与开放绕流截面，B3 公共位置速度/压力，B4 实际暴露面压力/剪切，B5 守恒/稳态/压降/力/加密，B6 原生速度场流线与方向，B7 横截面和切向牵引方向。A/B 分别有索引，并有共同 `OPEN_RESULTS.html`。没有新医学外壳、AI 生成图或无分析意义动画。

15. **未完成和停止。** A 低载 G_two/G_three 为 7.140647211281e-06/7.140644132718e-06 J/m²。旧小态高载恢复通过；宽薄高载候选的末十次位置变化及 1 pm 分量差分失败，原冻结门槛保留，没有真正首达阈值区间、高载 A1/A2、材料 Gc。只读 1 nm 差分支持舍入解释，但未重评接受。B 两个解及一次工程筛查完成；表面三角化收敛、接触线极值、到材料点载荷映射、超声 streaming、变形与双向 FSI 不在本轮完成范围。两线按范围停止。

16. **保护情况。** `PROTECTION_AUDIT.json` 对 1097 个事先登记的源码、材料卡、原生库和历史结果重新计算 SHA，结果 PASS。原生环境包版本完全一致，HEAD 未变。新代码只在 `analysis/dual_track_v1`，新数据只在本 run 的独立目录；没有修改 s00、生产载荷或历史接受态，没有 commit/push/reset/clean/stash，也没有停止其他任务。A 的 Y 周期独立能量/HVP 适配仅在隔离诊断目录。

17. **实际 WSL 文档。** `/home/lzy/projects/clot_rupture.md` 状态：`PASS_WRITTEN_AND_READ_BACK`。更新保留 ①–⑨、历史无栓背景流使用记录及各轮负面科研结论，新增含栓解和宽薄 PVCE 的已验证结论；不会把旧碎裂或输运轨迹改称使用了新含栓流场。

18. **备份、diff 与 SHA。** `document_edit/before.md`、`document_edit/after.md`、`document_edit/change.diff`、`document_edit/AUDIT.json` 记录最新原文、并发检查、写入/读回和前后 SHA；任务开始备份另在 `provenance/clot_rupture.before.md`。`SOURCE_MANIFEST.json`、`CODE_ADDITIONS.diff`、`PROTECTION_AUDIT.json`、`provenance/ENVIRONMENT_BEFORE.json` 与 `ENVIRONMENT_AFTER.json`、`REPRODUCE.txt` 提供源码/库身份、代码差异和条件复现说明。`SERVER_ARTIFACT_MANIFEST.json` 清点本轮服务器全部原始数据；`COMPACT_MANIFEST.json` 与 `LOCAL_RETRIEVAL_AUDIT.json` 记录 WSL 精简回收及 SHA 验证。完整 CFD 网格、原生全场、逐步状态和 LAMMPS 原始输入均保留服务器，不降低输出精度。

19. **下一轮唯一优先项。** 验证含栓 CFD 暴露面牵引向现有 LAMMPS 材料点的离线保力、保力矩映射，先核对单位、坐标、面身份和面积权重，再决定是否接入加载。 理由是 B 已得到含栓场及稳定的积分载荷，仍缺少到实际材料点的受力接口验证；A 未识别 Gc 不否定 B 的独立成果。建议只做保力/保力矩离线检查，不自动启动生产加载、材料升级、AI 大扫描或双向耦合。

加密 CFD 代表值：最大速度 10.745852134 mm/s；平均表压 0.516522370 Pa；平均/P95 切向剪切 0.016769856/0.037783332 Pa；总力 [5.9221682140260414e-08, 2.278610279267278e-08, 6.607382089731738e-07] N，模长 0.663778122 µN。压力是出口参考表压，力只含暴露流体面，不含壁面粘附反力。
