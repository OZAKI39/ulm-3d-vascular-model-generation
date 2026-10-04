# 当前 BraVa 开发进度同步清单

目标：`https://github.com/OZAKI39/ulm-3d-vascular-model-generation`，分支 `sync/roi-flow-mb1500-wss-stop-analysis-20260928`。本轮不向 main 提交、不 force push、不改写历史、不启用 LFS；仓库根 .gitattributes 保持原字节。

## 当前权威科学状态

**REFERENCE_PIPE_FIELD_FAIL_SCIENTIFIC + MIXED_FIELD_ERROR + P1_FIELD_TRAJECTORY_IMPACT_LARGE（DIAGNOSTIC_ONLY）**。

GPU svMultiPhysics/PETSc production 路线已在冻结Gmsh V2上完成steady，五个区间20→30→40→50→60→70通过，第71步正常保存，318次linear solve全部收敛、exit0。耗时5334.385s。局部场精度gate仍FAIL。

0.377/149.195µm支持P1表示限制，29.839/59.678µm更支持CFD−P1离散解差异；向量误差有负cross-term，不能用百分数直接相减。四sentinel exact-vs-P1诊断由A局部10µm通过时间差5.148892%触发LARGE，四组无径向分离、gap/contact类别改变。这不是正式CFD-vs-oracle或dt认证。

V3=PROPOSAL_ONLY/NOT_RUN；正式CFD-vs-oracle、0.5/0.25/0.125ms收敛、新BraVa CFD、新8/46轨迹验证、新population均NOT_RUN。本轮科学重算数为零，NO PRODUCTION PHYSICS CHANGE。

## 来源、分支与提交边界

本地独立算例：`/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin`，清理前HEAD `91ebc02162e7b0c774cf649b5669bcc5c8026e54`，无remote。实际GitHub同步checkout：`/home/lzy/projects/computation_examples/_archives/github_sync_20260928_latest`，原dev HEAD `0b808fb6c90e27674d79078c49358ced2e9258ff`；fetch后目标远端 `94c644686092158d7f9618c00a6e30cbb10f7f35` 是其祖先，ahead11/behind0。本次把目标本地分支fast-forward到已有dev，然后增加一个当前清理/同步提交；保留全部已有历史。算例自身另留一个本地清理提交，独立历史不推到GitHub。

映射：算例→GitHub `brava_flow_roi_18mlmin/`；共享粒子→`ulm_particle_formal_p9a5/particle_3d/src/`与`tests/`。144个源码/资源文件和189个测试文件逐字节一致；仅记录校验，不整库stage共享工作树。不触碰其它小鼠/RBC/几何项目源码与结果。

## 本次包含

- 当前生产particle/population与GPU backend、fixed-window协议，以及所有active tests。
- `scripts/`、`vendor/flow_solver_support/`、GPU环境wrapper、MPI launcher、RA header/policy、PETSc build recipe、host-sync patch及solver identities；不提交binary/runtime/第三方完整checkout。
- 冻结V2生成/转换/strict containment代码、exact cylinder、25/50/75µm shell/probes、cap correction、plan与QA。
- 最终V2原始MSH/solver mesh、GPU `result_071.vtu`、`stFile_071.bin`、12.26MB `FIELD_SAMPLES.npz`、960 probes及全部decomposition/particle diagnostic NPZ/JSON。
- 当前总报告、mesh/field/trajectory/dt报告、FIELD_VALIDATION、GPU identity、decomposition、sensitivity、V3 proposal。V3仅提案，不冒充已生成。
- 最终三张PNG及简洁HTML入口；原1565结果账本/身份/统计和最新四视角。旧结果仍包含141观察截止与2数值失败。
- 原运行execution.json、steady intervals、GPU compact evidence。2.82MB raw solver.log及195KB/5秒GPU samples是未改动certifier的直接依赖，作为必要identity证据保留。

逐文件SHA/size见 `GITHUB_SYNC_INCLUDED.json`；同步清单本身不嵌入自身哈希，最终Git tree负责完整身份。原绝对路径配置保持字节，不将此快照宣称为空白机器一键部署包。

## 清理与排除

净删除1631个文件/链接，94,766,770 bytes；另移除supervisor逐秒samples，节省2,414,984 bytes。包括V1/cap旧图和CPU monitor、V1冗余mesh、GPU10–70重复VTU/checkpoint、旧CPU瞬态VTU/重复checkpoint/raw log、PDF/SVG副本、pytest/Python缓存及GitHub旧第三方checkout。MMG/Omega/refine/MeshSim/P2/BC-B活动实现此前已删除，本轮不恢复；保留最终结论。

旧V1 BREP、端面、RAW_MESH、BOUNDARY_TOPOLOGY，以及CPU step10 checkpoint仍有当前代码依赖，保留并明确历史角色。UNKNOWN全部保留。引用复查发现3个历史删除路径已被当前权威报告复用，已恢复本轮开始时的准确字节并保护；旧名单不能据文件名直接重用。执行与复查见 `production_workflow/validation/sync_cleanup/`。

`GITHUB_SYNC_EXCLUDED.json`列出本地保留但本次未复制的文件；`EXCLUDED_LARGE_ARTIFACTS.json`列出其中>=5MB，20项>=100MB全部禁止普通Git上传（未压缩、未切片、未启用LFS）。主要为历史全域field/audit arrays、全壁CSV/NPZ和旧GPU cache。原始population accepted states和逐轨迹派生副本不整体加入Git；本地/服务器保留，账本/summary/identity保留在Git。复现优先取回原文件并验SHA，不自动重跑任何停止路线。

只读服务器核验759个大文件路径：758存在且size相符；20项>=100MB的服务器SHA与本地全部一致。1个旧流线动画副本仅本地存在，清单明确LOCAL_ONLY，不声称服务器可恢复。没有删除/修改服务器数据或服务。

发布后本地 `GITHUB_SYNC_RESULT_ZH.md`（或失败时FAILED报告）是操作回执，记录实际final SHA与push结果；为避免提交含自身SHA的递归问题，它明确作为本地排除项，不作为科学结果。GitHub中的本清单不提前声称push成功。

## 实际检查

38项快速unit tests通过；strict containment、serialized output、Gmsh几何gate、GPU identity守卫、解析gap/contact、向量cross-term与近壁权重测试，不运行trajectory。生产/解析reference/field-analysis/Gmsh generator/solver-support imports通过。Local Torch未安装，三个GPU backend import在既有服务器环境通过；142个对应服务器源码与本地一致，无kernel编译/计算。冻结mesh/field/checkpoint及15项field source SHA通过，40个contract core源码与本地/GitHub一致；科学JSON字节未改。

push前另检查diff、秘密模式和新增对象>=100MB限制；只记录扫描结论与路径，不泄露敏感内容。前后体积、LOC、tracked/untracked及最大20项见 `production_workflow/validation/sync_cleanup/WORKTREE_METRICS.json`。该统计不含.git数据库，清理腾出空间与新增同步文件体积分别计量。

push后fetch验证local HEAD与指定remote branch一致，然后停止。不得自动生成V3、运行CFD/particle或进入新BraVa。
