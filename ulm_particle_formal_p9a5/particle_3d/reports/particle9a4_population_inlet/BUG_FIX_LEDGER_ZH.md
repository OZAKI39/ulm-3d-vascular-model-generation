# Bug / change ledger

本轮指定科学逻辑审核中 CONFIRMED_BUG=0、fixed=0；existing scientific files 修改0。[code_change_ledger.json](data/code_change_ledger.json) 是空数组。115 baseline 全部通过才开始新实现。旧 deterministic clock、attempt retry、Method B/C 条件化属于 legacy 模型行为；P9-A.4 通过新增两个模块建立独立合同，没有重写旧模块。

新增 `continuous_infusion.py` 提供counter RNG、Poisson间隔、单次proposal、ordered ledger/checkpoint；新增 `population_inlet_p9a4.py` 显式加载NEW输入、浓度合同与CPU worker。仅复用 TruncatedSonoVue 类，没有迁移或改写其文件，保护旧 imports/SHA。新 tests 22个指定文件共35项，另有conftest及包命名文件。

开发阶段失败和处理（全部保留原日志，不计为已有科学代码 bug）：

| 事件 | 原因 | 最小处理与验证 |
| --- | --- | --- |
| new_initial | 新测试写错 BridgeParticle.radius_m 属性 | 改为原 radius；新35通过 |
| inlet_gate | 新gate将numpy.bool_交给标准JSON编码 | 新辅助脚本显式bool；成功gate与失败partial均保留 |
| smoke_analysis | 新分析脚本误把带identity的cohort字典当事件列表 | 读取events字段；同一已保存结果审核通过，未重跑动力学 |
| final_new_and_legacy | 新tests从全局conftest导入，与其他测试目录重名 | 新目录加包名并用相对导入；组合50项通过 |
| final_full_regression | 原portable工具在子进程临时cwd中解释相对--output，找不到JUnit | 使用绝对输出路径；完整115通过；未改原工具 |
| build_reports | 新报告汇总脚本误读deployment的remote_root键 | 改为实际server_root键；未影响任何计算或数据 |

最后一项是复现工具调用/相对路径可移植性限制，非Particle科学算法bug。本轮将该限制写入复现命令，不放宽科学身份，也不修改受保护旧工具。未执行所有历史模块的任意无效输入、多写者破坏性竞争或跨版本checkpoint迁移；这些边界不宣称已证明安全。

最终165项通过、0失败/0跳过。46,670个发布旧文件、40个原network输入、114个原科学源文件、114个server-bundle源文件及OLD/NEW场哈希全部匹配。92个生成pyc只在portable临时副本的旧源码合同中排除，原文件/合同没有编辑。保护日志：[final_protection_manifest.txt](logs/final_protection_manifest.txt)。
