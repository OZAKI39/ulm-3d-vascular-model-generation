# Transient artifact policy — prospective, pending review

本文件与 `SCIENTIFIC_PROTECTED_MANIFEST.json` 是完成 123 项分类后的新基线提案。**本轮恢复仍被旧测试依赖阻塞；没有批准继续清理、修改旧测试或替换旧 manifest。**

## 永久保护的内容

生产及候选科学源码、测试、配置、接受态、解释结论所需的诊断/失败态、最终科学 JSON/CSV、REVIEW、STATUS、科学图件、冻结源码快照和关键 provenance 均按 SHA256 保护。新清单采取保守范围；诊断态不能仅因未被接受而删除。

`results/pfpd_at1_c5_stability_v1/cache/` 中供当前测试和报告读取的物理 Hessian 等数据属于关键科学证据，继续保护。“cache” 目录名不决定数据能否清理。

## 今后不自动纳入永久保护的内容

可再生 Python/Numba/Matplotlib 缓存、未被最终结果依赖的 Hessian 分列中间件、Jacobian 装配中间件、worker scratch、下载临时包、pytest 临时目录，不应由递归 PRE/POST 清单自动变成永久科学基线。

只有在依赖审计证明没有实际读取者、没有唯一科学信息、已有可验证的最终替代物之后，未来单独获授权的清理才可处理这些项目。分类为 `UNKNOWN` 的条目必须保留。临时包中的唯一来源文献、版本 pin 或验证证据必须先单独保存。

## 历史兼容性约束

旧测试通过 manifest 遍历 `read_bytes()` 并比较 SHA256 时，对被读取文件形成真实测试依赖。此次 123 项全部具有这类依赖；其中 72 项尚无原字节，不能因为科学计算可再生就标成 `MANIFEST_ONLY` 或 `no-recovery-needed`。

所有旧 PRE/POST manifest 原样保存。新清单不回写旧清单，不修改旧测试断言，不隐瞒遗留失败。缓存再编译可能产生不同序列化字节；科学等价的 JSON 也不能满足旧 SHA 断言。

## 后续审核条件

当前只执行恢复审计和报告。审核者需先决定如何处置 `CRITICAL_UNRECOVERED.csv` 的未恢复依赖，再单独批准任何历史兼容性迁移。未获批准前，停止 repository cleanup、activation audit 和新的科学计算。

若未来开展迁移，必须同时保留旧清单、缺失异常记录、旧测试失败证据以及新保护清单的来源；不把迁移写成“旧完整性测试已恢复通过”。
