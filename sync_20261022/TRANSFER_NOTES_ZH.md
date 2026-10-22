# 发布过程说明（2026-10-22）

本分支为独立项目快照。源文件清单、逐文件SHA、排除原因与兼容链接见本目录；科学数据没有修改，也没有运行任何新计算。

发布时发现旧 VS Code credential helper 已失效，改为只对本次 Git 命令使用现有 gh 登录；没有更改全局认证设置或发布凭据。

原生第三方求解器测试目录中有未取得实际数据的 LFS 指针。GitHub 拒收后，发布副本改为排除这些资源文件，并将原始指针与来源保存在 UPSTREAM_TEST_RESOURCE_REFERENCES.json。全部项目冻结场和正式原始状态使用实际文件，没有指针替代。失败且未发布的本地提交保留，没有改写既有历史。

小鼠目录调整触发 Git sparse pack 对已存在对象的重复传输，HTTP 500 重试后改用 pack.useSparse=false。成功传输为404个对象、2.95 MiB；逻辑工作树大小仍约4.97 GB。该配置仅控制传输，不改变文件内容。BraVa 新原始状态分批提交，每个新增对象批次有界。

最终完整 Git 树、冻结场与远端提交核验回执在本地 github_sync_audits/flow_projects_20261022/DELIVERY_VERIFIED.json，报告为同目录 SYNC_REPORT_ZH.md。此文不内嵌最终自身提交SHA，避免递归。
