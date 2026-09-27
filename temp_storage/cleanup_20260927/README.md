# 本次 WSL 清理回执

已按 2026-09-27 用户指令删除旧工作树、历史开发内容及本地服务器设置副本。

- `result.json`：空间释放和最终检查结果；主清理约 39.1 GiB，另有少量收尾删除。
- `deletion_plan.json` / `deleted.json`：主批次的具体删除路径、原因和计数。
- `additional_deletions.json`：收尾删除及缓存清理。
- `retained_before.json` / `final_retention_verification.json`：清理前后的保护哈希。四个核验脚本按清理后的保留范围更新，一个全历史保护清单被当前科学数据清单替代；其余 8,290 个保护文件保持原字节。
- `science_verification_final.json`：7,760 个保留科学文件、500 条正式微泡、500 条点示踪核验通过。
- `current_checks/`、`network_final.xml`、`rbc_generation.xml`：131 项针对性检查通过。两个耗时的完整轨迹重算测试未在此次清理中执行。
- `rbc_style_verification.json`：从旧渲染器提取的配色与相机函数 AST 完全一致。
- `html_links.json`：保留的 6 个 HTML 入口均没有失效本地链接。

这是删除回执，不是旧数据备份；删除内容未转移到其他 WSL 归档目录。服务器、独立血管打印项目、当前共享环境和 Git 元数据保留。
