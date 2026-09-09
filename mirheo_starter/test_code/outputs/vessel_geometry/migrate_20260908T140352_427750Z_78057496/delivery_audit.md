# 血管迁移交付核查

正式数据包：`/home/lzy/projects/mirheo_starter/data/geometry/migrate_20260908T140352_427750Z_78057496`

HTML：`/home/lzy/projects/mirheo_starter/test_code/outputs/vessel_geometry/migrate_20260908T140352_427750Z_78057496/geometry_review.html`

真实迁移及回读 PASS：73,416 个原始点、67,262 个三角面；所有正式数组和身份精确一致，坐标差 0 m。

合成自动测试 57/57（原测试 40、新迁移测试 17）。真实包单独检查，并以读路径限制验证独立加载。

Windows Chrome 152.0.7977.82 软件渲染交互 PASS；完整三角网格、隐藏管壁、透明度、端口选择、鼠标旋转、滚轮缩放、恢复视角均通过，无远程请求或未捕获 JS 异常。

浏览器临时 UNC profile 出现缓存/数据库/Crashpad 告警，保留原 stderr；未单独测试 Edge。已关闭测试进程并清理临时 profile。

2,386 个受保护文件哈希不变，原血管 Git 状态不变，已安装 Mirheo 库与安装记录哈希一致。原 RBC 和 Hello World 没有重新运行，未宣称其功能回归通过。两个 starter 目录不是 Git 仓库，未进行提交、推送或远程覆盖。

重复执行实际核查命令直接复用：数据包和页面共 37 个文件的 SHA-256、mtime 均未变化。

完整清单、端口面积中心、容差和依赖记录见 delivery_audit.json、implementation_files.json、migration_checks.json 与中文 README。

本次仅完成已验收血管几何及标签向 mirheo_starter 的迁移和核查。
未修改原血管，未生成 SDF、粒子或数值边界，未运行真实血流。
本次迁移的用户人工核查状态为 PENDING。
