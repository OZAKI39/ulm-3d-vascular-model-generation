# 支撑研发输出清理计划（先生成、后执行）

本目录没有 Git 历史，不删除任何算法源码或测试。s1-8 系列降级为研究历史。只删除明确退休支撑目录内可重新生成的输出和 Python/pytest 缓存。所有其他生产、原始数据和共享源码保留；无法确定用途的文件标为 REVIEW_REQUIRED 并保留。

项目文件数：58223；逻辑字节：40334719985。预计删除 7125 文件，释放 6213972861 字节。

- KEEP_SHARED_CODE: 15960
- REVIEW_REQUIRED: 1554
- KEEP_PRODUCTION: 33287
- DELETE_REGENERABLE: 7125
- KEEP_COMPACT_HISTORY: 297

逐文件 KEEP/DELETE/REVIEW、导入引用和 SHA256 见 cleanup_inventory.json；生产保护见 research_retrospective/production_protected_before.json。相同 SHA 的候选副本只在明确退休目录删除，权威生产路径永久保留。旧报告内大文件链接属于历史证据，清理后不代表现存交付物。
