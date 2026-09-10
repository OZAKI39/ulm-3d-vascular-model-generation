# 单细胞基准归档校验

分支 `sync/hemocell-mirheo-single-rbc-benchmark-20260910T120318Z`，来源 campaign `rbc_shear_20260910`。基础分支 `sync/hemocell-mirheo-benchmark-20260909T220959Z`，基础提交 `bcac1d7cf8b6c596f663d047523febee82082b0a`。

- [仓库总入口](../../BENCHMARK_REVIEW.md)
- [完整文件清单](SYNC_MANIFEST.json)：源文件、目标位置、字节数、SHA-256、执行位及复制/继承身份。
- [相对路径映射](PATH_MAPPING.json) · [排除内容](EXCLUDED_FILES.json) · [大文件审查](LARGE_FILE_REVIEW.json)
- [原有 Git 状态](SOURCE_GIT_STATE.json) · [非凭据匹配审查](SENSITIVE_SCAN_REVIEW.json)
- [本次同步校验](SYNC_VALIDATION.json) · [只读校验脚本](verify_archive.py) · [轻量结果字段索引](RESULT_INDEX.json)
- [完整中文说明](../../mirheo_starter/test_code/README_single_rbc_benchmark.md)

在下载的分支中运行：

```bash
python3 sync_reports/20260910T120318Z/verify_archive.py
```

仅在原 WSL 机器上可加 `--check-sources` 比较源文件。脚本只读哈希、JSON/CSV/NPZ、Python 语法、运行封存及相对链接，不导入求解器，不重跑物理测试或浏览器。同步校验 PASS 不代表科学筛选通过；原结论仍为 PARTIAL / FAILED，速度比 null、人工 PENDING。

原始 HTML 和完整 results.json 保留精确字节，并核对以前的最终 Chrome 报告。GitHub HTML 页面用于查看源文件，应下载后本地打开。全部成功/失败记录与已有浏览器排版失败记录都保留；判断最终页面请使用 `browser_delivery_verified/browser_check.json`。

只排除 3 个运行锁以及本地编译程序/对象。所选 campaign 的数值原始数据、NPZ、分析缓存和图均保留。虚拟环境、完整 vendor 树、下载缓存、认证资料及无关项目不新增上传；所需历史依赖和版本、补丁、许可证由基础分支继承。原配置绝对路径不批量替换，路径映射不保证未安装环境的机器可直接执行求解。

`SYNC_MANIFEST.json` 与 `SYNC_VALIDATION.json` 不递归对自身或对方计算哈希；由本次 Git 提交确定身份。其余派生说明及脚本均纳入清单。推送后远程 SHA 和独立取回核查写入同步工作区 `.git/sync_work/REMOTE_VERIFICATION.json`，不反复 amend。

没有重装、编译、求解、重新标定、重新科学验证、原始证据修改、默认分支修改、旧分支覆盖、强推、PR 合并或 Pages 发布。
