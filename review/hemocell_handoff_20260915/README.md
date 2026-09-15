# HemoCell Stage4 → RBC Stage1 科研交接

入口：[NEW_CHAT_CONTEXT.md](NEW_CHAT_CONTEXT.md) → [CURRENT_PROJECT_STATUS.md](CURRENT_PROJECT_STATUS.md) → [CURRENT_STATE.json](CURRENT_STATE.json)。当前 **RBC Stage1 FAIL_GEOMETRY_GATE，RBC_TIMESTEPS=0，RBC wall interaction ABSENT**。下一阶段为 **RBC_GEOMETRY_AND_WALL_COMPATIBILITY_STAGE**。

## 范围与导航

| 目录 | 内容 |
|---|---|
| `pure_fluid/` | 旧介质 Step3C 正式长验证终态与独立审计 |
| `gpu_stage4/` | 原 RTX5090 Stage4 与新 RTX4090 恢复的紧凑证据 |
| `new_medium/` | PBS/BSA 生成合同、500/5000 步 smoke 数据与终审 |
| `rbc_stage1/` | 静态几何失败、构建、合同、独立复核及两个小型 VTP |
| `contracts/` | 冻结 Step3C 与 Stage4 合同；介质/RBC 合同保留在各自报告目录 |
| `source/` | 按实际执行 provenance 选择的项目源码、生成器、构建和评价工具 |
| `patches/` | Stage4、metadata、介质适配及原 HemoCell 兼容性补丁 |
| `provenance/` | 来源、独立核验、状态差异、审查限制和未上传数据身份索引 |

已有 `review_bundle/step3_review/` 与所有现有源文件不变。新目录是 Git 审查交接，不是免安装的可运行包；工具链、upstream、binary、长程全场数据均未携带。历史脚本保留实际使用的工作目录，后续恢复必须根据 source/README.md 核验和显式映射，不能直接在本交接目录执行。

## 保真与哈希

`SOURCE_MANIFEST.tsv` 列出本目录所有文件的来源、大小和角色。原件哈希及导出是否变换见 `provenance/ORIGINAL_FILE_HASHES.tsv`。科学源码、数值合同、报告均保持原件字节；本次新文档独立生成，未把旧报告的 PENDING/NEXT_STEP 改成新结果。

运行 `sha256sum -c SHA256SUMS` 可校验仅本目录的全部文件（排除校验表自身）。`SOURCE_MANIFEST.tsv` 自身的实际哈希在 SHA256SUMS；其自记录用 `SEE_SHA256SUMS` 标注。SHA256SUMS 自身用 `EXTERNAL_SEAL_ONLY` 标注，完整文件哈希保存在本地发布回执并由 Git commit 固定。两行的大小仍为实际字节数。这样避免伪造循环自哈希；其他文件记录普通 SHA256。

未上传的大型对象见 `OMITTED_LARGE_ARTIFACTS.tsv`。集合行的 sha256 明确是对应压缩文件身份索引的哈希，**不是目录内容的假想单文件哈希**；大文件行是原文件 SHA256。索引包含绝对本地保留路径、逻辑大小和原有 manifest 身份。集合行与大文件行重叠，不可直接相加，硬链接也会影响逻辑大小与实际磁盘占用的关系。范围见 `provenance/OMISSION_INDEX_SCOPE.json`。没有删除原数据。

局部 `.gitattributes` 将小型 VTP 保存为普通 Git blob，避免继承根目录 LFS 规则；不修改根属性。它保留原件换行，并将历史 CSV 的 CRLF 识别为合法行结束。仅 unified patch 的必需上下文空白、原版权文本空白，以及唯一 numerics generator 的既有 EOF 空行免除 whitespace 警告，保持合同记录的原始 SHA256。其他源码与新文档接受正常 `git diff --check`。

## 科学解释边界

旧介质 formal Step3C PASS 不代表 PBS/BSA multiplier 已长期校准；PBS/BSA 是未实验测量的开发假设。RBC 名义 50 µm³ 与实际 45.046330078 µm³ 分开记录。失败姿态不允许作为 spawn；0 timesteps 不等于 runtime safety PASS。原生库中的 API 存在不代表当前 Guo 应用已启用 RBC-wall repulsion。HUMAN_RBC_REVIEW 仍 PENDING。

本次仅归档、校验、单一新分支 commit/push；无编译、模拟、PR、merge、force push 或旧工作树修改。
