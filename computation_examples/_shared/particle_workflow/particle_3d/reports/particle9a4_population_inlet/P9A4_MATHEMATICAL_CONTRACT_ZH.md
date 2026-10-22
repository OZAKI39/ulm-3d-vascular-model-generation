# P9-A.4 数学与计算合同

设真实入口 S，n 指向血管内，q(x)=(u·n)₊，Q=∫S q(x)dA。尺寸 PDF f 是冻结 SonoVue 在 D≤4 µm 的条件 PDF，C 是该条件群体的 number/m³。稳态源是带独立标记的齐次 Poisson 点过程，强度测度为 `C f(D) q(x) dt dD dA`。

每个 source k：Δtₖ=−log1p(−Uₖ)/(CQ)，tₖ=ΣΔtᵢ；独立取 Dₖ~f、xₖ~q/Q。使用原真实 WALL 和 lower handoff，定义固定可行指示 A(D,x)。active={}，因此 A 不依赖其他粒子或过去事件。在这个假设下，确定性标记筛选仍为 Poisson thinning：

`Q_acc(D)=∫S q(x) A(D,x)dA`

`p_acc=(1/Q)∫f(D) Q_acc(D)dD`

`lambda_enter=C∫f(D)Q_acc(D)dD=lambda_source p_acc`

`p_enter(D,x)=f(D)q(x)A(D,x)/∫f(d)Q_acc(d)dd`

其尺寸边缘正比于 f(D)Q_acc(D)，不是原 f(D)。在失败后重抽位置或尺寸会重新条件化分布；若仍沿用 source CQ 时钟，则改变进入过程。因此一次 source 只调用一次 check，不调用 attempt 或 MethodCSource.event，不删除拒绝行。

正式生成不需要先知道 Q_acc。独立 MC 诊断只解释结果，不改变 source RNG、位置分布或 acceptance。无未来出口信息；无 outlet quota；不把有限尺寸 population split 强行等同 fluid split。

## 圆管解析 benchmark

Poiseuille u(r)=u_max(1−r²/R²)，Q_total=πu_maxR²/2。设 a=D/2、h=2 nm、R_c=R−a−h，x=max(0,min(1,R_c/R))。对半径 0…R_c 积分 2πr u(r)dr，得到 `Q_acc/Q=2x²−x⁴`。

全源径向 CDF `G(r)=2(r/R)²−(r/R)⁴`。固定 D 的 accepted 径向 CDF 在 r≤R_c 时为 G(r)/G(R_c)，随后为 1。离散尺寸概率 pⱼ 的进入边缘为 pⱼG(R_c,j)/ΣpᵢG(R_c,i)。tests 使用 R=2 µm、尺寸 0.6/2.4 µm、概率 0.4/0.6、lambda_source=5 s⁻¹、固定 seed 2594 和 40,000 个 proposals。

数值结果见 [synthetic receipt](data/synthetic_poiseuille_audit.json)：source 径向 CDF 最大误差 0.00495；两尺寸 accepted 比例 0.920753、0.284066，解析值 0.922048、0.293057；条件径向最大误差 0.001420、0.007332；小尺寸 entering 概率实测 0.683363、解析 0.677163；进入率实测 2.682858、解析 2.723267 s⁻¹。固定宽容差 tests 全通过。Fano 只作描述，不作为严格随机 gate。

## 数值实现与可恢复性

RNG 使用 PCG64/SeedSequence(master_seed,source_event_id,role)，独立角色 940/941/942/943 对应 arrival/D/x/orientation；原 STREAMS 顺序不变。worker 只生成独立 marks，主进程按 source ID 顺序 float64 累加时间并分配 accepted particle ID。身份绑定 NEW flow、浓度、源合同、原 114 源码 manifest 与两份新源码 SHA。

完整100k在1/3/6 worker下 ledger SHA 均为 `b89b019874b8d142f31f227ce45ed20f2d0b5301595b1b8dc3dd92cdad3c8cd7`；accepted births 也逐字节一致。这是相同软件/平台数值合同下的保证，不宣称跨任意 NumPy/VTK/CPU 版本必然 bitwise 相同。

checkpoint 保存全部 proposals、下一个 source ID、下一个 accepted ID、累计时间和身份哈希；完成 manifest 最后写入。恢复检查文件哈希并重放有序 reducer，不抽新 RNG。partial、篡改、丢失拒绝、身份不符、重复/缺失 ID、已有输出目录均被拒绝。这里只提供 population checkpoint；不迁移旧轨迹 checkpoint 的源码身份。

负/非有限浓度、非正/非有限正向 Q、rate under/overflow 显式错误。C=0 表示没有下一个 source event；极小正 rate 若超出有限 float64 时间范围则停止报错，不写 Inf。CDF 使用原冻结条件化实现，没有 4 µm 端点裁切堆积。时间精度耗尽会报错而不静默重复正间隔。

## 适用边界

恒定流、恒定条件源浓度、球形、独立单泡、既有 WALL+handoff 模型。P9-A.1 动力学及场插值不变。真实 continuous infusion 浓度、PK/clearance、多泡排斥或随时间变化流场不在本合同内；这些因素可能破坏当前齐次独立 thinning 前提。
