# 云端单红细胞完整验证

是否完成 Γ=4：NOT_COMPLETE

数值筛选：FAILED

云端库：/workspace/bloodflow/.venv/lib/python3.12/site-packages/libmirheo.cpython-312-x86_64-linux-gnu.so；SHA-256：d45b4fd1b4498b6f365a012e758eec35bc58a7844dabc66f22e617a7b8cf84ee

准备：NOT_COMPLETED；原判据未放宽。

原生正常返回步数：准备 None，剪切 0。

最后保存帧确认的应变下界：0.0；失败的连续调用不把计划步数记为完成。

几何完整帧 46；最大面积漂移 0.01962563474991197，最大体积漂移 0.00791998400044247。

确证自交帧 0；几何硬错误 0。近接触另列。

成员探针：192 点次，确证不符 0。每阶段初态及正常返回末态最多各 192 点；不证明所有粒子严格不渗透。

碰撞计数采样 365 条；容量溢出记录 1 条；原始报错匹配 5 行。没有修改原生容量。

退出码 255；停止原因 None；求解进程墙钟 36.9039248029876 秒。

壁粒子准备 4000 步不计入正式 Γ；setup 已包含壁准备，连续求解耗时包含原生输出，不能重复相加。

GPU 显存来自全设备采样；求解进程墙钟不是 Vast 租赁账单。

静壁到动壁沿用既有两协调器交接：坐标/速度/参考网格保留；ID、分类、oldPositions、力和随机状态重新生成。不是完整状态无缝延续。

DPD 温度目标 kBT=1；原生 Stats 含流动动能，末态 moments 使用扣除分箱平均流速的估计。

物理验证：NOT_MATCHED; material calibration, space/time convergence and HemoCell comparison NOT_TESTED this round；qualified_speedup=null。

下载校验由本地 LOCAL_ARCHIVE_VERIFIED.json 确认；人工验收 PENDING。
