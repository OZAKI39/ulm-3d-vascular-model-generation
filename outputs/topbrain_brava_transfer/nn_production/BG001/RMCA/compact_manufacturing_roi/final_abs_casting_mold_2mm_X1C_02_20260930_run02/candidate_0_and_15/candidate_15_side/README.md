# Candidate 15

侧放：+X 侧壁贴打印板，开口朝侧方。

请在 Bambu Studio 中打开 `BG001_RMCA_candidate_15_side_support_ON.3mf` 检查真实切片；`BG001_RMCA_candidate_15_side.stl` 是同一打印姿态的整体模型，`BG001_RMCA_candidate_15_side_support_ON.gcode` 与 3MF 中的 G-code 完全一致。

配置：Bambu Lab X1 Carbon、0.2 mm 喷嘴、Bambu ABS Vessel、0.10mm Standard Vessel；层高 0.10 mm，支撑 tree(auto) / tree_organic。预计 1507.5 分钟，总耗材 61.74 g，盒内支撑估计 15.15 g。

顶部清理通路代理判定有 0 个受阻区域。这个数值不等同于实物支撑一定被困；0 个也不代表已验证支撑可以拆除。两种姿态均需人工检查，尚未进行实物拆支撑验证。

原始切片内部对象名使用历史排序编号；当前文件名和 manifest 中的 candidate_id 才是本包的候选编号。本次 candidate 15 的切片排序为 rank 1，未修改 3MF 的内部对象名。

print_transform.json 是倒模坐标到打印板坐标的刚性变换；casting_restore_transform.json 为逆变换。侧放件在打印后需转回开口朝上再灌注。本次源模型的侧壁厚度 2.00 mm、底板厚度 2.00 mm；盒体按本轮配置重新生成并重新切片，血管、半径和端口沿用验收模型。导出阶段仅复制已校验结果。没有发送打印任务。
