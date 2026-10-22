# Candidate 0

正放：底板贴打印板，顶部开口朝上。

请在 Bambu Studio 中打开 `BG001_RMCA_candidate_00_upright_support_ON.3mf` 检查真实切片；`BG001_RMCA_candidate_00_upright.stl` 是同一打印姿态的整体模型，`BG001_RMCA_candidate_00_upright_support_ON.gcode` 与 3MF 中的 G-code 完全一致。

配置：Bambu Lab P1S、0.4 mm 喷嘴、Bambu ABS、0.20 mm Standard；树状支撑已开启。预计 110.7 分钟，总耗材 55.77 g，盒内支撑估计 2.22 g。

顶部清理通路代理判定有 28 个受阻区域。这个数值不等同于实物支撑一定被困；0 个也不代表已验证支撑可以拆除。两种姿态均需人工检查，尚未进行实物拆支撑验证。

原始切片内部对象名使用历史排序编号；当前文件名和 manifest 中的 candidate_id 才是本包的候选编号。Candidate 15 对应历史 rank 4，未修改原始 3MF 的内部对象名。

print_transform.json 是倒模坐标到打印板坐标的刚性变换；casting_restore_transform.json 为逆变换。侧放件在打印后需转回开口朝上再灌注。血管、半径、端口、盒体和原始文件均未修改，也没有发送打印任务。
