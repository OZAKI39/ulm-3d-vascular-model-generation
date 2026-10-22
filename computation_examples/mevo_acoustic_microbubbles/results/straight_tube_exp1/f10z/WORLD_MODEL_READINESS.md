# 世界模型数据准入与简单预测

WORLD_MODEL_DATA_NOT_READY；training_ready=false。
 
多个独立初态、动作、明显形变和短窗数值精度已经具备；本轮没有独立未见动作×未见初态的模型测试，也没有充分证明独立 XY+形状控制。两参数局部响应秩 2；离散动作有限差异秩 5。这些都不能自动升级成联合位置—形状模型准入。
 
简单基准训练来源 24 个 2 ms 初态/动作组；独立 validation/test 初态共 4 条（OFF/PAIR 均是已见动作）。看到前 0.5 ms 后预测剩余 1.5 ms；状态保持、前段恒速和训练动作平均响应。没有随机拆时间点，没有把步长/缓存复算当新独立样本。
 
| 轨迹 | split | 基准 | cx/cy/Rg 合成误差 µm | Rg 误差 µm |
| --- | --- | --- | --- | --- |
| medium_test__OFF__2000us__DT1 | test | state_hold | 0 | 0 |
| medium_test__OFF__2000us__DT1 | test | constant_velocity_first_0p5ms | 0 | 0 |
| medium_test__OFF__2000us__DT1 | test | training_action_mean | 0 | 0 |
| medium_test__PAIR__2000us__DT1 | test | state_hold | 0.348007 | 0.334227 |
| medium_test__PAIR__2000us__DT1 | test | constant_velocity_first_0p5ms | 0.0100095 | 0.0100054 |
| medium_test__PAIR__2000us__DT1 | test | training_action_mean | 0.0593179 | 0.0593033 |
| medium_validation__OFF__2000us__DT1 | validation | state_hold | 0 | 0 |
| medium_validation__OFF__2000us__DT1 | validation | constant_velocity_first_0p5ms | 0 | 0 |
| medium_validation__OFF__2000us__DT1 | validation | training_action_mean | 0 | 0 |
| medium_validation__PAIR__2000us__DT1 | validation | state_hold | 0.303281 | 0.287667 |
| medium_validation__PAIR__2000us__DT1 | validation | constant_velocity_first_0p5ms | 0.00876488 | 0.00876025 |
| medium_validation__PAIR__2000us__DT1 | validation | training_action_mean | 0.01275 | 0.0127432 |
 
这些有界基准不是已经训练的神经世界模型，也没有外推到秒级。联合模型未准入，因此 PyTorch 世界模型 NOT_RUN，MPPI_CONTROL_NOT_RUN；未生成声称模型/MPPI 成功的条件性报告。只用短时 X 平移做后续基线是可行候选，但仍需未见动作评测。
