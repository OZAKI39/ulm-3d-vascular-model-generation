# 全窗口数值误差与原始门槛

共 45 组完整窗口比较，45 组通过原连续指标准入。独立复核 299000 个接受步；失败前缀不冒称完整通过。逐泡对照检查所有共同物理时刻；额外形状指标每10接受步及终点检查。瞬时接触拓扑另列，不能由轨迹平滑或全局守恒替代。

| 对照 | 逐泡最大差异 / nm | Rg最大差异 / nm | 原始连续指标准入 | 事件拓扑一致 |
| --- | --- | --- | --- | --- |
| refinement__S0_r0__A7_X210__2000us__FINE | 0.003868 | 2.359e-07 | True | True |
| refinement__S0_r0__A7_X210__2000us__DT05 | 0.0002881 | 7.955e-07 | True | True |
| refinement__S2_r0__A7_X210__2000us__FINE | 0.03692 | 0.001891 | True | True |
| refinement__S2_r0__A7_X210__2000us__DT05 | 3.083 | 0.2412 | True | True |
| refinement__S2_r0__A5_X180__2000us__DT05 | 0.08737 | 0.003663 | True | True |
| refinement__S2_r0__A5_X180__2000us__FINE | 0.007198 | 6.575e-05 | True | True |
| refinement__S4_r0__A7_X210__2000us__FINE | 0.05606 | 0.004603 | True | True |
| refinement__S4_r0__A7_X210__2000us__DT05 | 19.13 | 0.6778 | True | True |
| refinement__S4_r0__A5_X180__2000us__DT05 | 0.06196 | 0.001761 | True | True |
| refinement__S4_r0__A5_X180__2000us__FINE | 0.006247 | 0.0008937 | True | True |
| refinement__S6_r0__A7_X210__2000us__FINE | 0.02201 | 0.0007673 | True | True |
| refinement__S6_r0__A7_X210__2000us__DT05 | 2.104 | 0.03391 | True | True |
| refinement__S6_r0__A5_X180__2000us__FINE | 0.0131 | 0.0003756 | True | True |
| refinement__S6_r0__A5_X180__2000us__DT05 | 0.05076 | 0.001002 | True | True |
| confirmation__S2_r1__A5_X180__2000us__FINE | 0.006881 | 0.0005159 | True | True |
| confirmation__S2_r1__A5_X180__2000us__DT05 | 0.06977 | 0.001591 | True | True |
| confirmation__S2_r2__A5_X180__2000us__FINE | 0.006504 | 0.000245 | True | True |
| confirmation__S2_r2__A5_X180__2000us__DT05 | 0.08536 | 0.003789 | True | True |
| sensitivity_checks__S2_r0__J_amp_m_half__2000us__FINE | 0.03708 | 0.001885 | True | True |
| sensitivity_checks__S2_r0__J_amp_m_half__2000us__DT05 | 2.941 | 0.2314 | True | True |
| sensitivity_checks__S2_r0__J_amp_p_half__2000us__FINE | 0.03677 | 0.001889 | True | True |
| sensitivity_checks__S2_r0__J_amp_p_half__2000us__DT05 | 3.234 | 0.2512 | True | True |
| sensitivity_checks__S2_r0__J_phase_m_half__2000us__FINE | 0.03544 | 0.001798 | True | True |
| sensitivity_checks__S2_r0__J_phase_m_half__2000us__DT05 | 2.853 | 0.1404 | True | True |
| sensitivity_checks__S2_r0__J_phase_p_half__2000us__FINE | 0.04976 | 0.001952 | True | True |
| sensitivity_checks__S2_r0__J_phase_p_half__2000us__DT05 | 5.647 | 0.2735 | True | False |
| sensitivity_checks__S4_r0__J_amp_m_half__2000us__FINE | 0.05581 | 0.004591 | True | True |
| sensitivity_checks__S4_r0__J_amp_m_half__2000us__DT05 | 18.28 | 0.6475 | True | True |
| sensitivity_checks__S4_r0__J_amp_p_half__2000us__FINE | 0.05695 | 0.004628 | True | True |
| sensitivity_checks__S4_r0__J_amp_p_half__2000us__DT05 | 20.32 | 0.7228 | True | False |
| sensitivity_checks__S4_r0__J_phase_m_half__2000us__FINE | 0.07507 | 0.003949 | True | True |
| sensitivity_checks__S4_r0__J_phase_m_half__2000us__DT05 | 4.081 | 0.3213 | True | False |
| sensitivity_checks__S4_r0__J_phase_p_half__2000us__FINE | 0.04603 | 0.003762 | True | True |
| sensitivity_checks__S4_r0__J_phase_p_half__2000us__DT05 | 42.35 | 2.074 | True | False |
| extension__S2_r0__A5_X180__3000us__DT05 | 0.136 | 0.004484 | True | True |
| extension__S2_r0__A5_X180__3000us__FINE | 0.01244 | 7.842e-05 | True | True |
| extension__S2_r0__A5_X180__5000us__DT05 | 0.136 | 0.005501 | True | True |
| extension__S2_r0__A5_X180__5000us__FINE | 0.03021 | 0.0006588 | True | True |
| target_confirmation__S2_r0__T_A700__2000us__DT05 | 0.03074 | 0.0005282 | True | True |
| target_confirmation__S2_r0__T_A700__2000us__FINE | 0.004576 | 5.841e-05 | True | True |
| target_confirmation__S2_r1__T_A700__2000us__DT05 | 0.01689 | 0.0002851 | True | True |
| target_confirmation__S2_r1__T_A700__2000us__FINE | 0.002972 | 0.0002061 | True | True |
| target_confirmation__S2_r2__T_A700__2000us__DT05 | 0.03093 | 0.0008049 | True | True |
| target_confirmation__S2_r2__T_A700__2000us__FINE | 0.003305 | 0.0001623 | True | True |
| target_confirmation__S2_r0__T_A700__2000us__DT1 | 0 | 0 | True | True |

形变准入要求 similarity RMS>max(5×数值差异,1nm)，有效F的偏伸展/剪切>max(5×差异,1e-4)。S0形变低于阈值。原50nm逐泡、接触/力学/墙面预算均未放宽。

适用范围：`EXPLORATORY` · `HOMOGENEOUS_WATER_CONTROL_PROXY` · `PHYSICS_UNCALIBRATED` · `ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED` · `training_ready=false`。Y 与四源静态结果另加 `HOMOGENEOUS_WATER_PROXY_SPECIFIC`。
