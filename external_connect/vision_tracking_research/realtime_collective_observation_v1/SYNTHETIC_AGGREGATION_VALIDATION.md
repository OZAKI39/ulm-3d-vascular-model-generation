# SYNTHETIC_AGGREGATION_TEST

开发种子94001用于工程调试。最终冻结种子95001/95002在运行前登记，阈值未按最终结果调整。10类场景×2种子；先用已知点集/速度验证数学趋势，再用真实录像图块、冻结OpenCV检测器、原CuPy PIV估计图像趋势。没有PZT声学模拟、没有训练大模型。

数学方向/拒判正确 20/20；图像方向/拒判正确 19/20；图像误报 0；图像拒判 5。这里将预期不可观测的漏检/失焦拒判算作正确防护，不能将这个比例称为真实聚集精度或独立跨视频泛化率。素材仍来自同一条真实录像，数学种子独立，光学域不独立。

| case / seed | math | measured image | first image convergence s |
|---|---|---|---|
| converge / 95001 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 0.5 |
| translation / 95001 | NO_CONSISTENT_CONVERGENCE | NO_CONSISTENT_CONVERGENCE | 0.6333333333333333 |
| converge_translate / 95001 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 0.5 |
| rotation / 95001 | NO_CONSISTENT_CONVERGENCE | NO_CONSISTENT_CONVERGENCE | None |
| appearance_loss / 95001 | NOT_OBSERVABLE | NOT_OBSERVABLE | None |
| local_converge / 95001 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 0.5 |
| defocus / 95001 | NOT_OBSERVABLE | NOT_OBSERVABLE | None |
| form_cluster / 95001 | VISIBLE_CONVERGENCE | NOT_OBSERVABLE | 0.5 |
| split / 95001 | VISIBLE_DISPERSION | VISIBLE_DISPERSION | None |
| onset / 95001 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 1.4333333333333333 |
| converge / 95002 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 0.5 |
| translation / 95002 | NO_CONSISTENT_CONVERGENCE | NO_CONSISTENT_CONVERGENCE | None |
| converge_translate / 95002 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 0.5 |
| rotation / 95002 | NO_CONSISTENT_CONVERGENCE | NO_CONSISTENT_CONVERGENCE | None |
| appearance_loss / 95002 | NOT_OBSERVABLE | NOT_OBSERVABLE | None |
| local_converge / 95002 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 0.5 |
| defocus / 95002 | NOT_OBSERVABLE | NOT_OBSERVABLE | None |
| form_cluster / 95002 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 0.5 |
| split / 95002 | VISIBLE_DISPERSION | VISIBLE_DISPERSION | None |
| onset / 95002 | VISIBLE_CONVERGENCE | VISIBLE_CONVERGENCE | 1.4333333333333333 |


完整每帧记录、输入图块SHA、clips与失败均保留。Onset出现后需时间窗积累；图像处理和数学窗口的首次提示时间分别记录。真实录像聚集准确率 NOT_VERIFIED。
