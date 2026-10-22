# POT匹配与数值复核

EXPLORATORY · PHYSICS_UNCALIBRATED · training_ready=false

POT 0.9.7官方Partial/Unbalanced熵正则求解器；空间平方距离、面积比、峰值差异组成代价，最大位移300 px/s×dt+2 px，超限代价100；reg=0.05、reg_m=1、partial mass=0.85，最多500步。开发选择distance_scale=8。行质量、匹配峰≥0.5、归一化熵≤0.6和真实有限性检查后交给原F10U拟合。

最终138例×2变体复核中，明确numerical errors共118例（Partial 96、Unbalanced 22），全部拒判。18个接受的Partial例存在中间divide警告，但最终矩阵有限、质量约0.85、求解日志误差接近机器精度；保留警告，不把所有warning混为有效收敛证据。数值错误拦截修正前后接受集合变化0例。冻结测试结果未被重新挑选。

Unbalanced最终87/138接受、0个已知错误；CPD风险门槛后58/138、1错误，PIV98/138、3错误。真实整片Unbalanced仅66个新鲜有效状态，不能把合成成功外推到真实录像。端点数量变化/孤立假点能使低正则核出现下溢，本轮没有在最终集上改用更容易通过的超参数。应在新开发数据上研究稳定/log域求解器及显式unmatched项，再用新测试集评价。

软对应不是长期材料ID；对称结构、漏检和固定纹理仍有歧义。tests/test_contracts.py验证已知平移和异常点下安全拒判。
