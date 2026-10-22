# 不确定性与有效性

SNR/peak2peak、局部中值、纹理、空间覆盖、条件数、残差与原六参数valid共同决定可发布性。它们是启发式风险指标，不是已校准概率或95%保证。TorchPIV只有原生ratio flag，未捏造SNR幅值。

Figure7在实际已保存局部矢量上复用F10U四象限block bootstrap，24次重采样，显示nominal95%区间。重叠PIV窗口不是独立样本；目前这些区间coverage没有完成独立校准，因此标为calibrated=false，不用于控制准入。

Xue/Charonko/Vlachos相关峰SNR不确定性研究（https://doi.org/10.1088/0957-0233/25/11/115301）说明该方向值得研究，本轮没有声称实现论文的标定不确定度模型。与其他后端一致或曲线平滑不构成真实准确率证据。

final已知程序运动仍有对称旋转超容差和非仿射错接。真实录像没有速度真值、比例尺和PZT日志；REAL_PHYSICAL_VELOCITY_ACCURACY=NOT_VERIFIED。schema中scientific_accuracy_verified=false。
