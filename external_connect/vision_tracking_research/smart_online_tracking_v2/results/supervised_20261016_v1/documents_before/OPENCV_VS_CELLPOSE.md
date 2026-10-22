# 两条识别路线的实际结果

OpenCV 复用原局部背景 / 对比度及连通域实现，Norfair 负责关联。额外保留较大连通域作为未审核区域；面积分成“小 / 大”只是候选分组，不是学到的单泡 / 团语义。原阈值保留，本轮没有在无真值测试集上调参。

Cellpose 3.1.1.3 的 cyto3 是本轮唯一 AI 检测模型，作为成本较低的 Cellpose 系替代。灰度 channels=[0,0]，diameter=30（内部尺度比为 1），无增强，flow_threshold=0.4、cellprob_threshold=0、min_size=15。预测 score 是 cell probability 的映射，不是已校准的微泡置信度。没有训练或微调。

以下为同一录像、相同原始像素、相同帧、相同 RTX4060 Laptop/CPU、同一计时方法的成对检测；每次 CUDA 同步后结束计时。全图为 1476×1120，ROI 为 x650 y350 w180 h350。这里只测 detection，不把全图 detection 时间误称完整在线时间。

| 范围 | 路线 | 帧数 | 平均检测 ms | p95 ms | 平均候选数 |
|---|---|---:|---:|---:|---:|
| FULL_FRAME | opencv | 12 | 33.47 | 55.44 | 624.25 |
| FULL_FRAME | cellpose | 12 | 863.11 | 942.63 | 6.08 |
| ROI | opencv | 12 | 2.16 | 2.39 | 38.58 |
| ROI | cellpose | 12 | 68.63 | 84.60 | 3.17 |

90 帧因果离线回放的完整引擎 p95：OpenCV 7.03 ms，Cellpose 81.16 ms。它是按顺序实际计算的离线对照，不是实时输入计时。

目视实际叠加：OpenCV 选出很多小对比区域；cyto3 候选更少，常偏向较大的模糊结构。这个差异**不能证明**后者背景误检少，也不能证明前者微泡召回高。当前没有完整人工真值，Precision、Recall、F1、定位 / 分割和 ID 误差均 NOT_EVALUATED。直播 Cellpose 在当前 ROI 没有可用候选，并频繁超出 100 ms 软件年龄门槛；失败保留。

当前工程选择：保留 OpenCV + Norfair 作为实时候选管线；Cellpose 暂不承担 60 fps 每帧输出。待人工标签确定实际漏检和误检后，才判断是否进行一次有限微调，不能仅凭“新模型”或候选数量宣布优胜。

许可证及权重身份见 [MODEL_AND_ENVIRONMENT.md](MODEL_AND_ENVIRONMENT.md)。
