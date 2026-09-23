### 图 05：05_sparse_resistance_matrix

![图 05](figures/05_sparse_resistance_matrix.png)

应该看什么：看 12 个球的非零块分布，以及候选筛选是否漏掉近场项。

实际看到什么：复用 P4 候选查询后，稀疏矩阵和求解结果与穷举完全一致。

有没有异常：查询包围盒扩展只服务本轮近场验证，不改变球尺寸，也不是生产邻居层。

数据：[05_sparse.json](data/05_sparse.json)、[05_sparse.csv](data/05_sparse.csv)。
