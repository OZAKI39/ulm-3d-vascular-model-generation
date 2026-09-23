### 图 06：06_resistance_dissipation_validation

![图 06](figures/06_resistance_dissipation_validation.png)

应该看什么：看固定种子的大量速度下，各类耗散是否出现负值。

实际看到什么：1024 个向量的 self、wall、pair 与总耗散均非负，self 缩放矩阵最小特征值为正。

有没有异常：舍入界与特征值原符号一起保存，没有对负特征值取绝对值。

数据：[06_dissipation.json](data/06_dissipation.json)、[06_dissipation.csv](data/06_dissipation.csv)、[05_sparse.json](data/05_sparse.json)、[05_sparse.csv](data/05_sparse.csv)。
