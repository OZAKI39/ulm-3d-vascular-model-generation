# 原作者数值保护补充（正式 guarded campaign 之前）

首次 probe 在迭代 16 完成后，下一次局部 Rayleigh 商为 -0.0114716792 而停止。
原始状态诊断、配置及首次代码保存于 initial_small_probe，未删除失败结果。
进一步核查 Madenci & Oterkus (2014) 教材官方附录，发现作者实际程序有正文未明说的保护：
商非正或位移分母零时 c=0；dtau=1 时 c>2 则改为1.9。
来源：[Springer 出版社配套材料](https://extras.springer.com/?query=978-1-4614-8464-6)，
extras_peridynamics.zip / extras/block_under_compression.f90，497–509 行；
代码和下载归档 SHA 见 resources/SPRINGER_CODE_IDENTITY.json。
只重新表达算法，没有复制或 vendor 其 Fortran 程序；不主张该附件具有开放源代码许可证。

新增 config streaming_erosion_equilibrium_guarded_v1.json 指定此已核实的保护。
不修改旧 config/plan；两次方案属于不同的明确版本。所有收敛容差、载荷和材料不变。
这不是绝对值处理，也不宣称是 2010 论文正文明确写出的分支。记录每一步的原始商及保护计数。
2010 eq19 使用总力；作者例子外力恒定故只差分内力。当前牵引依赖 x，因此继续使用总力差分。
为了保持来源明确，当前带作者保护版本只支持 dtau=1。
