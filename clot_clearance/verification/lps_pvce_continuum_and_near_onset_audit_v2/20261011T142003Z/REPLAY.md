# 最小复核与再运行限制

现有服务器：`source /workspace/y_clot/env.sh`，进入`/workspace/y_clot/project`，运行`python -m analysis.lps_pvce_v2.replay --verify-only`。此命令只核验清单和储能差，绝不重新求解。

WSL轻量复核只需Python标准库，无须安装NumPy或加载原生求解库。完整原生重算另需锁定LAMMPS库与energy.so，版本身份在保护清单，不能用其他库悄悄替代。FEM额外安装scikit-fem12.0.2到独立target目录（环境记录中给出），不得升级原环境。

实际数值启动命令与Supervisor配置在provenance/*_launch.sh、*_supervisor.conf；六个FEM ATTEMPT和一个LPS ATTEMPT对应七个独立求解，没有重复计算。旧run目录不可用于覆盖重放。新数值试验需新的路径、冻结方案及授权范围，本轮到停止条件后不自动追加。

三份重复body.data完整文本留服务器；它们是同一个完整未切节点/参考键表，可由旧V1 reference.npz和锁定write_typed_data重建。WSL已保存必要全精度网格、P2/P1自由度、节点状态、力、alive、m及全部数值审计。body.data并非唯一科学数据源。

FILE_MANIFEST.json不收录自己、集合收据及核验回执，以免循环哈希；WSL缺失仅允许上述body.data。最终清单核验同时在服务器和WSL执行。没有迁移其他项目、旧大轨迹、软件环境或视频。
