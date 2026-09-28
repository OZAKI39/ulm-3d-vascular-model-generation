# 代数配置参考及原生索引代码

本轮独立分块实验的选项依据：
- https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/
- https://petsc.org/release/manualpages/PC/PCFieldSplitSetSchurPre/
- https://petsc.org/release/manual/ksp/

文档只用于配置语法，不作为本项目数值通过的证据。实际结果以保存的同场CFD、true residual与源码为准。

原求解器petsc_impl.cpp 747–780自行构建速度/压力IS：`uindx[(dof-1)*i+j]=dof*ltg[i]+j`（j=0,1,2），`pindx[i]=dof*ltg[i]+dof-1`；速度IS blocksize3，随后PCFieldSplitSetIS("0",isu)、PCFieldSplitSetIS("1",isp)。它还设置PCFieldSplitSetBlockSize(pc,2)，因此不能再通过命令行设置strided blocksize4。第一次本轮实验发生这项配置冲突，已删除重复选项，未改求解器源码。
