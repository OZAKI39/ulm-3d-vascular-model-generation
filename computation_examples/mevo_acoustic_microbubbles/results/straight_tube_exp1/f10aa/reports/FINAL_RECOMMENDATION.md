# 推荐后端与后续科研路径

在原物理及原 Gate 不变前提下，本轮显著降低64/128泡完整计算成本：分别 6.012×、6.763×。32泡也加速 2.944×。

32/64泡推荐1线程WaterAccelerated或WaterPacked，原Warp声学+壁面。128单任务采用启动前固定4线程WaterAccelerated，完整1ms实测5.293s（相对原单线程8.384×）。多个独立任务仍每进程1线程；不在运行中调整BLAS线程。Packed在128没有额外明确收益。CuPy保留验证路线，不默认启用。四线程运行时切换的native异常未作为科学结果接受，启动前固定线程数的路线单列验证。

后续先做少量64/128泡更长窗口及多个独立初态/幅相的数值重复，优先确认接触事件与持续形变，而非立即扩到千泡。256预检PASS，需独立初态、多接触参考、dt/cache和墙面验证后才增加长轨迹。可增加完整逐泡研究数据吞吐，但不能将其标为真实PZT训练准入；下一轮动作—状态试验须仍保留固定ID、有限预算、失败及未校准标签。

EXPLORATORY · PHYSICS_UNCALIBRATED · HOMOGENEOUS_WATER_CONTROL_PROXY · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
