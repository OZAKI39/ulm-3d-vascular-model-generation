from common import *
entry='2026-10-06 重新执行结果：[中文主报告](population_v1/analysis/p2_linear_algebra_rerun_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)与[三张300dpi图](population_v1/analysis/p2_linear_algebra_rerun_20261006/OPEN_RESULTS.html)。新候选保留GPU主系统，将Schur/GAMG预条件器放到CPU，P1同检查点回归通过（124.4秒、Δu约1.86e−16 m/s）。一次P2真实系统setup因任务RSS达到25.22 GiB触发预声明24 GiB保护停止；无主机OOM、无CUDA分配失败。瞬时保守余量9.90 GiB仍高于所需8.41 GiB，但内存持续增长且setup未完成，不能称安全门通过。最终 **NO / STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT**，只限定本轮资源保护范围；完整P2 CFD、BraVa、轨迹和群体均0。曲面压力回填gate仍FAIL。所有本轮任务退出，不push；下方上一轮结果保留为历史。'
for name in ['README.md','CURRENT_RESULTS_ZH.md','ACTIVE_DEVELOPMENT_SCOPE.md']:
 p=BASE/name;s=p.read_text()
 if entry not in s:
  first=s.index('\n')+1;p.write_text(s[:first]+'\n'+entry+'\n'+s[first:])
