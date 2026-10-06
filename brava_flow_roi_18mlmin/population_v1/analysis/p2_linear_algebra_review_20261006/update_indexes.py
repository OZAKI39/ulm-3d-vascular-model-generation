from common import *
paragraph='2026-10-06 最新线性代数审计：[主报告](population_v1/analysis/p2_linear_algebra_review_20261006/P2_LINEAR_ALGEBRA_AND_RESOURCE_REVIEW_ZH.md)与[三张300dpi核查图](population_v1/analysis/p2_linear_algebra_review_20261006/OPEN_RESULTS.html)。223,119个压力边占位已定位；原生P2实测111.79M结构项，ASM2子域NNZ重复1.85748倍。ASM0和Schur/GAMG均通过P1单步解等价回归；选定Schur在P2真实装配后、selfp CUDA临时工作区分配失败，尚无线性迭代和完整PC峰值。另独立复现原生边压力回填的曲面反解异常。最终 **P2_PIPE_LINEAR_ALGEBRA_FEASIBLE=NO / STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT**；一次非预期CUDA分配失败已保留并补充拒绝未知显存预算的入口保护，无重试，host OOM=0。146项工程检查通过；完整P2 CFD/BraVa/trajectory/population均0。所有自身任务退出，不自动继续、不push。下方旧轮次仅为历史记录。'
for name in ['README.md','ACTIVE_DEVELOPMENT_SCOPE.md','CURRENT_RESULTS_ZH.md']:
 p=BASE/name;s=p.read_text()
 if '最新线性代数审计：' in s:continue
 lines=s.splitlines(True);idx=next(i for i,l in enumerate(lines) if l.startswith('# '))+1
 lines.insert(idx,'\n'+paragraph+'\n');p.write_text(''.join(lines))
print('Updated three current entry indexes; prior content retained')
