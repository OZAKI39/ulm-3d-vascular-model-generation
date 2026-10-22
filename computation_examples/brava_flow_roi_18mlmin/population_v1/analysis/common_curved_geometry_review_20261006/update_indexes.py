"""Append this review's entry without overwriting pre-existing local/server text."""
from common import BASE

ENTRY='''2026-10-06 common curved geometry 更新：[修复与验证总报告](population_v1/analysis/common_curved_geometry_review_20261006/COMMON_CURVED_GEOMETRY_MAPPING_REVIEW_ZH.md)及[三张300dpi图](population_v1/analysis/common_curved_geometry_review_20261006/OPEN_RESULTS.html)。独立、编译门控的 native mixed TET10 体装配补丁通过684项检查；四个P1/P1局部回归的R/K逐位不变。68项边界工程检查通过，但weak Dirichlet仍用四角点几何反解pressure：两个curved wall实际抛出原生异常，curved cap参考坐标不一致。按Phase F停止，最终 **COMMON_GEOMETRY_FIX_NOT_EVALUABLE**，不建议下一轮CFD；759项工程检查与10项边界语义/反解gate失败分别保存。初版Hessian失败及修正记录保留，容差未放宽。生产源码/binary、历史P1场和1565结果不变；新增CFD步、linear solve、trajectory、population均0。任务全部退出、无自动后续。'''

def main():
    for name in ['README.md','ACTIVE_DEVELOPMENT_SCOPE.md','CURRENT_RESULTS_ZH.md']:
        p=BASE/name;text=p.read_text()
        if ENTRY not in text:p.write_text(text.rstrip()+'\n\n'+ENTRY+'\n')

if __name__=='__main__':main()
