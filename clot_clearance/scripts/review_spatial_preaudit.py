"""Close the geometry hard gate with reviewable evidence, not a refinement claim."""
from pathlib import Path
import json
import os
import hashlib
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
from scripts.audit_spatial_geometry import ROOT, OUT


def read(p): return json.loads(p.read_text())
def write(p,d): p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')


def cases(path):
    result = {}
    for node in ET.parse(path).iter('testcase'):
        fail = node.find('failure'); skipped = node.find('skipped')
        result[node.attrib['classname']+'::'+node.attrib['name']] = dict(
            status='failed' if fail is not None else 'skipped' if skipped is not None else 'passed',
            failure_message=None if fail is None else fail.attrib.get('message'),
            traceback=None if fail is None else fail.text)
    return result


def figures():
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
    plt.rcParams.update({'font.family':'Arial','font.size':12,'font.weight':'normal',
        'figure.facecolor':'black','axes.facecolor':'black','savefig.facecolor':'black',
        'text.color':'white','axes.labelcolor':'white','axes.edgecolor':'#999999',
        'xtick.color':'white','ytick.color':'white','pdf.fonttype':42})
    fig,axs=plt.subplots(1,2,figsize=(12,6));fig.subplots_adjust(bottom=.22,wspace=.32)
    for ax,extent,color,title in [(axs[0],.1875,'#49b8ff','Full voxel support: 0.375 mm'),
                                  (axs[1],.125,'#ffb35b','Point hull proposal: 0.25 mm')]:
        ax.add_patch(Rectangle((-extent,-extent),2*extent,2*extent,fill=False,ec=color,lw=2))
        coords=np.array([-.125,0,.125]);x,z=np.meshgrid(coords,coords)
        ax.scatter(x,z,s=65,c=np.where(z.ravel()==0,'#b8bec9','#75dda5'))
        for plane in [-.125,.125]: ax.axhline(plane,c='#75dda5',ls='--',lw=1)
        ax.set(xlim=(-.21,.21),ylim=(-.21,.21),xlabel='x (mm)',ylabel='z (mm)',title=title)
        ax.set_aspect('equal');ax.grid(alpha=.2)
    fig.suptitle('Same point coordinates; different physical volume interpretations',fontsize=15)
    fig.text(.5,.075,'Green points / dashed lines: prescribed point planes at z = +/-0.125 mm\n'
             'No refined lattice generated. The 0.25 mm proposal has NOT been adopted.',ha='center',fontsize=11)
    for ext in ['png','pdf']:fig.savefig(OUT/'figures'/f'geometry_domain_ambiguity.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    rows=read(OUT/'comparisons/ISOLATED_COARSE_QUADRATURE_COMPARISON.json')
    fig,axs=plt.subplots(1,3,figsize=(13,5.5));fig.subplots_adjust(bottom=.28,wspace=.42)
    for ax,key,factor,title in zip(axs,['volume_m3','elastic_energy_J','nominal_reaction_N'],
        [1e9,1e12,1e6],['Represented volume (mm3)','Affine elastic energy (pJ)','Top reaction (uN)']):
        vals=[r[key]*factor for r in rows]
        ax.bar([0,1],vals,color=['#49b8ff','#ffb35b'])
        ax.set_xticks([0,1],['Legacy voxel','Hull proposal']);ax.set_title(title,fontsize=12)
        ax.set_ylim(0,max(vals)*1.24)
        for i,value in enumerate(vals):ax.text(i,value+max(vals)*.035,f'{value:.7g}',ha='center',fontsize=11)
        ax.grid(axis='y',alpha=.2)
    fig.suptitle('Isolated quadrature comparison at the same 1% affine strain',fontsize=15)
    fig.text(.5,.055,'Same 27 points, 185 bonds, material and physical horizon\n'
             'Proposal only: changed volume weights; not a new branch or a mesh convergence result',ha='center',fontsize=11)
    for ext in ['png','pdf']:fig.savefig(OUT/'figures'/f'coarse_quadrature_comparison.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)


def main():
    old=cases(ROOT/'verification/pfpd_at1_branch_stability_tracking_v1/tests/FINAL_FULL.xml')
    current=cases(OUT/'tests/BASELINE.xml'); new=cases(OUT/'tests/PREAUDIT.xml')
    changed=[k for k,v in old.items() if current.get(k)!=v]
    test=dict(old_passed=sum(r['status']=='passed' for r in current.values()),
        old_failed=sum(r['status']=='failed' for r in current.values()),
        old_skipped=sum(r['status']=='skipped' for r in current.values()),
        known_failure_IDs=[k for k,r in current.items() if r['status']=='failed'],
        previous_test_record_differences=changed,
        exact_original_failure_messages_and_tracebacks_unchanged=not changed,
        new_tests_passed=sum(r['status']=='passed' for r in new.values()),
        new_tests_failed=sum(r['status']=='failed' for r in new.values()),
        run_policy='Original complete suite rerun once; new isolated preaudit tests run separately. No redundant full-suite rerun.')
    write(OUT/'tests/COMPARISON.json',test)
    manifest=read(OUT/'provenance/BEFORE_MANIFEST.json');modified=[];missing=[];links=[]
    for name,record in manifest['files'].items():
        p=ROOT/name
        if not p.is_file() or p.is_symlink():missing.append(name);continue
        with open(p,'rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
        if p.stat().st_size!=record['bytes'] or sha!=record['sha256']:modified.append(name)
    for name,target in manifest['links'].items():
        p=ROOT/name
        if not p.is_symlink() or os.readlink(p)!=target:links.append(name)
    protection=dict(status='PASSED' if not(modified or missing or links) else 'FAILED',
        original_files=len(manifest['files']),original_links=len(manifest['links']),
        production_files_modified=sum(n.startswith('pd_clot/') for n in modified),
        historical_files_modified=sum(n.startswith(('results/','verification/')) for n in modified),
        existing_files_modified=modified,missing=missing,changed_links=links)
    write(OUT/'provenance/PRESERVATION_CHECK.json',protection)
    repro=read(OUT/'level_0/COARSE_REFERENCE_REPRODUCTION.json')
    facts=read(OUT/'preaudit/GEOMETRY_FACTS.json')
    status=dict(SPATIAL_DISCRETIZATION_AUDIT='INCOMPLETE',STOP_REASON='NEED_GEOMETRY_DECISION',
        PHYSICAL_DOMAIN_CONSISTENCY='NOT_VERIFIED',PHYSICAL_HORIZON='FIXED',
        PHYSICAL_HORIZON_SCOPE='Legacy reference and isolated same-point volume comparison only; no refined lattices exist',
        REFINEMENT_LEVELS_COMPLETED=['LEVEL_0_LEGACY_REFERENCE'],
        INTACT_MECHANICS_ACROSS_LEVELS='NOT_VERIFIED',STABILITY_INTERVAL_CONSISTENCY='INCONCLUSIVE',
        REACTION_AT_ONSET_CONSISTENCY='INCONCLUSIVE',CRITICAL_CONE_INSTABILITY='INCONCLUSIVE',
        ACTIVE_SET_BEHAVIOR='DIFFUSE_ALL_BONDS',ACTIVE_SET_SCOPE='Legacy level 0 only; cross-resolution behavior not assessed',
        CRITICAL_MODE_LOCALIZATION='INCONCLUSIVE',BOUNDARY_ASSOCIATION='INCONCLUSIVE',
        CURRENT_3x3x3_INSTABILITY='INCONCLUSIVE',MESH_OBJECTIVITY='NOT YET ESTABLISHED',
        PFPD_AT1_EVOLUTION_V1='REMAINS_NOT_VERIFIED',PRODUCTION_MODEL='UNCHANGED',
        LEGACY_REFERENCE_REPRODUCTION=repro['status'])
    write(OUT/'STATUS.json',status)
    figures()
    report=f'''# PFPD_AT1_SPATIAL_DISCRETIZATION_AUDIT_V1

**INCOMPLETE — NEED_GEOMETRY_DECISION。** 已完成几何预审、同一粗点云的隔离积分比较、旧临界区间复现和基线测试。根据用户提示词第3节hard gate，尚未生成细网格。该停止发生在物理问题定义阶段，不是数值求解失败。

## 1. 为什么做这一轮
检查原约1%应变失稳是否在固定物理问题、仅细化粒子离散后仍然存在。必须先明确“固定的物理问题”是什么。

## 2. h和delta
h是点间距；delta是非局部物理视界。旧h=0.125 mm，delta=0.25125 mm。后续主研究必须固定delta，不能按h缩放。当前尚无跨网格结果。

## 3. 是否保持相同物理问题
尚未确定。旧代码的体素材料域边长0.375 mm；点包围盒和夹持中心距离0.25 mm。直接固定其中一个并复用全部旧规则，会改变另一个。

## 4. Geometry / volume audit
详见[完整预审](preaudit/GEOMETRY_AND_QUADRATURE_AUDIT.md)。旧体积总和0.052734375 mm³，等于0.375³；点包围盒体积0.015625 mm³，二者相差3.375倍。顶面面积也存在0.140625与0.0625 mm²两种解释。代码约束的是z=±0.125 mm处的点，没有声明细化时将其延拓为精确边界面还是完整刚性体素层。

同27点的守恒包围盒积分规则已另做比较，**未采用**。固定仿射位形下，两种规则的能量密度均为约0.1010025 Pa，但总能量和反力改变。因此它不能被隐藏为普通细化。

## 5. 三档mesh信息
|档位|点数|h mm|delta/h|键数|平均family|
|---|---:|---:|---:|---:|---:|
|旧Level 0|27|0.125|2.01|185|13.7037037|
|Level 1|未生成|未确定|未确定|N/A|N/A|
|Level 2|未生成|未确定|未确定|N/A|N/A|
旧27点全部shape rank=3，全部horizon受边界截断，没有完整球形内部family。

## 6. Intact mechanics
隔离粗积分比较的Fbar/Ftilde仿射误差均小于2e-15；独立SVK能量密度检查通过。该检查没有替代每档完整mechanics gate。Julia/upstream核对尚未执行，因为几何hard gate早于细化及外部reference阶段。

## 7. Stability interval
旧Level 0在新目录按本轮11点加载表重新求驻点与曲率，复现区间{repro['reproduced_interval']}，宽{repro['width']:.12g}。185 strong→185 weak仍成立，位形、相场、反力、能量、残差和负曲率方向均通过与旧状态的比较。见[复现记录](level_0/COARSE_REFERENCE_REPRODUCTION.json)。

## 8. 临界区间是否收敛
INCONCLUSIVE。Level 1/2未生成，不能给出空间一致性PASS。旧解析阈值为0.01，尚无跨网格threshold检查。

## 9. Reaction / stress是否收敛
INCONCLUSIVE。旧1%应变反力约4.2822421875 μN；按两种候选面积归一化分别是30.4515 Pa和68.515875 Pa。材料域尚未选定，因此不把其中一个命名为本轮权威nominal stress。提议包围盒守恒积分同一仿射位形反力约1.2688125 μN，对0.0625 mm²为20.301 Pa；这不是新的分支稳定性结果。

## 10. Active-set switch
旧Level 0仍为100%键共同切换，其他档位N/A。

## 11. Negative-cone direction
旧参考负曲率方向已复现，原缩放下上端点Rayleigh约−0.01470214。该缩放与DOF数有关，没有冒充跨网格统一物理尺度。精确比较见复现JSON。

## 12. Mode局域化
跨网格PR、top-10%幅值占比、公共网格相关性尚未计算，不将缺失值写成0。

## 13. Boundary association
尚未进行跨网格判定；旧算例全部family截断不能单独证明boundary artifact。

## 14. 3×3×3分类
**INCONCLUSIVE**。旧失稳仍可复现，但其空间稳健性不能由一个分辨率或积分规则替换来证明。

## 15. Tests与保护
本轮实际重跑旧测试：{test['old_passed']} passed / {test['old_failed']} known failed / {test['old_skipped']} skipped；新增预审测试：{test['new_tests_passed']} passed。旧失败ID、消息及traceback与上一轮逐项比较：{'一致' if not changed else '存在差异，见COMPARISON.json'}。新增测试不修改旧失败。
保护核验：{protection['original_files']}个旧文件、{protection['original_links']}个链接；状态{protection['status']}。生产修改{protection['production_files_modified']}，历史修改{protection['historical_files_modified']}。

## 16. 当前状态
完整标签见[STATUS.json](STATUS.json)。主研究INCOMPLETE、物理域一致性NOT_VERIFIED、区间／反力一致性INCONCLUSIVE；模型仍REMAINS_NOT_VERIFIED。完成两个预审PNG/PDF图对：[几何定义](figures/geometry_domain_ambiguity.png)、[同点云积分比较](figures/coarse_quadrature_comparison.png)。提示词要求的九组跨分辨率图待几何决定和实际计算后生成，未用虚构网格或结果占位。

## 17. 唯一下一步
明确固定的材料域和夹持定义。建议选择0.25 mm边界节点立方体、精确上下夹持面、守恒边界体积；旧粗结果独立保留，修正coarse另名并重新通过全部门槛。另一选择是0.375 mm体素域与原0.125 mm厚夹持区，但需据此另行设计点云。等待用户对几何含义的回答，不能把沉默当作选择。
'''
    (OUT/'REVIEW_SPATIAL_DISCRETIZATION_AUDIT.md').write_text(report)
    print(json.dumps(dict(tests=test,preservation=protection,reproduction=repro['status'],status=status),indent=2))
    if changed or protection['status']!='PASSED' or test['new_tests_failed'] or repro['status']!='VERIFIED':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
