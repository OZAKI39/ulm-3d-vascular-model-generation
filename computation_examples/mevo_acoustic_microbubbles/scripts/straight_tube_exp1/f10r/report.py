"""Produce Chinese research reports from archived, separately admitted measurements."""
import numpy as np
from .common import ROOT,BASE,CONFIG,read,write,sha,identity


def table(headers,rows):
    cell=lambda value:str(value).replace('|',r'\|')
    return '\n'.join(['| '+' | '.join(map(cell,headers))+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(cell,r))+' |' for r in rows])


def main():
    an=BASE/'analysis';run=read(an/'RUN_ANALYSIS.json');gate=read(an/'FINAL_GATE.json');cmp=read(an/'COMPARISONS.json');cost=read(an/'COSTS.json')
    field=read(BASE/'diagnostics/SINGLE_PZT.json');phase=read(an/'X_PHASE_RESULTS.json');peaks=read(BASE/'diagnostics/AXIS_PAIR_PRESSURE_MAXIMA_REVIEWED.json')
    x=run['X_COMPAT_DT1'];y=run['Y_FOUR_X025_DT1'];pair=run['X_PAIR_DT1'];seg=run['X_SEGMENTED_DT1'];phase10=run['X_PHASE_+10_DT1']
    def pos(r):return ', '.join(f'{v*1e6:.9g}' for v in r['endpoint_m'])
    def delta_rg(r):return (r['rg_initial_final_m'][1]-r['rg_initial_final_m'][0])*1e9
    def error(r,axis='X'):return r['target_error_m'][axis]*1e6
    labels='`EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED`；`training_ready=false`。'
    progress_done=(BASE/'progress/UPDATE_AUDIT.json').exists() and read(BASE/'progress/UPDATE_AUDIT.json')['status']=='PASS'
    common='本报告来自实际保存的2 ms轨迹；单分辨率观察、完整数值准入、实验标定分别统计。'+labels+'\n'
    mainrows=[]
    for name,label,axis in [('OFF_DT1','OFF','X'),('X_COMPAT_DT1','原X全精度重放','X'),('X_PAIR_DT1','同轴X对','X'),('ONE_HOT_SEQUENCE_DT1','单片依次X−/X+','X'),('X_SEGMENTED_DT1','两个1 ms段','X'),('Y_FOUR_X025_DT1','新低压Y','Y')]:
        r=run[name];g=gate['groups'][name]
        mainrows.append([label,pos(r),f'{error(r,axis):.9g}',g['GROUP_POSITION_NUMERICALLY_ADMITTED'],g['GROUP_SHAPE_NUMERICALLY_ADMITTED'],g['PARTICLE_ID_NUMERICALLY_ADMITTED'],'是' if g['TARGET_REACHED'][axis] else ('仅单分辨率/不准入' if g['GROUP_POSITION_NUMERICALLY_ADMITTED']=='NOT_RUN' else '否')])
    result_table=table(['控制','终点COM (µm)','对应目标误差 (µm)','COM数值','形状数值','固定ID','目标准入'],mainrows)
    controls=[]
    for name in ['X_ORIGINAL','X_PAIR','X_SEGMENTED','ONE_HOT_X_SEQUENCE','Y_FOUR_X025']:
        for s in read(BASE/'RUN_SCHEDULES.json')[name]['segments']:
            controls.append([name,f"{s['start_s']*1e3:g}–{(s['start_s']+s['duration_s'])*1e3:g}",', '.join(f'{a:.9g}' for a in s['amplitude_Pa_equivalent']),', '.join(f'{v:.9g}' for v in np.rad2deg(s['phase_rad'])),','.join('ON' if b else 'OFF' for b in s['enabled'])])
    control_table=table(['控制','时间(ms)','Y−,Y+,X−,X+ 幅值(Pa-equivalent)','相位(°)','启停'],controls)
    phases=table(['仅X+变化(°)','团心Fx/原Fx','实际2 ms X(µm)','目标误差(µm)','证据范围'],[
        [row['delta_Xplus_deg'],f"{row['Fx_ratio']:.9f}",
         f"{next(v for v in phase['actual_2ms'] if v['shift_deg']==row['delta_Xplus_deg'])['endpoint_m'][0]*1e6:.9f}" if row['delta_Xplus_deg'] in [-6,6,10] else (f"{x['endpoint_m'][0]*1e6:.9f}" if row['delta_Xplus_deg']==0 else 'NOT_RUN'),
         f"{error(run[f'X_PHASE_{row["delta_Xplus_deg"]:+d}_DT1']):.9f}" if row['delta_Xplus_deg'] in [-6,6,10] else (f'{error(x):.9f}' if row['delta_Xplus_deg']==0 else 'NOT_RUN'),
         '完整dt/cache对照' if row['delta_Xplus_deg'] in [0,10] else ('单分辨率轨迹' if row['delta_Xplus_deg'] in [-6,6] else '仅声场')]
        for row in phase['field']['rows']])
    comparisons=[]
    for name in ['X_COMPAT_DT1','X_PAIR_DT1','X_SEGMENTED_DT1','X_PHASE_+10_DT1','Y_FOUR_X025_DT1']:
        for s in ['DT05','FINE']:
            c=cmp[name+'__'+s]
            comparisons.append([name,s,f"{c['COM_m']*1e9:.6g}",f"{c['fixed_id_m']*1e9:.6g}",f"{c['Rg_m']*1e9:.6g}",f"{c['Q_relative']:.6g}",f"{c['matched_RMS_over_R']:.6g}"])
    comparison_table=table(['控制','对照','COM差(nm)','固定ID差(nm)','Rg差(nm)','相对Q差','去质心匹配RMS/R'],comparisons)
    xdoc=f'''# F10R：X短距离导航、控制模式与有限鲁棒性

{common}
目标始终是绝对位置(0.2,0,0) µm，容差0.075 µm。初始COM仅有浮点舍入量，所有坐标与半径使用米计算。

{result_table}

## 全精度兼容和真实控制序列

原F10Q INPUT控制未经报告舍入，通过新接口重放后全部2001帧、32个ID、三个坐标逐位相同，最大差0 m。OFF逐泡严格静止。接口使用原场/APECSS/F2/Sequential QP，不重求声波传播。

{control_table}

表格仅方便阅读；实际运行参数的全部精度见RUN_SCHEDULES.json。Pa-equivalent不是Vpp或团心声压。单片依次开启是X−400 Pa-equivalent运行1 ms，再X+400运行1 ms，前后运动近乎抵消，终点不达目标；它不是失败的求解器。

两段控制前1 ms只把原X+相位减3°，后1 ms加3°。目标误差由{error(x):.9f}减至{error(seg):.9f} µm，改善约{(error(x)-error(seg))*1e3:.6f} nm，大于本次COM差异，但只是本次预设序列的小改善，不能据一例称分段控制普遍更好。Rg分别变化{delta_rg(x):.6f}/{delta_rg(seg):.6f} nm，均无有意义的指定形状塑造。

## 预登记的相位敏感性

{phases}

中心±1°差分得到3.728917%/度，约为“6%/度”粗算的62.15%。−6/+6/+10°轨迹位移与中心力比例预测相符，三个实际终点均在容差内；因此本例不支持±6°即越界的预测。+10°误差{error(phase10):.9f} µm，余量仅{(.075-error(phase10))*1e3:.6f} nm，不能据此向更大误差或长时间外推。−10°仅算声场，其缩放预测不冒充真实失败轨迹。

## 初态与形态

固定32ID初态绕Z旋转+5/−5°，COM不变且几何合法，没有强制投影。实际终点分别为({pos(run['X_rotate_z_plus5'])})与({pos(run['X_rotate_z_minus5'])}) µm；单分辨率观察均在目标容差内，Rg变化{delta_rg(run['X_rotate_z_plus5']):.6f}/{delta_rg(run['X_rotate_z_minus5']):.6f} nm，始终一个32成员几何图。未做这两条的独立步长/缓存重复，所以COM/形状/ID正式数值准入NOT_RUN，也不是独立实验样本。

{comparison_table}

数值差异是所有共同保存时间的最大值，不只终点。原X与OFF的半步/细场来自原F10Q完整轨迹，精确重放成立后才复用；其他列均为本轮新计算。形状判定同时检查Rg、Q、长短比、去质心匹配分布与可识别主轴，不能由COM通过代替。此初态最大两个特征值近简并，主要伸展方向不稳定可识别时标NOT_IDENTIFIABLE，没有用任意特征向量角度称控制成功。

## 场压力极大值与相位差

等幅同轴对δ=φplus−φminus，从−60°到+60°每10°计算轴线上局部声压极大值，随后跟踪同一连续分支。X支从相位0时x=0开始，拟合位移斜率{peaks['X']['branch_slope_um_per_degree']:.9f} µm/度，绝对值与约0.83 µm/度预期接近，负号由相位定义决定。Y强正侧支从y={peaks['Y']['initial_peak_m']*1e6:.6f} µm开始，斜率{peaks['Y']['branch_slope_um_per_degree']:.9f} µm/度，仅适用于当前近场代理。

初次Y“最近中心峰”选取混用了弱中心极值和外侧峰，错误4.23 µm/度拟合已明确作废；原始峰和错误拟合仍保存，修正版逐分支连续跟踪，不改场数据。压力极大值仅是候选陷阱，尚未证明是三维稳定声力平衡；长时间移动陷阱导航NOT_RUN。
'''
    single=[]
    for row in field['rows']:
        axis=1 if row['channel'].startswith('Y') else 0
        grad=row['gradient_complex_Pa_m'][axis];force=np.array(row['primary_force_N'])*1e12
        single.append([row['channel'],f"{row['pressure_magnitude_Pa']:.9f}",f"{grad['real']:.9g} {grad['imag']:+.9g}i",', '.join(f'{v:.9g}' for v in force)])
    h=read(an/'HISTORICAL_Y_FORCE_SUMMARY.json');s=h['runs']['DT1']
    ydoc=f'''# F10R：Y场效率、历史形态敏感性与新弱控制

{common}
**全部Y结论：HOMOGENEOUS_WATER_PROXY_SPECIFIC，不外推真实PDMS装置。**

## 新Y轨迹前的单圆片场诊断

各自单开1500 Pa-equivalent、2.5 MHz、7 mm圆片，团心在原点，以下是实际代理场及APECSS一次声力；不是实测PZT压力。完整三维复梯度在diagnostics/SINGLE_PZT.json，表中列沿圆片轴的分量，另两分量仅浮点量。

{table(['单片','团心|P|(Pa)','轴向∇P(Pa/m)','一次声力Fx,Fy,Fz(pN)'],single)}

Y轴近场零点为{field['axis_nearfield_zero_m']*1e3:.9f} mm，中心距Y片4.5 mm。声压幅值效率Y/X={field['amplitude_efficiency_Y_over_X']:.10f}，约1/{1/field['amplitude_efficiency_Y_over_X']:.3f}，证实用户粗算；压力幅值比不是功率比。Y−单片一次力朝−Y，Y+朝+Y，不能预设单片总把泡推离PZT。X−朝+X、X+朝−X。解析轴场及同一场有限差分梯度检查通过，详见原始误差记录。没有独立Z导航证据。

## 历史F10Q Y的真实失败仍成立

原控制2 ms只移动约0.051776 µm，目标误差0.148224 µm；COM对照通过，形状与固定ID不收敛。半步/细缓存的最大ID差395.608/764.640 nm、Rg差96.417/306.730 nm、Q差4.04%/12.16%，没有重新包装成成功。

最后实际保存的力采样时刻为1.9 ms，不能称2 ms端点力：原dt的逐泡一次力模平均{s['last_saved_force_pN']['primary_force_N']['mean_norm']:.6f} pN、二次力平均{s['last_saved_force_pN']['secondary_force_N']['mean_norm']:.6f} pN、接触反力平均{s['last_saved_force_pN']['contact_force_N']['mean_norm']:.6f} pN。最后一步活动接触对dt/half/fine为64/64/62，压力约1742–1800 Pa；固定2 ms几何重算F2增量平均{s['endpoint_F2_increment_pN']['mean']:.6f} pN。

相同端点几何、只换缓存时，一次/二次力最大差分别约0.000233/0.000255 pN；不同数值结果的端点几何放回同一缓存，二次力差可达91.323/133.583 pN、F2增量差49.834/76.101 pN。由此可见直接场误差很小，但内部接触重排和几何反馈敏感；步长和缓存改变都存在差异。此冻结状态检查不证明某一个唯一因果，也不证明物理接触模型已校准。

## 新低压力Y：稳定但不足以导航到目标

四个0.5 ms候选分别是单Y+、Y等幅对及原Y控制X对幅值的0.25/0.5倍。单Y+与Y对只移动0.0364/0.0589 nm，按预登记弱信号筛选不做完整窗口；0.25倍X候选形态变化最小，先登记选中再独立做dt/half/fine完整2 ms。0.5倍X只有短窗证据，NOT_RUN完整准入。

新Y仍保持原两Y幅值和相位，X幅值减至约119.334 Pa-equivalent，实际泡处声压{y['actual_pressure_min_max_Pa'][0]:.6f}–{y['actual_pressure_min_max_Pa'][1]:.6f} Pa；这是四路干涉控制而非只开Y片。COM终点({pos(y)}) µm，Y目标误差{error(y,'Y'):.9f} µm，未达目标。基础质量、团COM、形状、拓扑及固定ID都通过完整数值对照；最大ID步长/缓存差仅{cmp['Y_FOUR_X025_DT1__DT05']['fixed_id_m']*1e9:.9f}/{cmp['Y_FOUR_X025_DT1__FINE']['fixed_id_m']*1e9:.9f} nm。

Rg变化{delta_rg(y):.9f} nm，X/Y比从1到{y['xy_initial_final'][-1]:.12f}，三维长短比从{y['aspect_initial_final'][0]:.12f}到{y['aspect_initial_final'][-1]:.12f}；保持初始一个32成员图。微小Rg松弛超过本次差分，但不是指定团形状控制，也不是聚合形成或融合。主要方向近简并不作任意方向成功判定。不能为达到200 nm目标提高Y幅值或放宽形状门槛。

## 下一步

先利用等幅X对相位移动候选压力峰，验证力平衡与可承受相位误差，再有限扩大X距离。Y应先确认实验中真实近场和PDMS影响/电压声压标定；当前Y效率低由水代理特定几何零点主导，不能据它判真实装置失去Y控制能力。若继续代理探索，可在既有预算内研究干涉偏置或适当分段的低压Y，先通过形态与ID检查；任意团级形状控制仍证据不足。
'''
    report=f'''# F10R：仿真四PZT控制序列与XY导航交付

已完成三种仿真控制模式、原X精确重放、真实成对/分段实验、有限相位/初态对照、新Y稳定弱控制及四组图件。X短目标可达；新Y数值稳定但未达目标；没有真实硬件执行。

{common}

## 原始身份与历史保护

最新公开同步HEAD为467b17411deeb13a77370565b7382f071b41a774，本地完整科学源码HEAD为41ae8ed2b88f2cbbed834fffccaa0b5de6e91589，950个源码/配置/测试/依赖文件逐项SHA相同。完整本地F10Q工作区干净，公开包与科研包的不同来自归档范围；最新清单位于sync/20261014，根SYNC_MANIFEST为Stage5旧清单。独立分支dev/fast-pzt-f10r-simulation-control-navigation-20261014，无安全废弃对象可删，清理审计单独提交3e55371324701bb32e9eec229ebd226cff46eb25。所有旧实现、失败及轨迹保留，独立只读核验无旧跟踪文件修改或删除。

{result_table}

## 用户要求的22项回答

1. 是，从最新公开F10Q同步分支的准确HEAD建立独立worktree开发。
2. 已核实公开与本地完整科研版本的950个科学源码/配置等文件字节相同，提交不同不代表科学差异。
3. 无可证明安全删除对象，未强行删代码；审计和活跃入口索引已独立提交。
4. F1–F10Q历史证据完整保留；本轮170个远端文件回传并逐项SHA核验，原始失败不覆盖。
5. 已实现独立段号、起止、频率、启停、幅相、场身份及实际事件日志，非网格事件缩短步长准确对齐。
6. ONE_HOT_PZT、AXIS_PAIR_CONTROL、FOUR_CHANNEL_INDEPENDENT均已通过接口测试并实际用于轨迹。
7. 仍为2.5 MHz均匀水圆片代理与原合法APECSS，缓存复用，无频率外推或静默压力裁剪。
8. 严格区分Pa-equivalent与Vpp；原1500源上限不代表功放安全电压，未建立虚构换算。
9. 没有开发/运行真实仪器通信，也没有USB或实际电压输出。
10. 原X控制全精度重放，全部2001帧、32ID坐标逐位一致，最大差0 m。
11. OFF静止；原X、X对、两段X都达到原目标且通过完整分层检查；单片反向序列基本相消，未达目标，数据见上表。
12. ±5°团初始刚性旋转后观察到X约0.200894 µm，均在容差；仅两种确定性扰动、未做各自dt/cache重复，不泛化鲁棒性。
13. 找到低压力四路Y候选，COM/形状/拓扑/固定ID都通过实际半步与细缓存对照。
14. Y仅移动{y['endpoint_m'][1]*1e6:.9f} µm，目标误差{error(y,'Y'):.9f} µm，未达到0.2 µm目标。
15. 存在：历史F10Q Y的COM通过，但形状/ID失败。新弱Y通过并不抹去历史失败。
16. 本例两段控制比原固定控制目标误差少约{(error(x)-error(seg))*1e3:.6f} nm，形态代价相近；不声称分段普遍占优。
17. IDEAL_SCHEDULED_SEGMENTS；相位公共时钟和理想同步，USB延迟/同步偏差未测，不虚设仪器参数。
18. 最长实际完成2 ms；没有5/10/30 ms、长时间陷阱或任意Z导航。
19. 全部科学数值、分析、测试及独立核验{cost['total_scientific_s']:.6f} s；其中远端四进程{cost['remote_campaign_wall_s']:.6f} s、实际积分{cost['integration_s']:.6f} s（已包含在远端时间中，不重复相加）。CPU QP {cost['CPU_solver_s']:.6f} s、阻力装配{cost['CPU_resistance_s']:.6f} s均为积分内分项；峰值RSS {cost['peak_RSS_MiB']:.3f} MiB，纯GPU核时未独立剖析。
20. 四组PNG/PDF在本工作区results/straight_tube_exp1/f10r/figures/；完整绝对路径见文末。
21. {'已完成备份及最小插入更新，原文全部字节及标题结构保留' if progress_done else '等待工程和图件审查通过后运行update_progress'}；证据progress/UPDATE_AUDIT.json，备份名microbubble_aggregation.md.bak_F10R_20261014。
22. 优先验证相位移动X陷阱的三维力稳定性并有限扩大X距离；Y需避免把代理近零点结论外推实验，先标定/确认场再优化低风险Y。任意团形状控制证据不足，不优先扩大复杂形状任务。

## 补充要求的实际结论

单片场先于新Y轨迹完成：Y/X压力效率约1/{1/field['amplitude_efficiency_Y_over_X']:.3f}，Y零点距4.504167 mm、实际距4.5 mm，确认代理特有的Y低效。每片的复梯度及一次力完整保存在SINGLE_PZT.json和Y报告，全部Y标签HOMOGENEOUS_WATER_PROXY_SPECIFIC。

只改X+相位±1/±3/±6/±10°，局部一次声力敏感性3.73%/度，与6%/度预测不一致。−6/+6/+10°真正运行完整2 ms，均观察到容差内；+10°额外半步/细缓存通过，距边界仅约{(.075-error(phase10))*1e3:.3f} nm。−10°未运行轨迹，不能把场推算称实际越界。

等幅同轴对连续压力峰的位置曲线X约{peaks['X']['branch_slope_um_per_degree']:.6f} µm/度，Y强正侧支约{peaks['Y']['branch_slope_um_per_degree']:.6f} µm/度；都是场计算，未模拟长时陷阱导航，也不是已证明的三维稳定陷阱。

## 失败、未运行与工程修订

25条本轮轨迹均完成各自登记窗口并通过基础质量，无积分拒绝；但这不代表25条导航成功。单片相消序列和新Y未达目标；弱单Y/Y对未进入完整窗口；原Y形态/ID失败保留；旧强控制0.884 ms活动集循环是历史失败，本轮不重启。

原Y压力峰最近中心选取混淆分支，原错误结果保留并明确作废，改为连续强分支审查。另发现输出循环变量覆盖RESULT.source文本（在积分完成后），输入、控制及轨迹均正确；派生报告改读不可变INPUT.source，并恢复原OFF零响应准入规则。原RESULT和远端SHA不改，SOURCE_METADATA_REVIEW.json保留修订前后。没有修改科学轨迹或阈值掩盖失败。

未运行项包括旋转/±6°的完整数值准入、−10/±1/±3°完整轨迹、长时移动陷阱、任意Z或复杂团形状、实验标定与实际硬件。不存在把NOT_RUN当PASS的情况。

## 冻结执行与图件

四个独立Supervisor标签mevo_f10r_pilot_v1、screens_v1、navigation_v1、robustness_v1，均已退出；各自固定源码commit与tar SHA在execution/和logs/*_EXECUTION.json，运行没有重用可变工作目录。四场缓存/APECSS输入及预登记SHA一并保存。

工作区：`{ROOT}`。图件目录：`{BASE/'figures'}`。Fig_F10R_1_segmented_control、Fig_F10R_2_x_navigation_robustness、Fig_F10R_3_y_morphology_sensitivity、Fig_F10R_4_control_admission各有PNG/PDF。输入SHA、脚本身份和图件SHA在FIGURE_PROVENANCE.json，实际视觉审查在VISUAL_REVIEW.json。原始轨迹、所有逐步接触/真实互补性/事件、力采样与分析数组均可追溯。

图形渲染、网络传输与工程开发时间不计入科学600 s预算；没有把这些耗时称为GPU求解。所有物理仍未实验校准，数值可重复不等于实验准确或AI训练准入。
'''
    docs=ROOT/'docs'
    for name,text in [('FAST_PZT_F10R_DEVELOPMENT.md',report),('F10R_X_NAVIGATION_ROBUSTNESS.md',xdoc),('F10R_Y_SHAPE_STABILITY.md',ydoc)]:
        (docs/name).write_text(text,encoding='utf-8')
    inputs=[an/'RUN_ANALYSIS.json',an/'COMPARISONS.json',an/'FINAL_GATE.json',an/'COSTS.json',an/'X_PHASE_RESULTS.json',an/'HISTORICAL_Y_FORCE_SUMMARY.json',BASE/'diagnostics/SINGLE_PZT.json',BASE/'diagnostics/AXIS_PAIR_PRESSURE_MAXIMA_REVIEWED.json',BASE/'RUN_SCHEDULES.json']
    write(BASE/'REPORT_INPUTS.json',dict(identity=identity(),input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
          reports_sha256={p.name:sha(p) for p in docs.glob('*F10R*.md')}))
    print('Wrote development, X robustness and Y sensitivity reports')


if __name__=='__main__':main()
