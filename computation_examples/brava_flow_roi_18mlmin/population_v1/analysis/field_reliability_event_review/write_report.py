"""Evidence-linked Chinese delivery; classifications are per feature, not PASS."""
from common import *
from analyze import angle
import datetime

def main():
    m=read(HERE/'FIELD_METRICS.json');s=read(HERE/'SPATIAL_FEATURES.json');c=read(HERE/'FIELD_COMPARABILITY.json')
    balance={name:read(HERE/f'{name}_BALANCE.json') for name in NAMES};protocol=read(HERE/'EVENT_FIELD_PROTOCOL.json')
    supported='FIELD_FEATURE_SUPPORTED_ON_TESTED_RESOLUTIONS';sensitive='FIELD_RESOLUTION_SENSITIVE';defect='LOCAL_NUMERICAL_DEFECT_CONFIRMED'
    classifications={
      '1189':dict(background_u={'later_high_speed':supported,'early_speed_vector':sensitive},near_wall_direction_and_derivatives=[defect,sensitive],
        mechanism='进入几十 mm/s 背景流的空间特征重复；近壁强朝壁面结构仍存在，但方向占比、剪切和释放环境不稳定。不能据此验证原离壁时刻。'),
      '1329':dict(background_u={'later_high_speed':supported,'exact_local_values':sensitive},near_wall_direction_and_derivatives=[sensitive],
        mechanism='慢后快的切向背景分层在三解中重复；未出现1189那种近乎纯法向主导。法向分量及导数仍有可测差异，不能认证离壁时序。'),
      '642':dict(background_u={'later_high_speed':supported,'exact_local_values':sensitive},near_wall_direction_and_derivatives=[sensitive],
        mechanism='切向流向、整体进入快层趋势重复；后段代表速度原场偏低，近壁导数与控制体残差仍随网格改变。不能合并为1189的强法向缺陷。'),
      '822':dict(background_u={'predominant_streamwise_direction':supported,'original_36_228_mm_s_gain':sensitive,'rise_length':sensitive,'increment_vector_direction':sensitive},
        near_wall_direction_and_derivatives={'wall_release':'NOT_APPLICABLE: interior event; fixed path frame used','G':sensitive},
        mechanism='正的端点背景增量在三解中存在，但原36.228 mm/s幅值及96.672 um的固定10–90%尺度未被重复；中细档增量自身也未稳定。原场内部解释不能升级为可信的急增结构。')}
    for pid,e in classifications.items():
        e['physical_validation']='PHYSICAL_VALIDATION_NOT_PERFORMED'
        e['accepted_position_coverage']={name:c['fields'][name]['valid_event_probes'][pid]['accepted'] for name in NAMES}
    increments={}
    for name in NAMES[1:]:
        a=np.array(m['822']['representatives']['M0']['delta_u_mm_s']);b=np.array(m['822']['representatives'][name]['delta_u_mm_s'])
        increments[name]=dict(angle_to_M0_increment_deg=float(angle(a[None,:],b[None,:],1e-11)[0]),
            scalar_gain_difference_over_M0_gain=(m['822']['representatives'][name]['scalar_speed_increase_mm_s']-m['822']['original_event_gain_mm_s'])/m['822']['original_event_gain_mm_s'])
    result=dict(status='COMPLETED_CONTROLLED_EVENT_FIELD_REVIEW',new_CFD_count=0,reused_CFD_count=3,reused_CFD=NAMES,reference_field_count=2,
      new_particle_integrations=0,new_meshes=0,physical_validation='PHYSICAL_VALIDATION_NOT_PERFORMED',events=classifications,
      quantitative_metrics=m,spatial_features=s,local_closed_balance=balance,ID822_increment_comparison=increments,
      recommendation='先处理背景场局部可信度，尤其1189近壁环境和822入口邻域的网格敏感性；本轮不自动实施。旧dt敏感性事实仍成立，但调积分器不能修复共点CFD差异。',
      minimal_missing_control={'count':1,'description':'若要升级到可信轨迹，唯一最小科学缺项是同一原几何/历史BC事件坐标上的一个已证明局部场误差受控的参考；现有中细档尚不具备此资格。此项仅建议，不授权新网格或新求解。'},
      immutable_production='Formal M0, old1565(458/527/437,141censored,2failures), fixed-window NOT_RUN, pipe scientific FAIL, animations and compatibility links unchanged.',
      evidence_layers={'A':'Input/native/mapping verified for the admitted payloads','B':'Same-position sensitivity and local closed-CV residual measured; no universal accuracy PASS','C':'PHYSICAL_VALIDATION_NOT_PERFORMED'})
    put('FIELD_RELIABILITY_RESULTS.json',result)
    text=['# 原 BraVa 急变区域背景场数值核查','',
      '2026-10-05。**原场不能无条件视为可靠。三例近壁事件的后段快流得到有限的跨网格支持；1189 的离壁环境仍有具体局部缺陷和敏感性。822 的原急增幅值、增量方向及固定上升尺度没有得到空间稳定性支持。应先处理背景场，再评估轨迹修复。**','',
      '|事件区域|原场问题|实际可比参考|同位置实际差异（mm/s）|结论|图|',
      '|---|---|---|---|---|---|']
    for pid in ['1189','1329','642','822']:
        e=m[pid];r=e['representatives'];a=[r[name]['A']['speed_mm_s'] for name in NAMES];b=[r[name]['B']['speed_mm_s'] for name in NAMES]
        if pid=='822':desc='增量 '+ ' / '.join(f'{r[name]["scalar_speed_increase_mm_s"]:.3f}' for name in NAMES);issue='168.656 µm 内原增速36.228';conclusion='幅度及尺度敏感，原急增未获支持';fig='04_ID822_SPATIAL_RISE.png'
        else:desc='末速 '+ ' / '.join(f'{v:.3f}' for v in b);issue='近壁慢→快';conclusion='后段快流重复；'+('离壁环境显著敏感、有局部缺陷' if pid=='1189' else '近壁导数仍敏感');fig='02_NEAR_WALL_FLOW.png'
        text.append(f'|{pid}|{issue}|固定原PLC中/细档|{desc}（原/中/细）|{conclusion}|[查看]({fig})|')
    text+=['','新增：共点场查询、全部线段分段、独立切面/法向线、闭合通量、五组300 dpi图。**新增 CFD=0；复用 CFD=3（原场＋2个独立研究解）；新增网格、微泡积分、population均=0。**没有真实匹配测量，物理层结论为 `PHYSICAL_VALIDATION_NOT_PERFORMED`。','',
      '## 直接回答五个问题','',
      '1. **1189：**后段约49–53 mm/s的快流存在于三份实际解；慢段流向、切向剪切及局部连续性残差不能视为已准确。细档没有消除所有朝壁面分量；有限壁距处的非零法向流也不能直接叫作穿壁。',
      '2. **1329/642：**共享“进入更快切向流层”的空间趋势，但不共享1189近乎纯法向主导的严重局部结构。代表末点的中细档速率差较小，不代表整窗u/G、离壁时机已收敛。',
      '3. **822：**原/中/细固定线段增速36.228 / 4.507 / 16.858 mm/s。三解总体朝相近下游方向，但增量向量及幅度不同；参考均未完成原固定10–90%水平的穿越，不能计算一个伪造的可比上升长度。',
      '4. **原机制：**近壁三例“固定历史路径进入快背景流”获得新增独立场支持；具体离壁动力学没有被重积分或验证。822仅保留一般的空间变化事实，原36.228 mm/s急增仍是未稳定离散场中的内部解释。',
      '5. **优先级：**先处理场的局部可信度。上轮822的dt敏感性不被推翻，但相同坐标在独立CFD间已不同，单改时间推进不能消除该输入问题。唯一最小科学缺项是同几何、同历史BC事件区域上的局部误差受控参考；当前没有该参考，本轮不自动补算。','',
      '## 参考准入、版本与稳态','',
      '三份完整 WALL/INLET/O1/O2/O3 的**全部节点坐标与三角形无序集合精确相同**，不是只验几处gap。固定PLC中、细档体网格分别762096、895647个四面体，原场692313。已有网格合法性审计及当前原始文件哈希同时保留。独立求解从零初态开始，均实际保存第81步；原生终态u/p与冻结载荷逐节点精确对应，canonical点坐标和四面体连接集合核验通过。原生文件缺少GlobalNodeID时使用已验证的精确坐标/连接对应，不宣称读到了不存在的编号。','',
      '三份solver.xml逐字节相同；求解器可执行、PETSc库SHA及实际PETSC_OPTIONS相同：8 MPI、CUDA矩阵/向量、FGMRES和ASM/LU。保留原P1/P1稳定化语义、dt=0.01 s、generalized-alpha谱半径0.5、rho=1056 kg/m³、mu=0.00345312 Pa·s。入口Flat、Impose_flux、Zero_out_perimeter，外向Q=−3e−7 m³/s；相同端口三角剖分与求解实施保持相同节点归一化。自然牵引出口参数556.97065554884477 / 0 / 17.137427153562019 Pa。没有重新校准或改为中心线边界。','',
      '|场|mesh SHA256|flow SHA256|最后保存前区间E_u|', '|---|---|---|---|']
    for name in NAMES:
        r=c['fields'][name];text.append(f'|{name}|`{r["raw_payload"]["mesh_sha256"]}`|`{r["raw_payload"]["flow_sha256"]}`|{r["steady"]["last_five_intervals"][-1]["E_u"]:.4g}|')
    text+=['','三个场均达到既有连续5个区间稳态及物理门槛，真实执行exit_code=0，记录中无失败线性解；这些只支持稳态执行，不是局部精度证明。未构造忽略跨面黏性通量的P1强残差；含稳定化的局部分布弱残差为 **NOT_EVALUABLE**。详见 [FIELD_COMPARABILITY.json](FIELD_COMPARABILITY.json)。','',
      'GitHub所指定分支的实际远端头本轮只读查询为 `d9ece125ef92ac87b87626fdf2a9e474fdb071f2`。开始时同步工作树本地头为 `bac63597e`，已有上轮原因报告的本地提交；BraVa工作树为 `1ca0cc6`。因此“服务器/本地证据已完成”与“GitHub已发布”分开登记。本轮交付保留独立本地提交，不宣称已推送。','',
      '## 固定范围、空间权重和覆盖','',
      '协议在读取参考速度载荷之前冻结，见 [EVENT_FIELD_PROTOCOL.json](EVENT_FIELD_PROTOCOL.json)。四窗完全沿用上轮；822出生前不造点。近壁慢段截止于既有FACE→EDGE时刻，释放段截止于首个h≥a的原接受点；822分出生→A、A→B、B后。原时间只标历史坐标，不表示稳态场发生时间变化。ROI固定为原窗穿越单元加两层共面邻接，无额外扩展。','',
      '四窗3576个接受位置、3572条真实线段全部保留。每份网格枚举真实过面点，合并所有分段并在共点求原节点P1 u、完整G、curl；不平均G、不配对不同网格的cell编号。空间主统计按真实弧长加权，历史暴露另按原dt加权。两点Gauss对每个共同P1分段的速度向量平方差准确；速率、方向及含变化法向的派生量额外保存四点细化对照。','',
      '所有三场的历史接受点和共同路径查询均有效；没有删除或移动旧点。切面格点在原腔外的共同空白保留为无效，不混作参考缺失。法向线在原有限面片投影不再适用处截断，未延长成假想平壁；保留完整坐标、无效掩码和原因。沿路径d是中心壁距，h=d−a；2nm是旧球面gap，绝非CFD探针统一壁距。','',
      '路径使用同一原有限三角几何法向，并保存FACE/EDGE/VERTEX与并列面片；全部接受点和法向线另枚举等距面片及法向歧义。其它切面/积分点保留有限特征歧义标记，不宣称每个切面点的所有并列面片均已枚举。822用固定路径基底解释流向。图02平面法向分量是33200原面片固定法向投影；法向剖面分别使用各自原脚点法向。','',
      '## 代表位置：背景速度与近壁输入','',
      '|事件|场|A/B速度 mm/s|A法向分量 mm/s（822为路径分量）|A剪切模 s⁻¹|A散度 s⁻¹|A最小高度 µm|',
      '|---|---|---|---|---|---|---|']
    for pid in map(str,PIDS):
        for name in NAMES:
            r=m[pid]['representatives'][name];a=r['A'];text.append(f'|{pid}|{name}|{a["speed_mm_s"]:.6f} / {r["B"]["speed_mm_s"]:.6f}|{a["normal_or_path_component_mm_s"]:.6f}|{norm(a["shear_s_inv"]):.3f}|{a["div_s_inv"]:.3f}|{a["cell_geometry"][2]*1e6:.3f}|')
    text+=['','最小高度是当前所属四面体的最小面高，并不一律等于wall-normal高度；六条边的最小/最大长度及面高范围、体积和原壁距见逐点NPZ的cell_geometry及机器结果。s=(I−nnᵀ)Gn；完整G、curl、nᵀGn、向量u_t亦保存。','',
      '1189原A点朝壁面分量占速率约99.49%，中/细约76.96%/46.90%；强朝壁面方向没有全消失，但“几乎纯法向”结构不再稳定。三个零壁节点单元的η与法向占比是代数相关量，不当作两份独立证据。A射线有效至d=116 µm，三场穿过3/8/20个单元；切向分层和逐单元导数变化很大，不能把曲线连续或no-slip节点正确当局部精确。B射线及后段几十mm/s流速更一致。','',
      '1329和642的A点法向分量仅占速率约几个百分点，切向流占主导。中细档B点速率差分别约0.0059和0.2632 mm/s；但整个快段的向量RMS差仍分别约2.811和3.636 mm/s，不能用单个相近端点覆盖整窗差异。','',
      '## 分阶段中细档差异','',
      '|事件/阶段|弧长u向量RMS差 mm/s|时间u向量RMS差 mm/s|方向RMS差 °|G Frobenius RMS差 s⁻¹|s RMS差 s⁻¹|',
      '|---|---|---|---|---|---|']
    pair=NAMES[1]+'__'+NAMES[2]
    for pid in map(str,PIDS):
        for label,r in m[pid]['phase_comparisons'][pair].items():
            z=r['arc_weight_m']['metrics'];zt=r['time_weight_s']['metrics'];text.append(f'|{pid}/{label}|{z["u"]["RMS_difference"]*1000:.6f}|{zt["u"]["RMS_difference"]*1000:.6f}|{z["direction_deg"]["RMS_difference"]:.4f}|{z["G"]["RMS_difference"]:.3f}|{z["shear"]["RMS_difference"]:.3f}|')
    text+=['','这些是可比解的差异，不是对真值的绝对误差。相对L2的分母是同组首个场的加权向量范数，JSON显式保留分母；未逐分量除近零值，未引入统一5%或8%硬门槛。','',
      '独立平面（均匀空间权重，原域内点）的中细档速度向量RMS差如下；它们不会被长时间贴壁的接受记录密度支配。','',
      '|事件|原域内格点/总格点|u RMS差 mm/s|G RMS差 s⁻¹|', '|---|---|---|---|']
    for pid in map(str,PIDS):
        p=s[pid]['independent_plane'];z=p['reference_comparisons'][pair]['metrics'];text.append(f'|{pid}|{p["original_in_domain_points"]}/{p["grid_points"]}|{z["u"]["RMS_difference"]*1000:.4f}|{z["G"]["RMS_difference"]:.3f}|')
    text+=['','## 增长发生在哪里','',
      '以下只是固定原路径上的空间增量核算，单元编号仅在各自场内解释。逐单元贡献为该段末/初的速率差，有正有负，不是唯一物理因果百分比；所有单元的有符号和与A/B速率差核对。','',
      '|事件|A→B实际路径长度 µm|原场最大三个正贡献单元：起始弧长/覆盖长度 µm，速率增量 mm/s|', '|---|---|---|']
    for pid in map(str,PIDS):
        e=s[pid];cells=e['fields']['M0']['all_cells_signed_gain'][:3]
        desc='；'.join(f'{v["cell"]}: {v["first_arc_from_A_um"]:.2f}/{v["path_length_um"]:.2f}, +{v["signed_speed_gain_mm_s"]:.3f}' for v in cells)
        text.append(f'|{pid}|{e["event_path_length_um"]:.4f}|{desc}|')
    text+=['','1189的大幅增长并非全部发生在388769：最主要两个正贡献位于644093与644095，释放后几十至百余微米区间；之后还存在正负变化。1329、642的增长分散在较长切向路径。完整同位置各场分解见 [SPATIAL_FEATURES.json](SPATIAL_FEATURES.json)。','',
      '822原四单元为402304→400899→59816→212669，其中400899在约69.20–130.10 µm区间贡献+26.993 mm/s；59816有小幅负贡献，原段不是严格单调上升。Δu=ΣGΔx的原场向量误差约6.25e−17 m/s，只证明内部算术。','',
      '|822场|A→B向量增量 mm/s|标量速率增量 mm/s|端点流向转角 °|原固定10–90%长度 µm|', '|---|---|---|---|---|']
    for name in NAMES:
        r=m['822']['representatives'][name];rise=m['822']['spatial_rise'][name]
        length=f'{rise["length_10_90_m"]*1e6:.4f}' if rise['length_10_90_m'] is not None else '不可评价：本段没有完成原水平穿越'
        text.append(f'|{name}|{[round(v,6) for v in r["delta_u_mm_s"]]}|{r["scalar_speed_increase_mm_s"]:.6f}|{r["direction_change_deg"]:.6f}|{length}|')
    text+=['',f'参考增量方向相对原增量分别差{increments[NAMES[1]]["angle_to_M0_increment_deg"]:.3f}°、{increments[NAMES[2]]["angle_to_M0_increment_deg"]:.3f}°；不要把“流动总体方向相近”改写成“急增向量一致”。中细档标量增量相差12.351 mm/s，约为原事件增量的34.09%。原10%/90%固定水平为340.9356 / 369.9180 mm/s；两个参考从高于低水平开始、均未到达高水平，不能用各自端点重设阈值。','',
      '822 A/B距相同原INLET平面分别66.743 / 235.375 µm，位于声明的Flat入口人工接管附近；这只定位边界邻域，不把它解释成已测生理加速。图04同时给固定线段与独立局部平面。','',
      '## 同一闭合控制体的局部连续性残差','',
      '近壁控制体以各事件A所属原四面体为种子，822以关键四单元为种子，各扩展0/1/2层共面邻居。使用**同一个原物理闭合表面**比较三场，计入全部侧面与壁面；如触及端口，保留相同端口的精确节点P1 trace。它们不是两张小圆盘，也不宣称是封闭通道的完整横截面。','',
      '原场按三节点面均值做精确积分；参考场按每面1、4、16、64、256、1024、4096个等面积子三角重心逐级积分。最后逐面L1细化差是数值敏感性估计，不是严格误差上界；各级结果与逐面通量全保存。最终36个场/控制体组合没有缺失积分点。','',
      '|事件/层|原Qnet m³/s|中Qnet m³/s|细Qnet m³/s|原/中/细 Qnet / throughput|细档末级面L1变化 m³/s|', '|---|---|---|---|---|---|']
    for pid in PIDS:
        for layer in [0,1,2]:
            rr=[next(z for z in balance[name] if z['event']==pid and z['layer']==layer) for name in NAMES]
            q=[r['levels'][-1]['net_flux_m3_s'] for r in rr];v=[r['levels'][-1]['throughput_m3_s'] for r in rr];ratio=' / '.join(f'{a/b:.5f}' for a,b in zip(q,v))
            text.append(f'|{pid}/{layer}|{q[0]:.6e}|{q[1]:.6e}|{q[2]:.6e}|{ratio}|{rr[2]["final_face_L1_refinement_difference_m3_s"]:.3e}|')
    text+=['','throughput定义为0.5Σ|每个原CV三角面通量|，保留原面分组，未加epsilon。1189原单元Qnet/throughput≈−1.896，是严重局部连续性缺陷；扩大两层后仍约−0.618。细档分别降到−0.228与−0.0809，虽减弱但仍未消失；实际差远大于积分细化变化。其它事件的净通量可能变号、邻域抵消，不能从更大CV更小归一化值推断粒子所在单元输入已准确。','',
      'ΣQface=div·V在原场的精确对应只检验实现；连续性应接近零的偏离、分辨率变化、CV尺度变化分别报告。原138584四壁点全零负例身份仍成立，仅作回归，不新增失败轨迹研究。','',
      '## 执行、测试、可视化与限制','',
      '实际服务器资源重新读取并保存在RESOURCES.json。按场串行加载，使用原CUDA FP64 FEM/有限壁面查询，没有新GPU kernel；CPU执行几何裁剪、积分汇总和绘图。任务峰值GPU张量分配约0.35 GiB，不能将利用率/CPU-GPU一致当精度证据。所有新计算纳入项目监督器，失败尝试也计入90分钟诊断预算。没有新增CFD，未触发4小时补算条件。','',
      '工程问题及处理：VTK FindCellsWithinBounds在细档极短路径包围盒上停滞，三次受监督停止；FP64容器不能消除停滞，故不把问题归因为单精度。诊断改为scipy保守包围盒候选＋原重心裁剪，并与已有原/中档分段及暴力包围盒对照。另将已知原面/相同边界的积分改为精确P1 trace，内部临界包含使用扩展精度检查但不移动点。初次空间特征输出的整数/字符串键错误已修复。每次运行的日志、退出及耗时详见receipt。','',
      '针对性测试包括线性u/G/curl、SI单位与刚体逆变换、同场零差、解析闭合通量、分片积分细化、候选完整性、严格域外/缺失分母、真实过面Δu、图表源数与300dpi、原输入不变性。测试通过只验证诊断实现。详见 [TEST_RESULTS.json](TEST_RESULTS.json) 与 [EXECUTION_RECEIPT.json](EXECUTION_RECEIPT.json)。','',
      '[全部图文入口](OPEN_RESULTS.html) · [01定位](01_EVENT_REGIONS.png) · [02近壁](02_NEAR_WALL_FLOW.png) · [03固定路径](03_FIXED_PATH_FIELD_COMPARISON.png) · [04 822空间增速](04_ID822_SPATIAL_RISE.png) · [05闭合通量与证据](05_LOCAL_BALANCE_AND_EVIDENCE.png)。所有图已检查；相机、平移、切面、箭头比例和共同色标见 [FIGURE_METADATA.json](FIGURE_METADATA.json)。ParaView ROI在vtk/，u是原节点P1，G/div/curl是原单元量，未做cell-to-point平均。','',
      '保持原科学状态：旧1565条仍为458/527/437出流、141观察截止、2数值失败；固定1.240815660 s＋10 s协议没有新的正式大批量结果。Gmsh V2科学FAIL、旧动画、退休路线与兼容链接均不改写。参考不是解析真值；没有系统加密比和渐近证据，因此不计算收敛阶、Richardson或GCI。','',
      '方法层级参见 [NASA/NPARC Verification Assessment](https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html) 与 [Spatial Convergence](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html)：实现核验、计算精度证据和物理实测验证应区分；本轮的两解差异不冒充误差界或真值。']
    (HERE/'FIELD_RELIABILITY_REPORT_ZH.md').write_text('\n'.join(text)+'\n')
    print('REPORT COMPLETE')

if __name__=='__main__':main()
