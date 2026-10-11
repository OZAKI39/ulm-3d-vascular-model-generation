from analysis.dual_track_v1.common import *
import csv,shutil

def main():
 geom=json.loads((B/'GEOMETRY_AUDIT.json').read_text());base=json.loads((B/'baseline/CFD_RESULT.json').read_text());assert base['status']=='PASS_WITH_SCOPE';levels=['baseline'];finefile=B/'refined/CFD_RESULT.json'
 if finefile.exists() and json.loads(finefile.read_text())['status']=='PASS_WITH_SCOPE':levels.append('refined')
 chosen=levels[-1];d={l:json.loads((B/l/'CFD_RESULT.json').read_text()) for l in levels};t={l:json.loads((B/l/'CFD_CLOT_TRACTION_AUDIT.json').read_text()) for l in levels};c={l:json.loads((B/l/'CFD_BASELINE_COMPARISON.json').read_text()) for l in levels};meshes={l:json.loads((B/l/'MESH_QUALITY_AUDIT.json').read_text()) for l in levels};sens=json.loads((B/'MESH_SENSITIVITY.json').read_text()) if (B/'MESH_SENSITIVITY.json').exists() else dict(status='NOT_ACCEPTED_OR_NOT_RUN');r=d[chosen];tr=t[chosen];comp=c[chosen];old=comp['old_metrics'];metrics=r['metrics'];good=sens['status']=='PASS_ENGINEERING_SCREEN'
 for name in ['CFD_CLOT_TRACTION_AUDIT.json','CFD_MASS_BALANCE.json','TRACTION_INTERFACE.json','TRACTION_MATH_AUDIT.json','CLOT_SURFACE_TRACTIONS.npz']:shutil.copy2(B/chosen/name,B/name)
 save(B/'SELECTED_FLOW.json',dict(level=chosen,selection='Use accepted local-refined field when available; mesh sensitivity remains a separate scientific status',snapshot=r['snapshot'],snapshot_sha256=r['snapshot_sha256'],mesh_sha256=r['mesh_sha256'],surface_data_sha256=sha(B/'CLOT_SURFACE_TRACTIONS.npz')))
 states=dict(CLOT_PRESENT_GEOMETRY_VALIDITY='PASS_ORIGINAL_V3_ANALYTIC_DEFINITION',CFD_MESH_VALIDITY='PASS_CLOSED_LABELLED_POSITIVE_VOLUME',STEADY_FLOW_SOLVER='PASS_NATIVE_LINEAR_NONLINEAR_AND_STEADY_GATES',MASS_BALANCE='PASS_FROZEN_1E_MINUS5',CLOT_SURFACE_PRESSURE='EXTRACTED_GAUGE_PRESSURE',CLOT_SURFACE_SHEAR_TRACTION='EXTRACTED_P1_ADJACENT_ELEMENT_GRADIENT',INTEGRATED_FLUID_FORCE='PASS_EXPOSED_SURFACE_INTEGRAL_AND_MATH_AUDIT',MESH_SENSITIVITY=sens['status'],BASELINE_COMPARABILITY='PASS_SAME_SUPPLY_AND_NATIVE_Y_WITH_STATED_MESH_LIMITS',STREAMING_NOT_INCLUDED=True,TWO_WAY_FSI_NOT_INCLUDED=True)
 save(B/'FINAL_STATUS.json',dict(run_id=RUN,states=states,selected_level=chosen,formal_CFD_attempts=sum((B/l/'FORMAL_ATTEMPT.json').exists() for l in ['baseline','refined']),accepted_flows=len(levels),stopped=True))
 rows=[]
 keys=[('Q_in_m3_s',lambda mm,tt,cc:mm['Q_in_m3_s']),('max_velocity_m_s',lambda mm,tt,cc:mm['velocity_max_m_s']),('outlet_1_fraction',lambda mm,tt,cc:mm['outlet_fractions']['OUTLET_A']),('outlet_2_fraction',lambda mm,tt,cc:mm['outlet_fractions']['OUTLET_B']),('inlet_pressure_Pa',lambda mm,tt,cc:mm['area_average_pressure_pa']['INLET']),('mother_pressure_drop_Pa',lambda mm,tt,cc:cc['parent_pressure_drop'][1]['drop_Pa']),('surface_pressure_mean_Pa',lambda mm,tt,cc:tt['pressure_Pa']['mean']),('surface_shear_mean_Pa',lambda mm,tt,cc:tt['shear_Pa']['mean']),('surface_shear_P95_Pa',lambda mm,tt,cc:tt['shear_Pa']['p95']),('total_force_norm_N',lambda mm,tt,cc:tt['total_force_norm_N'])]
 for name,fun in keys:
  row=dict(metric=name,unobstructed='NOT_APPLICABLE_NO_CLOT_SURFACE')
  if name in ['Q_in_m3_s','max_velocity_m_s','outlet_1_fraction','outlet_2_fraction','inlet_pressure_Pa']:row['unobstructed']=fun(old,None,None)
  if name=='mother_pressure_drop_Pa':row['unobstructed']=comp['parent_pressure_drop'][0]['drop_Pa']
  for level in levels:row[level]=fun(d[level]['metrics'],t[level],c[level])
  rows.append(row)
 with (B/'CFD_BASELINE_COMPARISON.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=['metric','unobstructed']+levels);w.writeheader();w.writerows(rows)
 telemetry={}
 for l in levels:
  records=json.loads((B/l/'RESOURCES.json').read_text());gpu=[]
  for rr in records:
   value=rr.get('gpu_util_percent_memory_MiB')
   if value:
    try:gpu.append([float(x) for x in value.split(',')])
    except ValueError:pass
  telemetry[l]=dict(native_runtime_s=d[l]['runtime_s'],process_tree_peak_sampled_RSS_MiB=max(v['RSS_bytes'] for v in records)/1024**2,maximum_sampled_live_CPU_s=max(v['CPU_seconds_live'] for v in records),sampling_interval_s=15,CPU_scope='Process-tree sampled CPU lower bound; exited processes may not appear in the final sample',GPU_max_sampled_util_percent=max(v[0] for v in gpu) if gpu else None,GPU_mean_sampled_util_percent=float(np.mean([v[0] for v in gpu])) if gpu else None,GPU_peak_sampled_memory_MiB=max(v[1] for v in gpu) if gpu else None)
 save(B/'RESOURCE_SUMMARY.json',telemetry)
 summary=dict(selected_level=chosen,geometry=geom,flows=d,traction=t,comparison=c,mesh=meshes,mesh_sensitivity=sens,telemetry=telemetry,formal_CFD_attempts=sum((B/l/'FORMAL_ATTEMPT.json').exists() for l in ['baseline','refined']),total_accepted_native_runtime_s=sum(v['runtime_s'] for v in d.values()),scope='Rigid original synthetic V3, Newtonian water, original fixed inlet flux and zero traction outlets; no ultrasound/streaming/FSI')
 save(B/'ANALYSIS_SUMMARY.json',summary)
 mt='\n'.join(f'| {l} | {meshes[l]["points"]:,} | {meshes[l]["tetrahedra"]:,} | {d[l]["stop"]["step"]} | {d[l]["metrics"]["epsilon_mass"]:.6e} | {d[l]["runtime_s"]:.3f} |' for l in levels)
 sensitivity=('；'.join(f'{k}：{v*100:.4f}%' for k,v in sens['relative_errors'].items()) if 'relative_errors' in sens else '本轮未获得可接受的两网格对照；见停止证据。')
 nexttask='验证含栓 CFD 暴露面牵引向现有 LAMMPS 材料点的离线保力、保力矩映射，先核对单位、坐标、面身份和面积权重，再决定是否接入加载。' if good else '针对本轮未收敛的表面牵引指标设计一次有明确误差预算的近壁梯度/局部网格验证，先解决载荷分辨率，再进行材料点映射。'
 save(B/'NEXT_PRIORITY.json',dict(recommendation=nexttask,automatic_execution=False,basis=sens['status']))
 text=f'''# Y_CLOT_PRESENT_CFD_BASELINE_V1 科研审查

运行：`{RUN}`。已真正求解原生 Y 血管中的原 V3 含栓稳态背景流。独立身份为 `CLOT_PRESENT_RIGID_GEOMETRY_STEADY_CFD_BASELINE`。下列牵引数值来自 `{chosen}` 接受态；网格敏感性状态为 `{sens['status']}`。

## 用户要求的 12 项回答

1. **是否为原 V3。** 是。直接读取原 `y_clot_morphology_v3.profile`、冻结配置和原生 WALL 三角面，在原壁面沿径向向管腔内施加同一个多叶高度函数。没有从 1,927 个点做凸包，没有用球或椭球替代，也没有平滑形状。x 中心约 −6 mm，位于母管 +z 内壁；图中的主体仍是原材料点。

2. **体积与流体域。** 原连续表面体积 {geom['original_V3_continuous_volume_m3']*1e9:.12f} mm³，新 CFD 有限三角化体积 {geom['CFD_triangulated_clot_volume_m3']*1e9:.12f} mm³，相差 {geom['relative_to_original_V3_mesh']*100:.9f}%；材料点积分体积 0.240875 mm³，属于不同口径。没有为匹配体积调整参数。原生管腔体积 {geom['original_lumen_volume_m3']*1e9:.9f} mm³，含栓流体域 {geom['fluid_volume_m3']*1e9:.9f} mm³，差值等于进入腔内的血栓体积。高度函数从实际壁面向内构造，未发生另一个体积的墙外裁切。接触线共 293 个共享节点、间隙为零；闭合单连通流体体网格通过。最大向内突出约 {geom['max_inward_height_m']*1e6:.6f} µm，29 个真实截面记录于 `GEOMETRY_AUDIT.json`。

3. **暴露面。** `ModelFaceID=5` 的 3,611 个面元是流体真正接触的表面，面积 {tr['exposed_area_m2']*1e6:.9f} mm²；内部壁面附着面不计入流体边界或载荷。由相邻四面体确认流体外法向，取 n_clot＝−n_fluid；流体对血栓牵引为 σn_clot，σ＝−pI＋µ(∇u＋∇uᵀ)。法向符号、剪切与法向正交及力矩积分均有独立数学检查。

4. **网格与边界。** 使用已有 SimVascular 2023.05 / TetGen，保留完整原生 Y：母管半径 1.3 mm、长 14 mm，两支半径 1 mm、长 12 mm、±30°。物理标签为血管壁、入口、两个出口、血栓暴露面。全局目标边长 250 µm，血栓附近球形区域体网格目标 100/70 µm；两次体网格共享同一提前细分的 CFD 表面。所有四面体体积与符号质量为正，边界无遗漏或重复；没有另外生成棱柱边界层。入口保持原 Parabolic、Impose_flux=true、Zero_out_perimeter=true，Q＝2.6354060534101725e−8 m³/s，即约 1.581244 mL/min，平均速度 5 mm/s。出口保持原 Neu、Value=0（零牵引，不能解释成每个出口节点压力恰好为零）；两类壁面均 no-slip。ρ＝1000 kg/m³、µ＝0.001 Pa·s。

5. **是否实际求解。** 是。使用锁定 svMultiPhysics/PETSc/MPI 和原 GPU 选项，8 MPI、dt＝0.01 s、原广义 α 时间推进至稳态。速度、压力来自新的原生 VTU 和双精度重启记录，没有以旧场插值冒充含栓解。原无栓 CFD 未重算。原生速度与压力数组为 float64；原生 vtkPoints 序列化为 float32，按严格的一一对应和精确 float32 舍入身份映射到冻结输入的 double 坐标，所有积分使用输入几何。

6. **质量和稳态。** 冻结门槛为质量/入口误差≤1e−5，连续五个十步区间 E_u≤1e−5、E_Q≤1e−6，全部原生线性收敛，非线性相对残差≤1e−10，有限值和 no-slip 均通过。所选态质量不平衡 {metrics['epsilon_mass']:.9e}，入口误差 {metrics['epsilon_Q']:.9e}；两类壁面最大速度均为零，两个出口净流为正且无负向流出面积。入口负向外流比例为 1 是正常入流符号，不是倒流。

| 接受态 | 节点 | 四面体 | 稳态触发步 | 质量不平衡 | 原生墙钟 / s |
|---|---:|---:|---:|---:|---:|
{mt}

7. **有栓/无栓变化。** 旧最大速度 {old['velocity_max_m_s']*1e3:.9f} mm/s，新值 {metrics['velocity_max_m_s']*1e3:.9f} mm/s，变化 {(metrics['velocity_max_m_s']/old['velocity_max_m_s']-1)*100:.6f}%；新峰值坐标 {comp['new_peak_location_m']} m。旧两个出口比例为 {old['outlet_fractions']['OUTLET_A']:.9f}/{old['outlet_fractions']['OUTLET_B']:.9f}，新值 {metrics['outlet_fractions']['OUTLET_A']:.9f}/{metrics['outlet_fractions']['OUTLET_B']:.9f}。母管 x＝−9 至 −3 mm 的面积平均压降由 {comp['parent_pressure_drop'][0]['drop_Pa']:.9f} 增至 {comp['parent_pressure_drop'][1]['drop_Pa']:.9f} Pa；入口平均压力由 {old['area_average_pressure_pa']['INLET']:.9f} 变为 {metrics['area_average_pressure_pa']['INLET']:.9f} Pa。公共中心线和切面只做严格 P1 分析插值，域外留空不外推；几何/网格差别和插值误差与真正 CFD 解明确区分。

8. **表面压力与剪切。** 面积平均表压 {tr['pressure_Pa']['mean']:.9f} Pa，面积加权 P95 {tr['pressure_Pa']['p95']:.9f} Pa；面积平均切向剪切 {tr['shear_Pa']['mean']:.9f} Pa，P95 {tr['shear_Pa']['p95']:.9f} Pa，面元最大值 {tr['shear_Pa']['max']:.9f} Pa。压力为相对原出口参考的表压，不是绝对血压。黏性牵引来自相邻四面体的原生 P1 速度梯度，未做平滑或替换模型；压力法向、黏性法向、切向和总牵引分别保存。

9. **积分流体力。** 压力力＝{tr['pressure_force_N']} N，黏性力＝{tr['viscous_force_N']} N，总力＝{tr['total_force_N']} N，大小 {tr['total_force_norm_N']*1e6:.9f} µN。参考点为实际材料点质心 {tr['reference_point_m']} m，总力矩＝{tr['torque_Nm']} N·m。只统计流体暴露面，不含固体壁面黏附反力。接触线 100 µm 范围按实际边线段的精确距离识别，面积占比 {tr['contact_rim']['area_fraction']*100:.6f}%，牵引模面积积分占比 {tr['contact_rim']['traction_magnitude_integral_fraction']*100:.6f}%；全局剪切峰值是否位于该范围：{tr['contact_rim']['contains_global_shear_peak']}。

10. **一次网格对照。** 状态 `{sens['status']}`。局部体单元由 28,237 增至 38,127，中位体积由 3.921428e−14 降至 2.556330e−14 m³；表面几何、入口和出口条件保持同一冻结定义。相对变化：{sensitivity}。关键积分/统计量门槛为 15%，出口比例绝对变化门槛为 0.01；局部峰值单列。它只是一次局部体网格工程筛查，**不是表面三角化、渐近网格收敛或接触线极值已认证**；没有再加第三个网格。

11. **未来 LAMMPS 载荷来源。** `CLOT_SURFACE_TRACTIONS.npz` 保存全精度面元坐标、真实面积、相邻单元/节点 ID、法向、压力、分解牵引、积分力/力矩及积分点数据；`TRACTION_INTERFACE.json` 记录单位、坐标和表面身份。这些是候选背景流体载荷数据，尚未建立到 1,927 个材料点的保力/保力矩映射；没有把结果写入生产 LPS 载荷。

12. **为何不是 streaming/FSI。** 本轮是刚性 V3 障碍物上的牛顿水样背景流，没有微泡诱导 streaming、超声反馈、血栓变形或断裂后更新，也没有双向流固耦合。不能把静态剪切当作超声剪切，不能由 Pa 直接换算 PZT 功率，也没有患者级/临床认证。

## 接口问题与原始证据

正式 CFD 前的原生接口问题均保留在 `baseline_interface_failure_*`：int64 标签不兼容；旧 set_walls 逻辑在仅划分体网格时移除三个端盖；一次 SSH 同步中断导致未更新脚本被再次调用。通过原生 VTK 检查与源码定位后，改为 int32 标签、让体网格输入保留全部五类封闭边界；物理 no-slip 仍仅作用在壁面和血栓面。实际体网格导出单精度坐标后，以精确双向身份恢复冻结边界 double 坐标，重新检查全部几何门槛，未重划基线网格。

第一套正式 CFD 没有重启。旧只读监测器不能匹配原生 float32 坐标输出，使用任务独立辅助观察器按精确舍入身份读取同一运行，按原稳态门槛写 STOP_SIM，并在父进程回收前实际记录 MPI 启动器退出码 0。旧读取器的报告异常保留，最终状态来自实际原生日志、完整检查点、接受态及独立监测记录。第二套运行已直接使用修正读取器。冻结数值门槛未放宽。

## 审查与下一步

真实图索引为 `OPEN_RESULTS.html`：原生 Y/材料点、局部体网格/绕流截面、公共位置速度/压力、暴露面压力/剪切、质量/稳态/力/网格对照、三维流线与方向、横截面和剪切方向。没有新医学外壳或历史动画。

下一轮唯一优先建议：{nexttask} 本轮不自动执行。

原生数据全部位于本 run 目录；软件、输入 SHA、完整资源统计、精简回收清单及文档更新证据见共同目录。实际方法和接口依据见 [SimVascular TetGen API](https://simvascular.github.io/documentation/python_interface/modules/docs/meshing_TetGenOptions.html) 与 [原生接口源码说明](https://github.com/SimVascular/SimVascular/blob/master/Code/Source/sv/Mesh/TetGenMeshObject/sv_TetGenMeshObject.cxx)，运行时二进制身份另有锁定 SHA，不以当前上游源码代替原生验证。
'''
 (B/'REVIEW_Y_CLOT_PRESENT_CFD_BASELINE_V1.md').write_text(text)
if __name__=='__main__':main()
