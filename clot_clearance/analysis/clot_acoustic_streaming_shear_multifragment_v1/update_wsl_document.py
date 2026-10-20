"""Local stdlib-only update; exact historical bytes preserved outside new block."""
from pathlib import Path
import json,hashlib,difflib
root=Path(__file__).resolve().parents[2];task='clot_acoustic_streaming_shear_multifragment_v1';run='20261017T132706Z';out=root/'verification'/task/run
doc=Path('/home/lzy/projects/clot_rupture.md');folder=out/'wsl_document';folder.mkdir(exist_ok=True)
def read(p):return json.loads(p.read_text())
def sha(raw):return hashlib.sha256(raw).hexdigest()
assert read(out/'DELIVERY_VALIDATION.json')['status']=='PASS'
q=read(out/'ANALYSIS.json');status=read(out/'FINAL_STATUS.json');bridge=read(out/'CRACK_BRIDGE_AUDIT.json')['cases'];mapping=read(out/'SURFACE_MAPPING_AUDIT.json')['v3'];resources=read(out/'RUN_RESOURCES.json')
raw=doc.read_bytes();text=raw.decode();heading=f'**恢复声流／流体剪切主线（{run}，2026-10-17）**';assert heading not in text
backup=folder/'clot_rupture.before.md';assert not backup.exists();backup.write_bytes(raw)
parts=[heading+'：停止继续迭代刚性压头。此前单/双压头属于规定运动的机械法向接触代理，保留历史证据，但不是微泡声流或流体表面剪切。当前入口为acoustic_streaming_shear.py → 原始表面三角形 → native外力回调 → 不可逆κ/d → attachment-aware原始存活连接。',
 '实际重新核对F1/F10S/APECSS源数据：F1是其自身几何的声学解，F10S为2.5MHz均匀水控制代理；APECSS提供有效叶内体积谐波，质心轨迹只属于原坐标。没有当前血栓配准/间隙、近壁一阶速度与振荡边界、二阶平均流体解。STREAMING_LOAD_NOT_YET_DETERMINED，不能把声压、R(t)或气泡质心速度直接当成血栓表面声流应力。缓存黏度2mPa·s，本轮历史Stokes模板1mPa·s，二者并非同一套完整介质耦合。',
 '本轮唯一模式为PRESCRIBED_STOKES_SHEAR_TEMPLATE：u=(sin(2πy/700µm)v(z,t),0,0)，v_t=ν(v_zz−k²v)，25µm间隙、底壁静止、上边界4mm/s的连续反向剪切、初始静止、零压力增量。sinh稳态项及正弦扩散模态实际求解并以独立有限差分验证；该无限平行边界模型的无散和零对流严格成立，但没有泡驱动，也没有有限血栓边缘/变形反馈。',
 '保持1792材料点、336锚点、95666原始键、2736原预切，以及G/K/s0/sc不变。实际顶面原节点三角形以1/3非负权重将牵引转成N，保持总力/力矩/虚功。上升2ms，2.4ms下降，2.8ms停驱，自然黏性衰减并观察至3.6ms。冻结参考面一向加载会继续作用于移动或脱离原节点；大变形后只能称机械死牵引探索，不宣称当前真实Streaming。',
 f'真实V3面片只做映射审计：旧仿射映射权重最小/最大{mapping["signed_affine"]["min_weight"]:.6g}/{mapping["signed_affine"]["max_weight"]:.6g}，最大列绝对和{mapping["signed_affine"]["max_column_abs_sum"]:.6g}。新模式中{mapping["local_comparison"]["common_active_nodes"]}共同活跃点有{mapping["local_comparison"]["nodes_direction_difference_over_90deg"]}点与独立非负局部核合力夹角>90°，最大局部力预算比{mapping["local_comparison"]["signed_force_over_positive_local_force_budget_quantiles"][-1]:.6g}。这是实际合力比较，不能由负权重直接推断。非负替代牺牲一阶矩；没有简单裁零旧权重，没有完整运行V3。小试样P1三角积分的节点面积相对规则控制面积为2/3–4/3，按映射面积归一化剪切最大误差为峰值2.8112%，未做网格收敛认证。']
for name,label in [('stokes','A：Stokes近似'),('inverse','B：反向剪切需求')]:
 a=q[name];end=a['final'];m=a['final_metrics'];bridges=bridge[name]['final']['regions']
 parts.append(f'{label}：实际面片峰值{a["peak_shear_Pa"]:.9g}Pa，末态软化{int(end["softened"])}、本体/界面完整失效{int(end["bulk_failed"])}/{int(end["interface_failed"])}，真实游离{int(end["detached_nodes"])}点/{int(end["components"])}分量。历史最大游离块{a["history_max_free_nodes"]}点，20–200块历史最大/末态{a["history_max_N20_200"]}/{int(end["N20_200"])}，细屑/脱离质量={m["fine_debris_mass_fraction_detached"]}，最大块/脱离质量={m["largest_free_fragment_fraction_detached"]}。W={end["W_J"]*1e9:.9g}nJ、D={end["D_J"]*1e9:.9g}nJ，最大相对能量残差{a["max_relative_energy_residual"]:.9g}；三历史诊断区存活跨界键{[r["surviving_crossing_bonds"] for r in bridges]}，全部按原始存活拓扑核算。')
parts += [f'分层目标未改为最大断键或最大块：至少5个20–200节点自由分量、初始脱离质量≥10%、最大块≤脱离质量50%、细屑≤20%、内部真实存活连接和3.4→3.6ms有限近卸载观察。完整目标达成={status["five_fragment_target_met"]}；无脱离时比例保持null。B/A={status["template_traction_ratio"]:.8g}仅是同模式同几何时程的规定牵引倍数，不是声流不足倍数，更不是目标成功阈值。',
 'B的几何诊断在0.37725ms触发，早于0.88925ms首次软化；末态约99.65%外功转为动能，固定原始ID载荷继续加速已脱离材料，不能认证真实流体能量输入或碎裂效率。A未触发几何阈值也不代表运动壁面反馈或流体阻尼已验证。本轮建立了可复用的流体剪切到原生断裂路径，但声流源/边界、V3局部映射、大变形表面重算和材料校准仍限制物理结论。三维桥接连接能解释特定诊断区域尚未自由，不能从两工况唯一判定本构缺失混合模式。不得通过提高“预测声流”、人为删键或改材料达标。暂不启动微泡参数BoTorch闭环优化或真实V3；先补齐可信近壁流体源与耦合边界。',
 f'两组主计算合计墙钟{resources["native_wall_s"]:.3f}s、CPU{resources["native_cpu_s"]:.3f}s，峰值RSS{resources["peak_RSS_bytes"]/2**20:.3f}MiB；每例1rank×1线程。原生库无修改，dt0.25µs/事件最小0.03125µs，Δd≤0.1和1%能量门槛保留，R不作耗散；完整ID损伤重放、所有接受态分量与关键FP64状态通过回读。',
 f'见[本轮图表总览](temp_storage/clot_github_sync_20261015/clot_clearance/results/{task}/{run}/OPEN_RESULTS.html)、[20项科研问答](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{task}/{run}/REVIEW_CLOT_ACOUSTIC_STREAMING_SHEAR_MULTIFRAGMENT_V1.md)。本轮仅两例，无新完整V3、无大规模优化。']
block='\n\n'.join(parts)+'\n\n';i=text.index('# ⑥ 外加剪切');new=text[:i]+block+text[i:]
assert doc.read_bytes()==raw;doc.write_text(new);after=doc.read_bytes();assert after.decode()==new and new.replace(block,'',1)==text
(folder/'clot_rupture.after.md').write_bytes(after);(folder/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),new.splitlines(True),fromfile='clot_rupture.before.md',tofile='clot_rupture.after.md')))
(out/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(dict(status='PASS',path=str(doc),before_backup=str(backup),before_sha256=sha(raw),after_sha256=sha(after),history_preserved_exactly=True,readback_equal=True),ensure_ascii=False,indent=2)+'\n');print('WSL_DOCUMENT_UPDATED',sha(after))
