"""Run in WSL with standard Python after result synchronization."""
from pathlib import Path
import json,hashlib,difflib,datetime
root=Path(__file__).resolve().parents[2];out=root/'verification/clot_multifragment_crushing_mechanism_v1/20261016T222934Z';doc=Path('/home/lzy/projects/clot_rupture.md');folder=out/'wsl_document';folder.mkdir(exist_ok=True)
def read(p):return json.loads(p.read_text())
def sha(data):return hashlib.sha256(data).hexdigest()
assert read(out/'DELIVERY_VALIDATION.json')['status']=='PASS'
cases={k:read(out/'cases'/k/'ANALYSIS.json') for k in ['uniform','heterogeneous']};resources=read(out/'RUN_RESOURCES.json');raw=doc.read_bytes();text=raw.decode();assert '多块动态压碎对照（20261016T222934Z' not in text
backup=folder/'clot_rupture.before.md'
if not backup.exists():backup.write_bytes(raw)
if backup.read_bytes()!=raw:backup=folder/'clot_rupture.immediate_before.md';assert not backup.exists();backup.write_bytes(raw)
parts=['**多块动态压碎对照（20261016T222934Z，2026-10-17）**：本轮目标由最大单块转为多个20–200节点、保持内部连接的真实游离块。完成均匀材料和唯一冻结空间相关场的两组原生LAMMPS动态对照。沿用1792材料点、336固定底部锚点、95666原始键和2736原预切；没有为了得到目标碎块而增切、删键、换种子或修改行程。',
 '新增当前坐标下的球形压头接触：原生fix indent，半径182.818µm，投影面积等于上一轮局部压入42面片区域。初始间隙25µm，0–0.2ms接近，0.2–1.2ms压入240µm，保持至1.4ms，2.4ms撤回初始位置，继续观察到3.6ms。接触势能为ΣKindδc³/3，力为其空间负梯度；Kind约57973.6N/m²，使0.1dx重叠处切线刚度匹配原粒子接触。原生物理力转换后再施加压头力；未接触或释放后力为0，脱离碎片再次碰撞则可重新受力。该KINEMATIC_IMPACT_SURROGATE没有真实锤质量，也不是已验证的微泡冲击。',
 '非均匀性只作用于原始bulk键ID的s0：seed=20261017，lc=250µm，平滑随机场按bulk最大绝对值归一化，s0_e=0.15(1+0.10q_e)，范围0.135–0.1647022。sc=0.18、G/K、weighted volume、原始体积、界面和接触参数不变。lc/dx=5、lc/δ≈1.66；α是最大扰动而非标准差。这是未标定的机制假设，改变s0也可能改变耗能，不是固定Gc实验。新库按原始ID存储并保存restart；初始弹性、变分力、排序/restart、不可逆历史及移动接触功的关键检查通过，旧库未替换。']
for name,label in [('uniform','均匀'),('heterogeneous','非均匀')]:
 q=cases[name];z=q['final'];fraction='未定义（无脱离）' if z['debris12_detached_fraction'] is None else f"{100*z['debris12_detached_fraction']:.5g}%"
 parts.append(f"{label}材料：首软化{q['first_softening_s']}s，首完整失效{q['first_failure_s']}s；历史最多20–200点自由块{q['history_max_N20_200']}个，末态{int(z['N20_200'])}个，其占初始质量{100*z['mass_fraction20_200']:.5g}%；历史/末态最大块{q['history_max_free_nodes']}/{int(z['largest_nodes'])}点；末态最大块占初始材料{100*z['largest_fraction']:.5g}%；1–2点碎屑占脱离质量{fraction}。峰值压头反力{q['peak_reaction_N']*1e6:.6g}µN，最大重叠{q['max_overlap_dx']:.6g}dx，实际外功{z['W_J']*1e9:.7g}nJ。结果为{q['exploratory_target']}，分类{', '.join(q['classification'])}。")
parts.append('失效键网络按共享原始端点追踪出生与合并，排除初始预切。它与存活键定义的材料分量分离是两种诊断，不能把断键点云绘成未经证实的连续裂面。每个接受态保留分量原始ID、体积等效尺寸、内部存活率、平均损伤、质心速度和内部相对动能；接触对不作为材料连接。')
for name,label in [('uniform','均匀'),('heterogeneous','非均匀')]:
 q=cases[name];parts.append(f"{label}工况：撤回结束时中等自由块{len(q['withdrawal_medium_cohort'])}个，其中{q['withdrawal_medium_cohort_surviving']}个原成员组保持连通至末态；撤回后再接触={q['contact_after_scheduled_withdrawal']}，末尾连续无接触观察窗{q['actual_contact_free_tail_s']}s，撤回后原自由分量二次裂解{q['secondary_after_withdrawal']}次。仅称CONNECTED_THROUGH_OBSERVED_WINDOW，不称永久稳定。")
parts+=['两组是相同运动学轨迹，不是等反力或等能量对照；与旧1000Pa参考牵引也不构成单因素的等功比较。外功按接受态之间压头路径积分，账本为W=ΔU+ΔK+D+ΔUc−Jc+ΔUi+R；没有重复加入材料接触功或数值阻尼，残差仍保留为残差。压力压入通过侧向拉伸/几何剪切间接触发当前正键伸长损伤律；起裂附近保留局部Green应变、非仿射程度及材料键virial/参考体积诊断，不能据此宣称已认证多轴血栓强度。',
 f"两组原生累计墙钟{resources['native_wall_s']:.3f}s，CPU{resources['native_cpu_s']:.3f}s，单进程峰值RSS{resources['peak_RSS_bytes']/1024**3:.5f}GiB；1rank×1线程，顺序运行。未增加依赖、完整V3、CFD、随机种子或行程扫描。保留源码diff、SHA、必要FP64接受态、逐步损伤记录与失败预检。",
 '下一轮先审查局部损伤网络、实际接触重叠及卸载后演化，再判断是否需要有标定依据的剪切/压缩混合失效机制。单一行程、单一随机场不足以区分所有失败原因，也不足以直接迁移真实1927节点V3。20点碎片不等于充分解析的200µm实体；目标数量不是医学安全标准。真实微泡、微射流或超声到此压头轨迹、反力和功的物理联系仍缺失。',
 '见[八组主图与完整末态投影](temp_storage/clot_github_sync_20261015/clot_clearance/results/clot_multifragment_crushing_mechanism_v1/20261016T222934Z/OPEN_RESULTS.html)、[完整科研报告](temp_storage/clot_github_sync_20261015/clot_clearance/verification/clot_multifragment_crushing_mechanism_v1/20261016T222934Z/REVIEW_CLOT_MULTIFRAGMENT_CRUSHING_MECHANISM_V1.md)、[真实碎块分布](temp_storage/clot_github_sync_20261015/clot_clearance/verification/clot_multifragment_crushing_mechanism_v1/20261016T222934Z/MULTIFRAGMENT_SIZE_DISTRIBUTION.csv)。本轮完成后停止，等待审查。']
block='\n\n'.join(parts)+'\n\n';anchor='# ⑥ 外加剪切';index=text.index(anchor);new=text[:index]+block+text[index:];assert doc.read_bytes()==raw;doc.write_text(new);after=doc.read_bytes();assert after.decode()==new and new.replace(block,'',1)==text
(folder/'clot_rupture.after.md').write_bytes(after);(folder/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),new.splitlines(True),fromfile='clot_rupture.before.md',tofile='clot_rupture.after.md')))
(out/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(dict(status='PASS',path=str(doc),before_backup=str(backup),before_sha256=sha(raw),after_sha256=sha(after),history_preserved_exactly=True,readback_equal=True),ensure_ascii=False,indent=2)+'\n');print('WSL_DOCUMENT_UPDATED',sha(after))
