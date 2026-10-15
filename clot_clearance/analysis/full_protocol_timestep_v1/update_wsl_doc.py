"""Run in WSL only. Preserve every existing byte outside one inserted research block."""
from pathlib import Path
import hashlib,json,difflib,datetime
LR=Path('/home/lzy/projects/temp_storage/clot_github_sync_20261011/clot_clearance')
TASK='lps_native_softening_full_protocol_timestep_admission_v1';RUN='20261013T224838Z'
O=LR/'verification'/TASK/RUN
doc=Path('/home/lzy/projects/clot_rupture.md')
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
 before=doc.read_bytes();old=before.decode();assert RUN not in old
 checkpoint=sha(before)
 assert old.count('# ⑥ 外加剪切')==1
 # Read actual results before creating prose. Scientific numbers are assertions.
 a=json.loads((O/'FULL_PROTOCOL_NATIVE_TIMESTEP_AUDIT.json').read_text());assert a['executed_steps']==2562 and a['stop_reasons']==['MAX_DAMAGE_INCREMENT']
 assert a['maxima']['max_dd']==0.11363469187762688
 block=f'''**本轮0.25µs完整时程准入审查（{RUN}，2026-10-14）**：复用上一轮已验收的隔离原生库，SHA256仍为 **`866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd`**，没有修改或重新编译C++。341个节点、12047条原始键、49/47个夹持点及245个自由点、真实质量ρV、G/K/ρ、dx/horizon、s₀=0.15、s_c=0.18、α=0和完整参考m/r/taper均不变。按上一轮完全相同的半余弦解析位移，从**原始无损伤t=0**新运行唯一工况，dt=0.25µs，预注册2800步/700µs；没有使用105µs旧restart拼接，也没有其他时间步或新试样。

**完整协议仍为FAILED_ADMISSION。** 实际执行2562步，在卸载末段 **640.5µs** 因原始键 **1813** 的单步损伤增量 **Δd=0.11363469187762688>0.10** 首次越界立即停止；最后通过全部冻结门槛的状态是第2561步、**640.25µs**。越界步是原生已提交但未通过准入的停止现场，不能称为合格末态。已完成0–300µs升载及300–350µs最大位移保持，只覆盖350–640.5µs的部分卸载；尚未到达650µs卸载终点和650–700µs末保持。最终保持为NOT_RUN，不能写成完整加载—卸载基线通过。

| 实际观察或核查 | 本轮结果及科学边界 |
|---|---|
| 首次自然事件 | 第340步、**85µs**，原始键3393首次软化；第370步、**92.5µs**，2180、2194、3368、3393、3904、9197、10264共7键首次完全失效。 |
| 越界局部状态 | 键1813两端为自由节点24、245；κ从0.150007021399变为0.152903004587，d从0.000280842817变为0.113915534694。原生与冻结独立模型逐值一致，越界时两端额外接触力均为0；它是惯性运动下接近s₀的局部新软化事件，尚无证据归因于本构或接触代码错误。 |
| 时间步、历史和动量 | 每个实际物理步恰提交一次，无κ/d回退、重复失效或负损伤耗能。最大速度**0.792842052m/s**，最大步位移/dx **0.003964225**、相对端点步位移/dx **0.006704868**；动量余额最大**1.51198×10⁻²⁵N·s**。这些项目在实际窗口内通过，不抵消Δd越界。 |
| 独立完整力学核算 | 所有真实步核对完整端点变分模型，材料力最大差**4.2352×10⁻²¹N**、储能最大差**2.2127×10⁻²³J**、每步损伤耗能最大差**1.6750×10⁻²³J**；κ/d差均为0。原始键ID、固定m和停止时存活端点对一致。 |
| 实际动态闭合接触 | 第930步、**232.5µs**，节点247、278间已失效原始键11121对应材料点首次进入1.0dx接触。1576个实际步有接触，同时最多76对，停止时45对；已失效原键拉力为0、没有存活键与pair力重复计力。独立接触力/能量核算一致，但最大软罚力重叠为**0.246649dx（约12.3325µm）**，不等于已认证硬壁不可穿透或生物接触精度。 |
| 实际拓扑 | 停止时4834条原键完全失效、7213条仍存活，依然是**一个341点连通体**；没有可辨识完整分離区域或已测新裂面。停止时最大自由位移109.3576µm、最大自由速度0.5627503m/s，不能称准静态平衡。 |
| 旧0.5µs对照 | 只复用旧0–122µs存档，245个相同物理时刻直接比较，无插值、无旧工况重跑。最大坐标分量差**2.123462nm**、最大速度分量差**0.000149696m/s**；共同末态520个完全失效ID集合一致。新时程到122µs的最大Δd为0.05164741，余额−3.601696779pJ，旧值−7.261265927pJ。局部改善不能认证后续全程或连续时间极限。 |

完整账本仍采用 $W=\\Delta U+\\Delta K+D+\\Delta U_c-J_c+R$，外功由真实边界路径梯形反力功加精确夹持动能变化，D为**同一新坐标**的完整U_old−U_new；−D只作带符号拓扑项展示，绝不再次累计。640.5µs停止现场：**W=4.589016074nJ，U=0.422481930nJ，K=0.889773093nJ，D=3.200673638nJ，Uc=0.098860217nJ，Jc=0**，R=**−22.772804649pJ**。使用与旧工况完全相同的累计尺度，最大相对余额 **0.372426823%**，低于冻结2%门槛；这是已观察窗口的能量通过，完整700µs仍未通过。升载保持期间仍有51条新失效，卸载已观察段又有538条，表明边界停住或回退不意味着自由点已经停止。

累计时间层项为：Verlet功求积缺陷+0.020985435pJ、固定损伤空间求积−0.009172868pJ、损伤力切换功+22.775050094pJ、接触求积缺陷−0.014058012pJ；其负和重建R，最大独立余量约2.23×10⁻²²J。损伤力切换仍为主要数值余额来源，没有被隐藏为材料耗散。自然新裂面面积、材料Gc为 **UNAVAILABLE**，宏观峰值强度及严格时间收敛均未认证；历史 **G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²**、旧运动学173.156623/230.776817nJ及所有前轮记录继续保留。

计算采用1 MPI/单线程，主模拟及逐步在线核查约 **27.864s**，含最终压缩保存约29.703s，进程峰值约660.1MiB。仅限定核查上一轮动态入口和既有manifest，没有重扫31001个历史文件；安全删除0。真实位置、速度、力逐步以FP64保存一次，稀疏损伤变化与关键完整状态均可精确重建。实际生成六组PNG/SVG及 [本轮图表索引](temp_storage/clot_github_sync_20261011/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)，并保存 [36项中文报告](temp_storage/clot_github_sync_20261011/clot_clearance/verification/{TASK}/{RUN}/REVIEW_LPS_NATIVE_SOFTENING_FULL_PROTOCOL_TIMESTEP_ADMISSION_V1.md)、越界前后状态及SHA。

下一轮唯一优先项是**同一冻结341点工况的局部损伤事件分步与力—功—能量一致时间积分研究**：当前数据支持局部事件时间分辨率不足，尚不支持更换材料。先独立设计和预注册算法；不是自动把全局dt降到0.125µs，也不调s_c、s₀、幅值或阻尼。本轮P4缺口试样为NOT_RUN，未部署生产、未运行完整V3、CFD、微泡或AI，交付后停止等待审查。

'''
 anchor='# ⑥ 外加剪切';new=old.replace(anchor,block+anchor,1);after=new.encode();assert new.replace(block,'',1)==old
 for val in ['0.07170070833333','0.06969623734234','173.156623','230.776817']:assert new.count(val)>=old.count(val)
 for sym in '①②③④⑤⑥⑦⑧⑨':assert old.count(sym)==new.count(sym)
 (O/'clot_rupture.before.md').write_bytes(before)
 assert sha(doc.read_bytes())==checkpoint,'Concurrent WSL document modification; do not overwrite'
 doc.write_bytes(after);assert doc.read_bytes()==after
 (O/'clot_rupture.after.md').write_bytes(after)
 (O/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(doc)+' (before)',tofile=str(doc)+' (after)')))
 result=dict(status='PASS_ACTUAL_WSL_FILE_UPDATED',actual_path=str(doc),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),before_sha256=checkpoint,after_sha256=sha(after),before_bytes=len(before),after_bytes=len(after),concurrent_hash_checked=True,readback_identical=True,old_text_byte_for_byte_preserved_after_removing_single_insert=True,original_sections_1_to_9_preserved=True,backup='clot_rupture.before.md',after='clot_rupture.after.md',diff='clot_rupture.diff')
 (O/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(result,indent=2,ensure_ascii=False));print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
