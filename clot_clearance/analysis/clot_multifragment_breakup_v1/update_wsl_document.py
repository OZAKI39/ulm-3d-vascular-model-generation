"""WSL-only, history-preserving research log update after validation."""
from pathlib import Path
import json,hashlib,difflib
root=Path(__file__).resolve().parents[2];task='clot_multifragment_breakup_v1';run='20261017T110849Z';out=root/'verification'/task/run
doc=Path('/home/lzy/projects/clot_rupture.md');folder=out/'wsl_document';folder.mkdir(exist_ok=True)
def read(p):return json.loads(p.read_text())
def sha(raw):return hashlib.sha256(raw).hexdigest()
assert read(out/'DELIVERY_VALIDATION.json')['status']=='PASS'
cases=read(out/'ANALYSIS.json');resources=read(out/'RUN_RESOURCES.json');afterbridge=read(out/'CRACK_BRIDGE_AFTER_LOADING.json')['cases'];status=read(out/'FINAL_STATUS.json')
raw=doc.read_bytes();text=raw.decode();heading=f'**多区域时空加载与真实多块碎裂（{run}，2026-10-17）**'
assert heading not in text
backup=folder/'clot_rupture.before.md'
if not backup.exists():backup.write_bytes(raw)
if backup.read_bytes()!=raw:
    backup=folder/'clot_rupture.immediate_before.md';assert not backup.exists();backup.write_bytes(raw)
parts=[heading+'：本轮正式以MULTI-FRAGMENT BREAKUP为目标，区分表面侵蚀、单大块脱离和多个中等内部凝聚块体。沿用1792材料点、336固定锚点、95666原始键、2736原预切和均匀LPS变分软化，不新增缺陷、不人为删除桥接键。',
    '探索性目标在正式仿真前冻结：共同2.8ms群组及末态≥5个20–200节点真实自由分量，脱离质量≥初始10%，最大块≤脱离质量50%，1–2节点细屑≤脱离质量20%，至少5个群组原始成员保持连通至3.6ms。独立保存各项约束，不用最大块或加权总分代替目标。六个尺寸组为1–2、3–13、14–33、34–113、114–200、≥201节点，均非医学粒径分类。',
    '历史离线诊断：均匀末态上半部896点诊断区有3461条真实跨区存活键，3423条位于切口尖端右侧，2987条长度≥2dx，2105条当前受压，3441条历史κ未越过s0；存在避开|Y|≤50µm中心观察层的真实锚点路径。表面帽/切口唇尚有2961/3439条桥接。非均匀对应3467/2959/3384条。失效簇出生与合并不等于材料分离；有界固定源核最小割只是完整度诊断代理，未作为删键或断裂能。',
    '新原生加载为MULTISITE_MECHANICAL_SURROGATE：两个外部作用点(−325,−75,175)µm和(325,75,−175)µm，分别沿+X/−X压入；R=182.818µm、Kind=57973.6N/m²、间隙25µm、压入240µm。每源0–0.2ms接近、1.2ms到最大、保持至1.4ms、2.4ms撤回；同步同时启动，错时第二源延后0.4ms，延迟按半个点间剪切波传播时间取整。仅改变时序，两组均至3.6ms；实际反力和外功并不相等。无摩擦径向接触不是切向摩擦牵引，也不是已验证微泡/射流。']
for k,label in [('simultaneous','同步'),('staggered','错时')]:
    q=cases[k];v=q['final'];m=q['final_metrics'];bridges=afterbridge[k]['final']['regions'];fine=m['fine_debris_mass_fraction_detached'];dominance=m['largest_free_fragment_fraction_detached']
    parts.append(f"{label}：首次软化/完整失效{q['first_softening_s']*1e3:.6g}/{q['first_failure_s']*1e3:.6g}ms；末态脱离{int(v['detached_nodes'])}点，历史/末态最大块{q['history_max_free_nodes']}/{int(v['largest_nodes'])}点，历史最多/末态20–200节点块{q['history_max_N20_200']}/{m['fragment_count_target_band']}个；末态≥34/≥114点块{m['N_ge34']}/{m['N_ge114']}个。细屑/脱离质量={fine}，最大块/脱离质量={dominance}，总脱离/初始质量={m['total_detached_mass_fraction']}。实际输入功{v['W_J']*1e9:.8g}nJ，两源峰值反力{[f*1e6 for f in q['peak_site_reactions_N']]}µN，最大重叠/dx={q['max_site_overlap_dx']}；最大相对能量残差{q['max_relative_energy_residual']:.7g}。冻结目标通过={q['exploratory_target_met']}。三诊断区残余桥接={[r['surviving_crossing_bonds'] for r in bridges]}条。")
    parts.append(f"{label}共同卸载群组：2.8ms目标带块{len(q['common_cohort_medium_ids'])}个，原成员连通至末态{q['common_cohort_medium_surviving']}个；共同窗口后二次裂解{q['secondary_after_common_cohort']}次，末尾实际连续无源接触{q['actual_contact_free_tail_s']*1e6:.6g}µs。撤回完成不等于实际不再接触；无中等块时不能声称其稳定。完整原始成员、内部存活率、体积/质量等效尺寸、三维跨度、回转半径、PCA主轴、形态秩及内部存活键密度逐接受态保存，不造平滑壳面。")
parts += [status['SCIENTIFIC_CONCLUSION'],
    '从连接主体形成多个中等块，需要多个有实际体积区域的完整分离边界：既切断所有跨界bulk/interface存活路径（含长horizon键和三维绕行），又保留块内连接。仍未达阈值的桥接同时受加载覆盖与正伸长损伤规则影响；本轮不能唯一归因于本构错误。后续可把有文献/标定依据的混合模式损伤作为独立假设验证，不能以更容易得到碎块作为正确性标准。',
    '原生修改限于两个命名接触源的物理单位力累加；材料核及历史身份未改。新接口无接触/释放、方向、势能梯度、反力/外功分项及叠加检查通过，正式VV和不可逆历史逐步核验，原始ID损伤事件重放与全部接受态分量回读通过。能量W0+W1=ΔU+ΔK+D+ΔUc−Jc+ΔUi+R，无阻尼；R不当作耗散。一次固定材料球心穿越预检的长步积分失败和修正fixture完整保留；正式协议未改。接触罚刚度和离散分辨率尚未做收敛认证。',
    f"两组累计墙钟{resources['native_wall_s']:.3f}s、CPU{resources['native_cpu_s']:.3f}s，峰值RSS{resources['peak_RSS_bytes']/1024**3:.6f}GiB；1rank×1线程顺序运行。只运行两组，无追加对照、完整V3、CFD、材料扫描。微泡Streaming/微射流到实际接触面积、方向、时间谱、冲量和外功的物理映射仍缺失，结果不足以直接认证真实V3多块机制。",
    f"见[可视化总览](temp_storage/clot_github_sync_20261015/clot_clearance/results/{task}/{run}/OPEN_RESULTS.html)、[18项科研问答](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{task}/{run}/REVIEW_CLOT_MULTIFRAGMENT_BREAKUP_V1.md)、[真实多块指标](temp_storage/clot_github_sync_20261015/clot_clearance/verification/{task}/{run}/MULTIFRAGMENT_RESULTS.csv)。本轮结束后停止，等待审查。"]
block='\n\n'.join(parts)+'\n\n';index=text.index('# ⑥ 外加剪切');new=text[:index]+block+text[index:]
assert doc.read_bytes()==raw;doc.write_text(new);after=doc.read_bytes();assert after.decode()==new and new.replace(block,'',1)==text
(folder/'clot_rupture.after.md').write_bytes(after)
(folder/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),new.splitlines(True),fromfile='clot_rupture.before.md',tofile='clot_rupture.after.md')))
(out/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(dict(status='PASS',path=str(doc),before_backup=str(backup),before_sha256=sha(raw),after_sha256=sha(after),history_preserved_exactly=True,readback_equal=True),ensure_ascii=False,indent=2)+'\n')
print('WSL_DOCUMENT_UPDATED',sha(after))
