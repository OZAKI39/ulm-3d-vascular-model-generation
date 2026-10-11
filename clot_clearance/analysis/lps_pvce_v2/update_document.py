"""Apply supported new conclusions to the actual WSL document, with concurrency guard."""
from pathlib import Path
import json,hashlib,difflib,os,tempfile,time
C=json.loads(Path('/home/lzy/projects/.pvce_v2_current.json').read_text());B=Path(C['local']);p=Path('/home/lzy/projects/clot_rupture.md');before=p.read_bytes();h=hashlib.sha256(before).hexdigest();initial=hashlib.sha256((B/'provenance/clot_rupture.before.md').read_bytes()).hexdigest();assert h==initial,'Concurrent edit: reread and reconcile before writing';text=before.decode();assert '本轮连续体与近起裂审查（PVCE V2）' not in text
addition='''**本轮连续体与近起裂审查（PVCE V2）**：新增六个三维小应变线弹性混合有限元状态，裂纹两侧使用真正独立的自由面，材料 G、K 与 LPS 相同。保留 600×600×1000 µm 名义实体；实际轴向控制面在 z=±475 µm，间距 950 µm，两外侧各保留 25 µm 弹性端帽；只用三个横向点自由度消除刚体零模。这个映射没有混淆节点中心与体积外边界，但连续控制面与非局部节点的半单元传力仍不完全等价。

2 µm 开度下，局部加密混合 FEM 的两点/三点能量释放率分别为 **6.0930588×10⁻⁶、6.0104634×10⁻⁶ J/m²**；LPS 相对它分别高约 **1.5246%、2.6594%**。一次裂尖带加密的 G 变化为 0.1561%、0.2217%，通过预注册 5% 差分稳定性筛查；不过，同一位移场的未投影体积应变能重建比混合变分能仍高 26.86%–27.39%，**完整连续体场与厚度方向尚未认证收敛**。这些数值只作为有限网格参考，不能把 2.66% 差异全部归为非局部效应，也不能用它宣称 LPS 已获物理验证。

近起裂搜索先重读 V1 的真实 A0 最终状态：最大伸长是 **0.00407827**，此前约 0.00535 是优化试探历史最大值。由真实最终值预测的 73.5606 µm 仅作为首个加载候选。实际候选最大伸长 **0.143025**、夹具利用率 **0.953499**，缺口侧表面/中央为 0.832050/0.694233；四个最大键位于参考中点 (275,±275,±450) µm 的端面角点，顶部反力约 72.0964 µN。**该候选未通过静态验收，反力不是临界拉力，也不能说夹具已经首先断裂。**

候选原生相对自由力残差 6.19×10⁻⁷ 虽通过，最大绝对自由力仍为 **1.1291×10⁻¹² N**，超过运行前冻结的 1×10⁻¹³ N 门槛。原生力与独立梯度相对差约 9.35×10⁻¹⁶、接触为零；L-BFGS-B 因目标下降达到浮点分辨率而停止。这是数值验收失败，不能解释成已证明材料失稳。按冻结停止规则不放宽门槛、不追加候选；**真实首次阈值开度区间与位置均未确定，高载 A1/A2 未运行，G_high 不可用而不是零**。最后严格有效的 A0 仍是 V1 的 2 µm 平衡态。

本轮实际新增 6 个 FEM 求解和 1 个 LPS 候选。低载 G_drive 不等于临界撕裂阻力；未验收 A0 的储能比平方外推高约 5.69%，也不能代替高载三态 G。材料固有 Gc 仍未识别；不足以证明“LPS 比真实血块脆多少倍”、解释 Y 血栓细碎屑的唯一原因，或支持修改 s00 / 引入新耗能机制。生产材料卡、原生 LAMMPS 与旧 V1–V3/PVCE V1 结果保持不变。

'''
marker='---\n\n# ⑥ 外加剪切';assert text.count(marker)==1;text=text.replace(marker,addition+marker)
old='下一轮唯一优先项是为同一微型试样、裂纹和位移边界建立独立线弹性连续体 G_drive 参考，评估非局部及有限域偏差；不自动启动。V1–V3 的自然裂纹面积和 Gc 未识别结论保持。本轮到此停止。';assert text.count(old)==1
new=f'本轮 [连续体与近起裂审查]({B}/REVIEW_LPS_PVCE_CONTINUUM_AND_NEAR_ONSET_AUDIT_V2.md) 和 [四组 V2 审查图]({B}/OPEN_RESULTS.html) 记录了有限网格参考及高载静态门槛失败。下一轮唯一优先项是保持材料、原生力和验收门槛不变，修复并独立验证 A0 高载静态求解器的浮点停机问题，达标后再恢复首达阈值区间搜索；不自动启动。V1–V3 的自然裂纹面积和 Gc 未识别结论保持。本轮到此停止。'
text=text.replace(old,new);backup=B/'provenance/clot_rupture.before_write.md';backup.write_bytes(before);diff=''.join(difflib.unified_diff(before.decode().splitlines(True),text.splitlines(True),fromfile=str(p)+'.before',tofile=str(p)));(B/'provenance/clot_rupture.diff').write_text(diff)
assert p.read_bytes()==before;mode=p.stat().st_mode;fd,tmp=tempfile.mkstemp(prefix='.clot_rupture_pvce_v2_',dir=p.parent)
with os.fdopen(fd,'wb') as fp:fp.write(text.encode());fp.flush();os.fsync(fp.fileno())
os.chmod(tmp,mode);assert p.read_bytes()==before;os.replace(tmp,p);after=p.read_bytes();assert after==text.encode()
r=dict(status='ACTUAL_WSL_UPDATED_AND_READ_BACK',file=str(p),before_sha256=h,after_sha256=hashlib.sha256(after).hexdigest(),same_as_start_backup=True,backup=str(backup),diff=str(B/'provenance/clot_rupture.diff'),changed_only='Section⑤ new PVCE V2 evidence and final review/one nextpriority paragraph; all historical paragraphs retained.',utc=time.strftime('%FT%TZ',time.gmtime()));(B/'provenance/DOCUMENT_UPDATE.json').write_text(json.dumps(r,indent=2));s=json.loads((B/'FINAL_STATUS.json').read_text());s['DOCUMENT_UPDATE']=r['status'];(B/'FINAL_STATUS.json').write_text(json.dumps(s,indent=2));print(json.dumps(r,indent=2))
