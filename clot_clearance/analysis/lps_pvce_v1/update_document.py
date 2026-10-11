from pathlib import Path
import json,hashlib,difflib,os,tempfile,time
C=json.loads(Path('/home/lzy/projects/.pvce_current.json').read_text());B=Path(C['local']);p=Path('/home/lzy/projects/clot_rupture.md');before=p.read_bytes();sha=hashlib.sha256(before).hexdigest();start=hashlib.sha256((B/'provenance/clot_rupture.before.md').read_bytes()).hexdigest();assert sha==start,'Concurrent document edit: reread and reconcile before writing';text=before.decode();assert '本轮虚拟裂纹能量释放（PVCE）' not in text
addition=r'''
**本轮虚拟裂纹能量释放（PVCE）**：按 Zhang & Qiao（2020）的固定位移能量差思想，单独计算 a=200、250、300 µm 的三个预定裂纹。固定端面总开度 2 µm，材料、点距和作用半径保持不变；复用独立 LPS 储能及其梯度，通过静态 L-BFGS 求平衡，并与原生实际力核对。三个状态均通过运行前冻结的自由力、绝对力、储能平台、位移稳定和反力平衡条件；原生自由力相对残差为 7.81×10⁻⁷、9.04×10⁻⁷、8.91×10⁻⁷，没有额外自然断键，最终接触能与接触力均为零。静态优化不产生动力学 K/U 曲线。

三个裂纹累计切除 2536、3302、4068 条真实参考键，后两次各增加 766 条。统一采用 x<裂尖的半开几何规则，修正了旧 M5 中 34 条前缘浮点交点的选择；三个状态都固定使用完整未切参考邻域的加权体积 m，膨胀量按存活键重算，因此本轮 A0 不等同旧 M5 的“预删后重新归一化”状态。这两项定义在求解前已冻结，旧 V1–V3 结果不改。

用完整三维独立储能 $U_0=1.8055120142$、$U_1=1.6199333620$、$U_2=1.4334157154$ pJ，取单个裂纹推进面积 $B\Delta a=600\times50\ \mu\mathrm m^2=3.0\times10^{-8}\ \mathrm m^2$，不乘双面系数，也不用断键数量换算面积。固定位移下 $G_{\mathrm{two}}=-(U_1-U_0)/(B\Delta a)=6.1859551\times10^{-6}$ J/m²；三点式 $G_{\mathrm{three}}=-(-3U_0+4U_1-U_2)/(2B\Delta a)=6.1703052\times10^{-6}$ J/m²。相对三点值的差异约 0.2536%。储能差约为经验求解误差尺度的 6500–6900 倍；三点 G 的传播数值分辨尺度约 ±1.82×10⁻⁹ J/m²，但这不是严格误差界，未包含虚拟增量、非局部边界、网格和有限尺寸误差。

**当前新增能力是试样相关的 $G_{\mathrm{drive}}$，不是材料固有 $G_c$**。只求了三个平衡态，三态只给一个局部能量导数估计；未验证临界自然扩展、网格收敛或完整 G(a) 曲线。原生力与独立能量梯度的相对差约 10⁻¹⁵；不完整 native PE、累计外功、动能及数值阻尼均未进入 G 的能量分子。Liu 的实验总断裂能不作为这个微型低开度试样的调参目标。现有新证据仍不足以认定必须更换失效机制或增加软化。

'''.replace(chr(92)*2,chr(92))
marker='---\n\n# ⑥ 外加剪切';assert text.count(marker)==1;text=text.replace(marker,addition+marker)
old='下一轮唯一优先项是基于已有 M5/M5B/M5C 与 Liu 原文实验，开展只读的断裂观测量可辨识性与本构需求审查，明确哪些实验量可约束 LPS 和软化需求；结束同类夹具试错。本轮到此停止，没有自动开始后续计算。';assert text.count(old)==1
new=f'本轮 [PVCE 能量释放可行性审查]({B}/REVIEW_LPS_PVCE_ENERGY_RELEASE_FEASIBILITY_V1.md) 和 [四组真实科研图]({B}/OPEN_RESULTS.html) 记录了三态平衡、能差与适用范围。下一轮唯一优先项是为同一微型试样、裂纹和位移边界建立独立线弹性连续体 G_drive 参考，评估非局部及有限域偏差；不自动启动。V1–V3 的自然裂纹面积和 Gc 未识别结论保持。本轮到此停止。'
text=text.replace(old,new);data=text.encode();backup=B/'provenance/clot_rupture.before_write.md';backup.write_bytes(before);diff=''.join(difflib.unified_diff(before.decode().splitlines(True),text.splitlines(True),fromfile=str(p)+'.before',tofile=str(p)));(B/'provenance/clot_rupture.diff').write_text(diff)
assert p.read_bytes()==before,'Concurrent write detected';mode=p.stat().st_mode
fd,tmp=tempfile.mkstemp(prefix='.clot_rupture_pvce_',dir=p.parent)
with os.fdopen(fd,'wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
os.chmod(tmp,mode);assert p.read_bytes()==before;os.replace(tmp,p);after=p.read_bytes();assert after==data
record=dict(status='ACTUAL_WSL_UPDATED_AND_READ_BACK',file=str(p),before_sha256=sha,after_sha256=hashlib.sha256(after).hexdigest(),same_as_start_backup=True,backup=str(backup),diff=str(B/'provenance/clot_rupture.diff'),changed_only='Section⑤ fracture/energy addition and final review/one nextpriority paragraph',concurrency_check='Exactbyte comparisons before editing and immediately before atomic replace; start hash unchanged.',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));(B/'provenance/DOCUMENT_UPDATE.json').write_text(json.dumps(record,indent=2))
s=json.loads((B/'FINAL_STATUS.json').read_text());s['DOCUMENT_UPDATE']=record['status'];(B/'FINAL_STATUS.json').write_text(json.dumps(s,indent=2));print(json.dumps(record,indent=2))
