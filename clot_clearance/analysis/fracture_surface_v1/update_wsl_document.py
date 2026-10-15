"""One guarded, insertion-only WSL document update after scientific acceptance."""
from pathlib import Path
import json,hashlib,difflib,os,stat,re
from datetime import datetime,timezone

C=json.loads(Path(__file__).with_name('config.json').read_text()); R=Path(C['local']); O=R/'verification'/C['task']/C['run']; F=R/'results'/C['task']/C['run']; doc=Path('/home/lzy/projects/clot_rupture.md')
def sha(b):return hashlib.sha256(b).hexdigest()
def write(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2))
assert json.loads((O/'IMPLEMENTATION_CHECKS.json').read_text())['status']=='PASS'
assert json.loads((O/'HISTORICAL_PROTECTION_AFTER.json').read_text())['status']=='PASS'
dest=O/'wsl_document'; dest.mkdir(exist_ok=True)
assert not (dest/'clot_rupture.before.md').exists(), 'Preserve existing backup; update is one-shot'
before=doc.read_bytes(); original=before.decode('utf-8'); before_sha=sha(before); before_stat=doc.stat()
assert before_sha=='bef3f226f641cad5f9b14e6d3813d80d647069efd7e424fb8fad1f7195e116cf', 'Latest document changed since full manual reading; reread before proceeding'
anchor='# ⑥ 外加剪切：模拟微泡对血块表面的“刮”'
assert original.count(anchor)==1
addition=r'''**本轮三维裂面标称能量识别（20261012T215617Z）**：实际读取 Zhang & Qiao（2018）完整原文，核对其三维 ordinary state-based 公式 Eq.(2)–(4)、(9)–(12)、(17)–(22)、(26)–(27) 与 Fig.2。对论文影响函数 $\omega=\delta/r$，$\beta=\delta/5$、$\beta'=26\delta/405$，有

$G_{0,\mathrm{theory}}=[(9K-15G)\beta'+15G\beta]s_{00}^{2}=\delta s_{00}^{2}\big[(26/45)K+(55/27)G\big]$。

保持实际输入 $K=33333.333333333336$ Pa、$G=940$ Pa、$s_{00}=0.15$、$\delta=150.5$ µm、$\Delta x=50$ µm、$\alpha=0$，结果为：

| 口径 | K 体积项（J/m²） | G 剪切项（J/m²） | 标称裂面能量（J/m²） |
|---|---:|---:|---:|
| 理想连续三维理论 | 0.065216666667 | 0.006484041667 | **0.071700708333** |
| 真实点阵、z=0 中央窗口 | 0.063433272566 | 0.006262964776 | **0.069696237342** |

理论 K 项约占90.9568%，离散 K 项约占91.0139%；展开式中的 $K-5G/3$ 是混合系数，不能整项称为体积贡献。实际原生影响函数是 $1/r$，同时缩放参考 $m$ 和全部积分分子后，与论文相差的常数 $\delta$ 抵消。原生半单元线性 taper 同时进入完整 $m$、跨面 $m$、膨胀量和偏伸长积分，没有切换成全权重邻域。

使用此前保存、与历史 manifest SHA 一致的 **16000点、909392条完整未切参考键**，没有以 A0/A1/A2 存活键代替。周期 Y 跨缝键131196条，ID、最小镜像、原始体积及固定参考 $m$ 均核对。先冻结 X 中央1.5 mm、覆盖完整 Y 周期0.4 mm的窗口，投影面积 **0.6 mm²**；保留窗口节点的全部邻居。全试样跨 z=0 的唯一键为55408条；窗口4800个节点中，1440个有非零贡献，两侧各720个，均具完整122邻居。计算 $\lambda_i=m_{\mathrm{cross},i}/m_i$，按两侧节点的 $W_iV_i$ 求和后除以一次投影面积，不再乘2。非零窗口节点的 $\lambda$ 约0.00664311–0.39280596。

直接 $\lambda$ 表达式与逐邻居计算膨胀量、体积能和偏能的表达式，总能量相对差 **3.16×10⁻¹⁶**；偏能包含那些形式扩展为零的非跨面邻居。该恒等只认证本轮完整参考邻域。主结果/理论为 **0.9720439165**，低约2.7956%；有限壳层采样、原生 taper 和三层单侧距离求积均被保留，没有调整材料或面积制造吻合。唯一一次 z=50 µm 相邻间隙对照的差为2.78×10⁻¹⁷ J/m²：不同原键ID集合在平移后等价，是同网格相位的内部位置检查，不能称为相位独立或网格收敛。当前 $\delta/\Delta x=3.01$ 的限制仍在。

**0.069696237342 J/m² 是理想 point-state 标称裂面能量，不是已经测得的材料固有 $G_c$。** 形式上令跨面相互作用同时达到 $s_{00}r_0$，未必存在全局可实现位移；本轮没有真实裂纹传播、确定的新生裂面和完整分离功账本。它也不是历史 PVCE 的 $G_{\mathrm{drive}}$ 或逐键独立弹簧耗能。与 Liu 2021 人体全血5.90±1.18 J/m²、Liu 2024 牛血栓最低速率总撕裂能18.27±1.45 J/m²及模型反演内禀功约6–8 J/m²只能比较尺度；材料、速率和耗能定义不同。当前仅支持 `NOMINAL_ENERGY_SCALE_BELOW_EXPERIMENTAL_REFERENCES`，不构成实验误差认证，也不足以要求立即增加渐进软化。

本轮两项离线积分合计约0.0681 s、进程峰值约108.8 MiB，未新建LAMMPS实例，未改生产C++、参数或显式theta顺序。**下一轮唯一优先项调整为受控 Mode-I 分离试验**：在明确的自然扩展面积上核对外功、弹性储能、动能、接触、数值阻尼和带符号拓扑项，判断实际等效扩展阻力与本轮标称值的关系；本轮不自动运行。

> 理论来源：[Zhang & Qiao 2018, International Journal of Fracture 211:217–235](https://doi.org/10.1007/s10704-018-0285-8)。完整计算、文献限制和28项审查见文末本轮链接。

'''
tail=f'\n\n本轮 [三维裂面标称能量识别（28项回答）]({O}/REVIEW_LPS_3D_FRACTURE_SURFACE_ENERGY_IDENTIFICATION_V1.md) 与 [四组真实点阵、局部能量、理论对照和方法审查图]({F}/OPEN_RESULTS.html) 保存两项冻结积分、原始参考键、全精度结果及身份检查。上述各轮历史结论保持；本轮之后的唯一优先建议是受控 Mode-I 分离以测量实际扩展阻力，本轮仅冻结标称能量，停止等待后续科研审查。\n'
after_text=original.replace(anchor,addition+anchor)+tail; after=after_text.encode('utf-8')
assert after_text.replace(addition,'',1)[:-len(tail)]==original
assert re.findall(r'^# .+',original,re.M)==re.findall(r'^# .+',after_text,re.M)
(dest/'clot_rupture.before.md').write_bytes(before);(dest/'clot_rupture.after.md').write_bytes(after)
diff=''.join(difflib.unified_diff(original.splitlines(True),after_text.splitlines(True),fromfile='clot_rupture.before.md',tofile='clot_rupture.after.md'))
assert not [l for l in diff.splitlines() if l.startswith('-') and not l.startswith('---')]
(dest/'clot_rupture.diff').write_text(diff)
assert sha(doc.read_bytes())==before_sha and doc.stat().st_mtime_ns==before_stat.st_mtime_ns,'Concurrent document change detected'
temporary=doc.with_name('.clot_rupture.'+C['run']+'.tmp')
with temporary.open('xb') as fp:fp.write(after);fp.flush();os.fsync(fp.fileno())
os.chmod(temporary,stat.S_IMODE(before_stat.st_mode))
assert sha(doc.read_bytes())==before_sha,'Concurrent change before atomic replacement'
os.replace(temporary,doc);readback=doc.read_bytes();assert readback==after
audit={'status':'PASS_ACTUALLY_UPDATED_AND_READ_BACK','document':str(doc),'updated_utc':datetime.now(timezone.utc).isoformat(),'before_sha256':before_sha,'after_sha256':sha(readback),'before_bytes':len(before),'after_bytes':len(after),'concurrency_checks':2,'atomic_replace':True,'original_text_preserved_exactly':True,'original_headings_and_numbering_unchanged':True,'existing_tables_unchanged':True,'inserted_under_section':'⑤; latest references appended','removed_original_lines':0,'added_lines':sum(l.startswith('+') and not l.startswith('+++') for l in diff.splitlines()),'backup':str(dest/'clot_rupture.before.md'),'after':str(dest/'clot_rupture.after.md'),'diff':str(dest/'clot_rupture.diff'),'readback_sha256':sha(readback)}
write(dest/'DOCUMENT_UPDATE_AUDIT.json',audit);write(O/'WSL_DOCUMENT_UPDATE.json',audit);print(json.dumps(audit,ensure_ascii=False,indent=2))
