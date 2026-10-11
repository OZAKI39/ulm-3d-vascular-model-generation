"""One serial, guarded update of the real WSL research document."""
from pathlib import Path
import hashlib,json,difflib,re,time
C=json.loads(Path('/home/lzy/projects/.traction_failure_current.json').read_text());ROOT=Path(C['root']);out=ROOT/C['common']/'document_edit';doc=Path('/home/lzy/projects/clot_rupture.md')
assert not out.exists(),'Do not overwrite document audit'
before=doc.read_bytes();old=before.decode();h=lambda b:hashlib.sha256(b).hexdigest();out.mkdir();(out/'before.md').write_bytes(before)
s=old.replace('尚未映射到 LAMMPS 材料点，更未改动下面的工程剪切脉冲。','该轮尚未映射到 LAMMPS 材料点；后续离线映射见下文，下面的工程剪切脉冲没有改动。',1)
added_b='''**本轮离线保守映射（20261011T184653Z）**：重新读取上面的加密含栓解，把真实暴露面的压力、切向黏性和法向黏性分别积分，再按邻近非负权重分到现有材料点，最后做最小保力、保力矩修正。1927 个材料点中有 **580 个受载、1347 个外载为零**，1561 个锚点全部排除。有效面积来自实际三角面积分配，不把每点面积直接设为 Δx²，也不按两种体积之比放大力。输出已是 **N**，接入时不能再乘一次材料点体积。

最初有 10 个很小的接触边缘面被切平面条件排除，占暴露面积约 0.0233%；本轮只做一次局部邻接支持修正，随后全部面覆盖。邻接图弦长不是严格曲面测地距离，仍属于离散工程筛查。压力/剪切/总项的 L2 修正分别为 **0.618% / 9.015% / 0.772%**；剪切接近预设 10% 警示尺度。各分项总力、总力矩的归一化代数误差约 10⁻¹⁶，未出现压力或剪切反向节点力；这不表示物理误差也达到该数量级。50 µm 点距、曲面平均和边缘支持仍有限制。

这些节点力可以用于下一轮原生接口验证；**本轮没有把它施加到生产碎裂计算**。旧无栓背景流、断裂轨迹、200 Pa 工程脉冲和碎片输运保持原来的来源。新牵引仅对应固定完整血栓，不是微泡 streaming，也不是双向流固耦合；未来还要避免与已脱离碎片的背景阻力重复计力。

'''
marker='> **出处**：';pos=s.index(marker);s=s[:pos]+added_b+s[pos:]
added_a='''**本轮表面修正与能量失效审查（20261011T184653Z）**：对真实 M5C 缺口试样和 V3 血栓分别计算初始邻域与 Shojaei/perifit 权重，没有重新捏形状或恢复缺口键。M5C 权重约 0.577–3.190，均为正；V3 约 −0.330–49.161，出现 **15 个负权重点、59 条非正平均权重键**，算子残差也没有普遍改善。因此不能把“算出了权重”写成“表面力学已修好”，V3 权重不进入生产计算。

只在正权重 M5C 上做了一次保留原核、初始归一化的变分适配诊断：同一储能的负梯度与节点力通过差分检查，误差小于 10⁻⁸，但同一形变下储能和力明显变大，尚无物理边值真值证明这是改善。它不是 Shojaei 论文逐式复现。Foster 2011 全文也不能直接提供当前 LPS 的独立逐键弹簧能：实际闭合形变路径上，总相互作用功与总储能一致，单条相互作用的累计功却可不回零，因为膨胀量包含多键交叉项。数值阻尼另算，不充当生物耗散。

所以本轮仍**不启用能量断键或渐进软化**。本体 s₀₀=0.15 和原瞬断规则保留为定性工程基线，不宣称已匹配材料强度或韧性；独立 Gc 仍不可用。若以后研究新的耗散机制，还缺同材料、同速率的无缺口强度、稳定裂纹扩展能、指定血管壁界面粘附以及循环/松弛实验。

'''
pos=s.index('# ⑥');boundary=s.rfind('---',0,pos);assert boundary>0;s=s[:boundary]+added_a+s[boundary:]
oldnext='下一轮唯一优先项是验证含栓 CFD 暴露面牵引向现有 LAMMPS 材料点的离线保力、保力矩映射，先核对单位、坐标、面身份和面积权重，不直接启动生产加载或双向耦合；不自动启动。'
assert oldnext in s
s=s.replace(oldnext,'当时建议的离线映射已在本轮完成，仍未进入生产动力学。',1)
s+='\n本轮 [含栓牵引映射与 LPS 失效可行性共同报告]('+str(ROOT/C['common']/'COMBINED_REVIEW.md')+') 和 [七组真实审查图、节点力与分项报告]('+str(ROOT/C['common']/'OPEN_RESULTS.html')+') 保存完整证据。下一轮唯一优先项是：用现有 V3 完整态、禁止断键，做一次原生静态受载接口验证，核对本轮节点力的单位、总载荷、约束反力/力矩及小变形响应，并排除碎片阻力重复作用；不自动运行。\n'
headers=lambda text:re.findall(r'^# .*[①②③④⑤⑥⑦⑧⑨].*$',text,re.M)
assert headers(s)==headers(old) and len(headers(old))==9
for text in ['GC_UNAVAILABLE','未保留缺口优先起裂','没有血栓','Gc 仍未识别','高载']:
 assert text in s
assert doc.read_bytes()==before,'Concurrent modification detected; no write'
doc.write_text(s);readback=doc.read_bytes();assert readback==s.encode()
(out/'after.md').write_bytes(readback);(out/'change.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile=str(doc)+'.before',tofile=str(doc))))
audit=dict(path=str(doc),written_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),sha256_before=h(before),sha256_after=h(readback),read_latest_before_write=True,concurrent_content_check=True,readback_match=True,chapters_changed=['①背景流新增离线映射','⑤新增权重/失效限制','文末证据与下一项'],all_nine_headers_unchanged=True,historical_negative_results_preserved=True,old_unclotted_CFD_and_trajectories_not_relabelled=True)
(out/'AUDIT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n');print(json.dumps(audit,ensure_ascii=False))
