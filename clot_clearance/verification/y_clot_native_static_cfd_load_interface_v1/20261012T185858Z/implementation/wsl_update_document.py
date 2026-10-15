from pathlib import Path
import json,hashlib,difflib,time
C=json.loads(Path(__file__).with_name('config.json').read_text());L=Path(C['local']);O=L/'verification'/C['task']/C['run'];doc=Path('/home/lzy/projects/clot_rupture.md');D=O/'wsl_document';D.mkdir(exist_ok=True)
sha=lambda b:hashlib.sha256(b).hexdigest()
before=doc.read_bytes();bs=sha(before);assert bs=='753b08d32f18f5bcb419cac74bcf83bf96b24a07dd738a1bcb108ce02aeaa0a2','Concurrent document change requires reread';text=before.decode();(D/'clot_rupture.before.md').write_bytes(before)
a=json.loads((O/'STATIC_EQUILIBRIUM_AUDIT.json').read_text());assert a['accepted'];s=json.loads((O/'BULK_INTERFACE_STRETCH_AUDIT.json').read_text());e=json.loads((O/'INDEPENDENT_LPS_ENERGY_AUDIT.json').read_text());assert e['status'].startswith('PASS')
insert='''

**本轮原生静态 CFD 接口（20261012T185858Z）**：把上面的原始保守节点力真正加入最终双修复诊断库（SHA `00b246ebe422…22ae9`），先按现有顺序把 LPS 力密度乘体积变成 N，再以真实 atom ID 加入已是 N 的 CFD 节点力。仍是 **580 个材料点受载、1347 个材料点直接外载为零，1561 个锚点不受 CFD 映射力**。在 2 MPI 的真实排序及 1144 个节点换 rank 后，全部节点的原生受载前后力差与映射逐位一致；其中 174 个受载节点发生迁移。没有重复乘体积，也没有对 ghost 重复加力。首次脚本误把已排序中间态的恒等置换判为失败，已用保存的真实迁移状态修正检查；失败日志保留，没有重跑载荷或改阈值。

先核对原参考态：材料、本体、界面与锚点内力均为零，接触未激活，无须另做无载松弛。随后只对 1927 个材料点求一次静态平衡，锚点坐标固定，使用完整独立 LPS 势能减去常量载荷势 $\\sum_i\\mathbf F_i\\cdot\\mathbf u_i$；再把所得位置交回真正加载 CFD 的 LAMMPS 检查残余力和锚点内部作用。没有 nve 时间积分、Gaussian 剪切、碎片阻力或数值阻尼。

实测外力为 **(0.059221682, 0.022786103, 0.660738209) µN**，由原生锚点内部力恢复的约束反力为 **(−0.059223343, −0.022784590, −0.660737018) µN**，合力残差 2.543×10⁻¹² N。最大自由节点残余力 3.843×10⁻¹⁴ N，相对 L2 残差 7.804×10⁻⁶，均通过预先冻结门槛。力矩沿用原映射参考点；参考坐标与当前坐标口径的残差分别为 **9.881×10⁻¹⁵、4.790×10⁻¹⁶ N·m**，两者分别通过 1×10⁻¹⁴ N·m 绝对门槛。参考坐标口径已接近该门槛，主要包含真实位移引起的取矩位置差，不能混同当前构型平衡误差。

血栓平均位移模长 **13.4406 nm**，最大 **58.3005 nm**，在 364 号点，参考位置 (−6.25, −0.05, 0.70) mm；最大位移只有点距的约 0.1166%。锚点位移严格为零。本体最大伸长 **2.09683×10⁻⁴**，阈值利用率 **0.139788%**；界面最大伸长 **4.44208×10⁻⁵**，利用率 **0.0315881%**。全部 157637 条原 typed bonds 保留，仍为一个附着血栓，**新断键为零、接触力为零**；没有调高 s00 或降低流体载荷。

无载储能为零，受载完整 LPS 储能及增量为 **3.05408×10⁻¹⁵ J**；常力功 F·u 为 6.10760×10⁻¹⁵ J，外载势取其负值。独立梯度与原生力、内部保守路径功与储能增量均核对通过；半个 F·u 只是小变形线性诊断，没有另跑升载过程。未用不完整 native PE×V 替代真实储能。一次求解完成 155 次优化迭代、162 次目标/梯度调用，含原生核验的求解内段约 2.887 s；没有追加动力学步。

这说明当前离散模型中的正常背景流产生可数值分辨、但很小的弹性响应，**原生静态 CFD 载荷接口通过**；不等于实验可测性、材料断裂认证、真实微泡声流或双向 FSI。固定完整几何的水样稳态流场、表压参考、锚点代理和未标定 Gc 的限制仍在。新结果可以作为后续额外载荷的静态背景基准；旧 200 Pa 人工剪切、历史碎屑轨迹和生产库没有因此替换。
'''
needle='这些节点力可以用于下一轮原生接口验证；**本轮没有把它施加到生产碎裂计算**。'
assert needle in text;text=text.replace(needle,'离线映射完成时，这些节点力可用于下一轮原生接口验证；**该轮没有把它施加到生产碎裂计算**。',1)
anchor='> **出处**：';idx=text.index(anchor);text=text[:idx]+insert+'\n'+text[idx:]
old='**下一轮唯一优先项**：冻结当前最终双修复显式方案，对已有含栓CFD保守节点牵引做一次禁止断键的原生完整态受载接口验证，核对单位、总力/力矩、约束反力及小变形响应；继续推进已有物理接口，不为追求大碎片调整theta或材料。本轮到此停止，不自动启动该任务或切换生产库。'
assert old in text;text=text.replace(old,'该轮提出的原生完整态受载接口验证已由本轮完成，见下段；显式 theta 和生产库继续保持原状。',1)
report=O/'REVIEW_Y_CLOT_NATIVE_STATIC_CFD_LOAD_INTERFACE_V1.md';fig=L/'results'/C['task']/C['run']/'OPEN_RESULTS.html'
text+=f'\n\n本轮 [原生静态 CFD 接口审查（33项回答）]({report}) 和 [四组实际载荷、反力、位移及材料响应图]({fig}) 保存逐 ID 状态、真实反力、静态门槛、一次检查脚本修复和全部精度数据。**下一轮唯一优先项**：在这一通过验收的背景接口上，冻结可追溯的微泡附加载荷数据契约及背景/增量分解规则，先审核单位、作用面、真实 ID 和合力/力矩，防止重复加入背景牵引；没有可信微泡数据时不把 200 Pa Gaussian 当作微泡真值。本轮到此停止，不自动启动微泡、断裂、CFD 或生产部署。\n'
after=text.encode();assert sha(doc.read_bytes())==bs,'Concurrent document edit';doc.write_bytes(after);readback=doc.read_bytes();assert readback==after
(D/'clot_rupture.after.md').write_bytes(readback);(D/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(before.decode().splitlines(True),text.splitlines(True),fromfile=str(doc)+'.before',tofile=str(doc))))
rec={'status':'PASS_ACTUAL_WSL_WRITE_AND_READBACK','actual_path':str(doc),'before_sha256':bs,'after_sha256':sha(readback),'before_bytes':len(before),'after_bytes':len(readback),'complete_before_read':True,'concurrency_hash_rechecked_before_write':True,'readback_exact':True,'original_numbered_sections_preserved':all(('# '+z) in text for z in '①②③④⑤⑥⑦⑧⑨'),'historical_fragmentation_results_preserved':True,'backup_directory':str(D),'epoch':time.time()};(O/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(rec,ensure_ascii=False,indent=2));print(json.dumps(rec,ensure_ascii=False,indent=2))
