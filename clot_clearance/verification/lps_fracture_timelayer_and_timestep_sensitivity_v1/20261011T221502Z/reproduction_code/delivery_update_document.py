from pathlib import Path
import json,hashlib,difflib,re
C=json.loads(Path(__file__).with_name('config.json').read_text());O=Path(C['local'])/'verification'/C['task']/C['run'];p=Path('/home/lzy/projects/clot_rupture.md');before=p.read_text();oldsha=hashlib.sha256(p.read_bytes()).hexdigest()
assert oldsha=='76f5de1f75e518345c6d6a49044b4b108e360085ff5dc7a2558f6e7e33ba7408','Concurrent document change: reread and adapt, never overwrite'
dest=O/'document_provenance';dest.mkdir(exist_ok=False);(dest/'clot_rupture.before.md').write_text(before)
paragraph='''**本轮断键时间层与时间步审查（20261011T221502Z）**：先确认历史341点、12047键的“五条断键”是直接改坐标的静态诊断，没有实际Velocity-Verlet速度或连续加载时程。复用其同坐标状态，正常拓扑力变化最大8.673 µN，断前theta造成的力差最大4.411 µN，带符号储能变化为−59.133 pJ；这些量不能全部当成额外数值冲量。

经明确允许后，只建立一套动力学适配：参考点阵、原子ID、键、G、K、s00及15 µm局部位移幅值不变，从历史T0预应变态及零速度出发；90号点在10 µs内按半余弦沿x移动15 µm，随后保持到50 µs，其余340点自由积分。原小试样没有接触和阻尼，本轮也不新增；这不是旧静态事件的原配置动力学复现，也不是准静态平衡初态。分别完成0.5 µs的100步和0.25 µs的200步，每步直接读取原生半步与整步速度。

两组首先断开的都是原键202（4–90），实际首断为5.00/4.75 µs；末态都是同五条断键，仍为一个341点连通体，没有出现额外级联。首断时全点最大时间层力差为0.5933/0.5486 µN，自由点最大值为0.3662/0.3396 µN；最大全点差位于被规定运动的90号点，不能直接换成该点的自由速度增量。首断自由点的完整两次半更新冲量代理为1.8312×10⁻¹³/8.4906×10⁻¹⁴ N·s，比值0.4637；全窗口累计最大代理比值0.5067。真实动量另用质量乘整步速度差核对，余额约3.3×10⁻²⁶ N·s。**代理沿实际轨迹计算，不是另一种theta算法的实测速度差。**

独立能量账本保留初始储能、边界功、自由点动能和同坐标带符号拓扑变化；累计拓扑储能变化为−40.2948/−39.2861 pJ。最大余额残差为0.7569/0.3008 pJ，占运行前固定的222.9515 pJ参考尺度0.3395%/0.1349%。实际路径上的theta滞后力净功为−1.6804/−0.8599 pJ，拓扑切换时间项为+2.3095/+1.1089 pJ，其余为常规积分与功求积项；不能把残差小理解为每一项都为零，也不能把这些项加成材料Gc。两组自由点动能峰值仅差约0.0384%，最大自由点速度约0.23187/0.23194 m/s。

本小试样支持有限时间层影响随步长减小，**仍未完成全局时间收敛证明，不能推广为完整V3的稳定步长认证**。保留现有显式lagged-theta方案；没有足够证据要求立即开发两阶段断键，也不能解释V3的1–2点细碎屑主因。生产库、材料参数与此前199点/190游离组件结论保持原样。实际完成两次动力计算，另一次输入表达式解析失败只保存初态、没有接受动力学步；失败和修正记录均保留。

'''
anchor='# ⑥ 外加剪切：模拟微泡对血块表面的“刮”';assert before.count(anchor)==1
after=before.replace(anchor,paragraph+anchor)
old='**下一轮唯一优先项**：冻结最终双修复库，做受控小试样的断键时间层与时间步敏感性审查，区分正常拓扑能量释放与显式theta滞后的数值影响。本轮到此停止，不自动部署、调材料或启动新计算。'
review=O/'REVIEW_LPS_FRACTURE_TIMELAYER_AND_TIMESTEP_SENSITIVITY_V1.md';fig=Path(C['local'])/'results'/C['task']/C['run']/'OPEN_RESULTS.html'
new=f'上述小试样审查已完成，见[时间层与时间步审查（28项回答）]({review})和[四组实际审查图]({fig})。当前只支持这套局部瞬态的有限工程筛查，不宣称完整V3时间收敛。**下一轮唯一优先项**：冻结当前最终双修复显式方案，对已有含栓CFD保守节点牵引做一次禁止断键的原生完整态受载接口验证，核对单位、总力/力矩、约束反力及小变形响应；继续推进已有物理接口，不为追求大碎片调整theta或材料。本轮到此停止，不自动启动该任务或切换生产库。'
assert after.count(old)==1;after=after.replace(old,new)
assert re.findall(r'^# [①-⑨].*$',before,re.M)==re.findall(r'^# [①-⑨].*$',after,re.M)
assert hashlib.sha256(p.read_bytes()).hexdigest()==oldsha
p.write_text(after);readback=p.read_text();assert readback==after
(dest/'clot_rupture.after.md').write_text(readback);(dest/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile=str(p)+'.before',tofile=str(p))))
audit={'status':'PASS_ACTUAL_WSL_UPDATED','actual_path':str(p),'before_sha256':oldsha,'after_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'read_back_equal':True,'concurrent_write_guard':'SHA checked immediately before write','original_numbered_titles_retained':True,'historical_199_190_and_native_bugfix_retained':True,'full_V3_convergence_claimed':False,'backup_directory':str(dest),'scope':'New verified small-specimen results and limitations; historical conclusions retained'}
(O/'DOCUMENT_UPDATE_AUDIT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2));print(json.dumps(audit,ensure_ascii=False))
