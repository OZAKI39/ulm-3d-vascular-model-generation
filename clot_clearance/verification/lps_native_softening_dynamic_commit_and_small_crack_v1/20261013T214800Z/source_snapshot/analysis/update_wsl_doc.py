from pathlib import Path
import json,hashlib,difflib
root=Path('/home/lzy/projects/temp_storage/clot_github_sync_20261011/clot_clearance');task='lps_native_softening_dynamic_commit_and_small_crack_v1';run='20261013T214800Z';o=root/'verification'/task/run;doc=Path('/home/lzy/projects/clot_rupture.md');before=(o/'clot_rupture.before.md').read_bytes();meta=json.loads((o/'WSL_DOCUMENT_BEFORE.json').read_text());assert hashlib.sha256(before).hexdigest()==meta['sha256'];assert doc.read_bytes()==before,'Concurrent document change: re-read and merge required'
text=before.decode();marker='# ⑥ 外加剪切：模拟微泡对血块表面的“刮”';assert text.count(marker)==1
addition='''**本轮真实动力学接入与停止审查（20261013T214800Z）**：先执行P0，实际审计28个候选源码文件，确认可安全删除的废弃代码为0，故未删除文件。现行数学模型、复现代码、历史source_snapshot和用途不确定的编译探针均保留；新建精简的 `ACTIVE_DEVELOPMENT_ENTRYPOINTS.md`。开发入口为服务器 `analysis/dynamic_softening_v1/`，在独立源码中更新 `bond_bpm_peri_variational.h/.cpp` 并新增内部fix `fix_lps_variational_dynamic.h/.cpp`。冻结端点邻域自由能、原参考m/r/taper、G、K、密度及s0不变，主动态仅用sc=0.18。实际计算共享库SHA256为 **`866065394c48ebc8a6d32f147bde992a325c2a503a49d50392243eee10fea9bd`**；生产库及原两套诊断库保持不变。

真实调用采用LAMMPS `Verlet::run → fix nve`。自由点先用旧力完成第一次速度半步和位置推进；新协调fix在 `POST_INTEGRATE`、邻居重建之前，以新坐标提交一次最大历史伸长κ和不可逆d，保存固定坐标下完整自由能降幅。新完全失效键在本步删除owned拓扑及1–2排斥屏蔽；同一步pair/bond计算使用新状态。在 `POST_FORCE`、第二次速度半步前，把原BPM力密度乘V换成N，质量保持真实ρV。`run 0`与试探坐标求力不提交损伤，试探坐标的储能允许改变；恢复接受坐标后力及全部历史一致。这是损伤算子分裂，**不能沿用纯弹性Verlet的二阶精度宣称整个断裂过程二阶收敛**。

主工况复用341个真实点和12047条原始键，左右端面分别49、47点规定x向平滑运动，245点的三个分量由原生力自由积分。预注册每侧幅值60µm、0.3ms升载、0.05ms保持、0.3ms卸载及0.05ms末保持，总0.7ms；dt=0.5µs，无阻尼。完整初始切线的线性Verlet临界dt约13.9074µs，只用于初始尺度筛查，不替代损伤和接触安全验证。

| 真实观察或核查 | 本轮结果及边界 |
|---|---|
| 首次自然软化 | 第170步，**85µs**，原始bond ID **3393**。 |
| 首次完全失效 | 第185步，**92.5µs**；同一步7条原键：**2180、2194、3368、3393、3904、9197、10264**。均由实际伸长触发，不人为删键制造事件。 |
| 主时程停止 | 第244步、**122µs**，最大单步Δd=**0.1025542483**，超过冻结0.10门槛；立即保存并停止。实际失效520条、存活11527条，仍为一个341点连通体。没有完成预注册700µs时程或卸载段。 |
| 自由节点及积分 | 自由点最大位移**23.3556µm**；全点最大速度**0.444338m/s**；最大步间位移/dx为**0.00443339**。原生两次半步方程最大位置误差约1.36×10⁻²⁰m、速度误差5.47×10⁻¹⁷m/s；动量余额最大**1.09×10⁻²⁵N·s**。 |
| 力、储能与历史 | 全245个保存态独立重建，节点力最大差**3.3881×10⁻²¹N**，完整储能最大差**1.0960×10⁻²³J**；κ/d逐值相同、无回退，每个物理步只提交一次。原生完整存活端点对与原始ID表核对一致，固定m差为零。 |
| 失效后传力与接触 | 完全失效原键的非接触拉力为零；当前自然张开窗口没有触发闭合接触，接触能为零。**不能把上一轮49对静态闭合接触记作本轮动态闭合已验证**。 |
| 唯一dt/2窗口 | 从同一第210步接受restart出发，仅检查105–122µs，dt=0.25µs，68步。最大Δd降至**0.05165050**；末态520条失效ID集合与粗步一致，最大坐标分量差约**1.226nm**。窗口能量余额由−6.05133pJ降至−3.00260pJ，比值**0.49619**。这是局部改善，不是完整dt/2时程或全局收敛证明；没有dt/4。 |

按实际离散时间层建立 $W=\Delta U+\Delta K+D+\Delta U_c-J_c+R$：边界功由真实夹持点位移与梯形约束反力求积，并加入规定边界动能的精确变化；它是数值功估计。停止时外功 **5.013269603nJ**，完整可回收储能 **2.039926979nJ**，实际动能 **2.518907725nJ**，不可逆耗散 **0.461696165nJ**，数值阻尼、接触储能及接触拓扑注能均为零。带符号材料拓扑变化为−D，单独显示但不再次计入右侧。总余额 **−7.261265927pJ**，按冻结累计能量尺度最大相对误差 **0.144631%**，没有把未知残差塞入D。

独立时间层分解得到：Verlet功求积缺陷+0.0462424pJ、固定损伤空间求积−0.0118855pJ、损伤后力切换功+7.2269090pJ；它们的负和重建上述余额，差约2.82×10⁻²³J。主要误差来源可解释，仍存在有限时间离散误差，**不称严格能量守恒或全程步长认证**。当前能量尺度不能与历史全节点规定运动学的173.156623nJ直接作材料韧性优劣比较。

两项开发中断均保留：标量min/max表达式不被该LAMMPS接受，修正为等价ternary前为0个接受步；第170步只读测试误要求不同试探坐标的U不变，修正Python断言后，从同一保存接受态继续，位置、速度、κ/d、D及提交计数逐值一致，未重跑或修改物理工况。其后第244步的损伤增量越界是本轮真正的科学停止原因。编译及增量生命周期修正合计约40.101s，主时程与短dt/2窗口合计约3.863s（含观测和保存），动力进程峰值约657.6MiB；未进行新MPI动态扩展测试。

**P4为NOT_RUN**：主时程未通过冻结损伤时间步门槛，短窗口结果不足以准入新的缺口试样。没有新建缺口加载工况，也没有可辨识完整裂面；在仅作拓扑诊断的x=20µm切面仍有2041条存活跨面键。新增裂面面积、材料Gc为 **UNAVAILABLE**，宏观强度和生产伸缩性未认证。新、旧二值失效模型历史力向量L2差 **8.814180703270×10⁻⁶N**继续保留；在本次停止坐标及相同二值掩码下，两种自由能给出的力L2差为1.239512741×10⁻⁵N，这是本构重分配差异，没有经验修正强迫吻合。

本轮在已记录窗口内证实了真实原生积分、不可逆自然损伤、完全失效及力学/动量/历史的一致性，能量余额已明确量化；**整体动力学准入仍未通过，不允许据此直接进入更大非均匀裂纹或宏观标定**。历史 **G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²**及上一轮173.156623/230.776817nJ均保持。下一轮唯一优先项是：对**同一341点冻结工况完成损伤时间步准入与全时程能量验证**，先预注册完整时程的更小步长，不改材料和夹具；通过后再审查P4。本轮到此停止，不自动续算、部署生产、运行V3、CFD、微泡或AI。

'''
new=text.replace(marker,addition+marker)
new+='\n\n本轮 [真实原生软化动力学审查（51项回答）]('+str(o/'REVIEW_LPS_NATIVE_SOFTENING_DYNAMIC_COMMIT_AND_SMALL_CRACK_V1.md')+') 与 [五组真实动力学、损伤及停止证据图]('+str(root/'results'/task/run/'OPEN_RESULTS.html')+') 保存P0零删除依据、隔离源码/库、244步主时程、唯一短dt/2窗口、明确停止状态和能量分解。P4未运行，生产库及历史科研证据保护核对通过。\n'
for token in ['0.07170070833333','0.06969623734234','143.8594','173.1566','230.7768']:assert token in text and token in new
assert doc.read_bytes()==before;doc.write_text(new);after=doc.read_bytes();assert after==new.encode();(o/'clot_rupture.after.md').write_bytes(after);diff=''.join(difflib.unified_diff(text.splitlines(True),new.splitlines(True),fromfile=str(doc)+'.before',tofile=str(doc)));(o/'clot_rupture.diff').write_text(diff)
result=dict(status='PASS',actual_WSL_path=str(doc),before_SHA256=hashlib.sha256(before).hexdigest(),after_SHA256=hashlib.sha256(after).hexdigest(),before_bytes=len(before),after_bytes=len(after),concurrent_change_check='UNCHANGED_BEFORE_WRITE',readback_exact=True,old_text_preserved=new.replace(addition,'').startswith(text),historical_parameters_and_energies_preserved=True,backup=str(o/'clot_rupture.before.md'),after=str(o/'clot_rupture.after.md'),diff=str(o/'clot_rupture.diff'))
(o/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(result,indent=2,ensure_ascii=False));print(json.dumps(result,ensure_ascii=False))
