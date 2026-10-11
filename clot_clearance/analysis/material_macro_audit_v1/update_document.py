"""Task-specific conditional WSL edit, with immediate hash check and text backup."""
import datetime,difflib,hashlib,json,os,stat,tempfile
from pathlib import Path

C=json.loads(Path('/home/lzy/projects/.macro_audit_current.json').read_text())
b=Path(C['local']);p=Path('/home/lzy/projects/clot_rupture.md')
checks=json.loads((b/'IMPLEMENTATION_CHECKS.json').read_text())
assert all(checks[n]['independent_gradient'] and checks[n]['energy_residual_below5percent'] for n in ['M1','M2','M3','M3P','M4','M5','M6','T1','D1'])
raw=p.read_bytes();before_hash=hashlib.sha256(raw).hexdigest()
assert before_hash==json.loads((b/'provenance/local_start.json').read_text())['document_sha256'],'Concurrent edit: reread and merge before writing'
s=raw.decode();before=s
def replace(old,new):
 global s
 assert s.count(old)==1,old[:100]
 s=s.replace(old,new)

replace('这两个数交给 LPS 模型（之前讲过：K 管“胀”，G 管“歪”）。',
 '这两个数交给 LPS 模型（之前讲过：K 管“胀”，G 管“歪”）。它们来自牛全血凝块的未损伤小应变映射，不等于完整 Yeoh 或损伤模型。\n\n**本轮检查**：内部点自由、只约束加载边界的小试样测得 $E=2.681$ kPa、$G=0.759$ kPa、$\\nu=0.4831$。E 比理论映射低约 4.0%，G 低约 19.2%，后者没有通过预设的 15% 筛查。原侧边固定试样的中心仍有横向应变，不能直接拿来算 K；补做固定横向周期盒的压缩，得到 $K=33.55$ kPa（用测得 G 换算；用输入 G 则为 33.31 kPa）。这些是有限试样结果，生产参数没有因此改动。')
replace('**门槛 s₀₀ 怎么定**：做一个小方块（8×8×8 个点）的数值拉伸实验，量出它最大能承受多大应力；调 s₀₀，直到这个强度落进 Liu 2024 的实验范围。',
 '**门槛 s₀₀ 怎么定**：旧标定在 8×8×8 点小方块中给全部节点规定仿射运动，得到约 6.78 kPa；旧界面约 1.96 kPa 来自规定刚体拉脱。这些历史结果保留，但不能当成内部自由试样已验证的材料强度。目前仍用本体 $s_{00}=0.15$、界面 $s_{00}\\approx0.140624859375$；界面参考是猪皮代理，不是人体血管内壁粘附。')
start=s.index('反过来算，$s_{00}=0.15$');end=s.index('\n\n---',start)
replace(s[start:end],
 '这个公式属于键型模型，不能直接给当前 state-based LPS 指定断裂能；**0.2 J/m² 不是本模型已知的断裂能**。本轮重读 Liu 2024 正文和原始补充材料：18.27±1.45 J/m² 是最低加载率下的实验总断裂能；约 6–8 J/m² 是特定模型与强度条件下反演的 intrinsic work，两者不能混为同一个范围。\n\n**本轮实际拉伸**：600×600×1000 µm、点距 50 µm、作用半径 150.5 µm，内部点自由。无缺口、缺口和三层锚点拉脱的名义峰值分别为 **0.488、0.303、0.475 kPa**，都已越过峰值；前两者首断发生在夹具过渡角点，拉脱是本体先断、混合损伤，所以不能直接称为血块固有强度或纯界面粘附强度。加载采用 24 ms 平滑时程，不是论文慢速实验复现。\n\n断前/断后拓扑的储能变化已独立计算，能量残差在预设门槛内；但部分变化为负，缺口试样约 87% 的新断键不跨预设裂面，单裂纹新增面积无法可靠识别，因此 **Gc 仍未识别**。时间步 0.5→0.25 µs 的峰值差约 0.003%、整条应力曲线差约 4.78%；点距 50→40 µm 时曲线差约 44.2%、净拓扑能量差约 48.1%，不能宣称网格收敛。现在没有足够证据写成“强度对了但非常脆”，也不能据此认定必须加软化。')
replace('意思是让力作用得够久，整个血块都能“感受到”它，而不是只有表面抖一下。',
 '意思是让力作用得够久，整个血块都能“感受到”它，而不是只有表面抖一下。这个 200 Pa 脉冲仍是工程测试输入，尚未由微泡或 PZT 功率标定。\n\n**微泡补充核查**：平行项目当前 F7 路线采用 Euler 预测和 12 次接触投影。对一个可追溯的 32 泡成熟群集、0–1 s 工况，仅将实际步长从 10 µs 减到 5 µs（保存间隔是 500 µs）。粗粒度群集比例接近、接受穿透低于既有 0.01R 几何门槛，但同 ID 粒子位置最大差约 6.85R，内部接触轨迹未建立时间收敛。近场物理也未标定，不能据此认证群集 AI 标签或本节载荷。')
replace(r'$F=15\,k\,\bar V\,\frac{d-r}{\delta}\quad(r<d),\qquad k=\frac{18K}{\pi\delta^4},\qquad d=1.0\,\Delta x$',
 r'$F_i=15\,k\,V_i\bar V\,\frac{d-r}{\delta}\quad(r<d),\qquad k=\frac{18K}{\pi\delta^4},\qquad d=1.0\,\Delta x$'+'\n\n这里写的是 N：原生数组先给力密度，现有程序乘本点体积 $V_i$ 后才加外载、壁力和阻尼。本轮修正的是文字中漏写的体积因子，生产程序的转换顺序已核对正确。')
old=next(t for t in s.split('\n\n') if t.startswith(r'$\Delta t=0.5\ \mu$s，共'))
replace(old,r'$\Delta t=0.5\ \mu$s，共 4,660 步，即 2\.33 ms。单弹簧振子的 Verlet 稳定条件是 $\omega\Delta t<2$；当前壁弹簧 $\omega\Delta t=0.316$，材料纵波速度约 5.88 m/s，$c_L\Delta t/\Delta x=0.0588$。这些是有用的尺度筛查，不能单独证明整个非局部网络与断裂接触的稳定性。'+'\n\n**能量核对**：原生 PE 即使乘体积，仍未完整计入 LPS 偏能；不能只看变量名 U_J 就拿来算断裂能。本轮增加与锁定源码的力相容的只读势能评估，并独立记录外功、自由点动能、接触储能、数值阻尼、带符号拓扑能量和残差。没有改动生产受力与历史轨迹；依赖旧 PE 的能量解释需要复审。')
replace('质心落在流场网格的某个四面体里，就用四个顶点的速度按距离加权平均：如果质心落在网格外面，不外推，这一步不给力。',
 '质心落在流场网格的某个四面体里，就用四个顶点的速度按四面体重心坐标作线性加权（P1 插值）：如果质心落在网格外面，不外推，这一步不给力。')
replace('新版碎块的速度和血流基本一致，代价是计算时间变成约 3 倍，因为每一步都要调用一次 Python。',
 '新版末态碎块速度接近血流，但这不代表整个过程中都很慢。当前全点力记录有 2,104 个保存时刻，重建完整步速度后，已保存节点最大值为 **3.91997 m/s**（0.3795 ms，节点 364）；脱离节点的 4,068 个逐步记录中，游离组件质心最大值为 **3.72918 m/s**。没有每步保存的附着节点峰值仍未知，不能用末态中位数代替。4,661 个原生壁日志时刻中无活动罚接触，最小中心净距约 25.80 µm。计算时间变成约 3 倍的历史记录保留，原因是每一步都要调用一次 Python。')
replace('- **定理**：这是线性常微分方程的精确解，所以就算步长比 τ 大很多也不会算炸。代码规定每步不超过 5 ms，且一步位移不超过 20 µm。',
 '- **定理的范围**：这是单步内冻结 $\\mathbf u$ 和 τ 后的精确解，不是空间变流、接触和整段轨迹已收敛的证明。实际 232,276 条组件积分记录中，步长最大 5 ms、中位约 3.463 ms，接受位移不超过 20 µm；包含交接初态的质心速度最大约 0.02994 m/s。')
replace('- **管壁处理**：如果碎块里某个点要钻进管壁 5 µm 以内，就沿法向把它推回来，并去掉朝向管壁的速度（不弹、不摩擦）。',
 '- **管壁处理**：如果碎块里某个点要钻进管壁 5 µm 以内，就沿法向把它推回来，并去掉朝向管壁的速度（不弹、不摩擦）。实际有 1,025 次约束事件，最大位置修正约 0.221 µm；原生段的 20 µm 罚接触已切换成这里的 5 µm 约束。\n\n- **交接范围**：质心位置、速度、同态拖曳律和绝对时钟已核对连续；内部成员偏移与附着主体随后冻结，内部速度最多截断约 0.117 m/s。这不等于内部转动、变形也连续。')
report=b/'REVIEW_LAMMPS_MATERIAL_MACRO_CALIBRATION_AND_TIMESTEP_AUDIT_V1.md'
replace('所以结论是 **附着主体 \\+ 表面碎屑脱落**，还没有出现大块断裂。',
 '所以结论是 **附着主体 \\+ 表面碎屑脱落**，还没有出现大块断裂。\n\n本轮 [材料与时间步审查]('+str(report)+') 保留了实际曲线、点阵状态和适用限制；数值检查通过不等于材料匹配或生物实验验证。下一轮唯一优先项是重审自由内部试样与旧仿射小单元的宏观强度映射及边界影响，再决定是否调整本构。')
assert s!=before
backup=b/'provenance/clot_rupture.immediate_before.md';assert not backup.exists(),'This task edit has already been attempted'
backup.write_bytes(raw)
diff=''.join(difflib.unified_diff(before.splitlines(True),s.splitlines(True),fromfile=str(backup),tofile=str(p)))
(b/'provenance/clot_rupture.diff').write_text(diff)
assert hashlib.sha256(p.read_bytes()).hexdigest()==before_hash,'File changed before write; backup and proposal retained'
mode=stat.S_IMODE(p.stat().st_mode)
with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=p.parent,prefix='.clot_rupture_update_',delete=False) as f:
 f.write(s);f.flush();os.fsync(f.fileno());tmp=Path(f.name)
os.chmod(tmp,mode)
assert hashlib.sha256(p.read_bytes()).hexdigest()==before_hash,'Concurrent edit just before replacement'
os.replace(tmp,p)
after_hash=hashlib.sha256(p.read_bytes()).hexdigest();assert p.read_text()==s
(b/'provenance/clot_rupture.after.md').write_text(s)
receipt=dict(status='UPDATED_ACTUAL_WSL_FILE',path=str(p),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),before_sha256=before_hash,after_sha256=after_hash,immediate_backup=str(backup),diff=str(b/'provenance/clot_rupture.diff'),concurrent_change_detected=False,condition='Independent implementation/energy evidence checked; reliable mismatches and identifiability limits retained',scope='Existing sections IV-IX and final review link; original headings and historic observations preserved')
(b/'DOCUMENT_UPDATE.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False))
f=b/'FINAL_STATUS.json';d=json.loads(f.read_text());d['document_update']=receipt;f.write_text(json.dumps(d,indent=2,ensure_ascii=False))
print(json.dumps(receipt,indent=2,ensure_ascii=False))
