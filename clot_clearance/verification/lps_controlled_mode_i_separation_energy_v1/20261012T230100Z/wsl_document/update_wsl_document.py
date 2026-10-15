from pathlib import Path
import json,hashlib,difflib,os,datetime
c=json.loads(Path(__file__).with_name('config.json').read_text());root=Path(c['local']);o=root/'verification'/c['task']/c['run'];p=Path('/home/lzy/projects/clot_rupture.md');before=p.read_bytes();s=before.decode();sh=lambda b:hashlib.sha256(b).hexdigest()
assert sh(before)=='ba8b5875e71be49b06ddfc2b1783d4733a262d08d4244ce241c9fae156785780','Concurrent document change: read latest before applying'
assert json.loads((o/'IMPLEMENTATION_CHECKS.json').read_text())['status']=='PASS_IMPLEMENTATION_AND_FAILURE_REPORTING_ONLY'
marker='# ⑥ 外加剪切：模拟微泡对血块表面的“刮”';assert s.count(marker)==1 and '20261012T230100Z' not in s
insert=r'''**本轮受控 Mode-I 分离能准入（20261012T230100Z）**：复用 X×Y×Z=5×0.4×1 mm 的16000点周期层、909392条完整参考键，初始 z=0 缺口长度1.50 mm、前缘 X=−1 mm；材料、点距、作用半径及最终双修复原生库不变。先在完整拓扑建立参考 $m$，再去除明确的16416条初始预切键，剩余892976条。原生逐ID $m$ 相对差约2.19×10⁻¹⁵；零载初态原生力、膨胀量、储能、动能和接触能为零，初始态可信。

针对历史端面统一轴向位移的夹具竞争，本轮只选择一次**上下外表面均匀法向牵引**候选：上下各800点，每点±50 nN，等效±20 Pa、单面合力40 µN；横向及轴向位移均由受力自由响应，Y保持周期。没有内部裂尖强制运动、特殊强度保护或额外删键。一次独立无损伤静态预检的端点，裂尖/加载端最大 $s/s_{00}$ 为 **0.924699/0.824168**，伸长排序有所改善；不过，这不是已观察到的自然起裂，也不能认证夹具失效已经排除。

**该候选未通过运行前冻结的预检门槛，正式断裂没有启动。** 最大伸长 **0.138705** 超过本轮预检的小变形上限0.03；最大残余力 **1.954×10⁻¹² N** 超过1×10⁻¹³ N。优化器因目标下降达浮点条件而结束，不能把其 success 当作物理平衡验收。零载初态通过与已加载端点不通过是两个不同结论；未放宽门槛、换载荷或追加夹具试错。

保存端点平均开度约 **140.609525 µm**，独立可回收 LPS 储能 **2.508694009 nJ**；原生力与独立储能梯度相对差约1.73×10⁻¹³，保留完整固定 $m$ 的力相容展开式通过核对。终点 $\sum_i\mathbf F_i^{\rm ext}\cdot\mathbf u_i=5.624381003$ nJ 只是当前外载势能收缩，**不是实际从零加载的路径外功**。优化迭代没有物理时间，两个原生静态坐标核对没有积分器，不能重建动能、阻尼、接触功、lagged-theta功或断裂能余额时程。

55408条原始跨面键中，16416条为初始预切、**38992条仍存活**；静态检查没有新拓扑变化，材料仍是一个16000点连通体。面积识别保留全部跨Y周期及非局部长键，候选每个50 µm推进列对应 $B\Delta a=2\times10^{-8}$ m²，初始6×10⁻⁷ m²预切面积排除。交点与完整X投影两种检查在初始前缘有非局部影响带差异，不能用“断键数×点距平方”或局部损伤比例代替真正新生裂面。正式动态新增面积、真实首断时间和能量账本均为未测量。

因此 **G_SEPARATION_PROTOCOL 与 GC_EFFECTIVE 均为 UNAVAILABLE，而不是零**。上一轮 **G0_DISCRETE = 0.06969623734234 J/m²** 的标称结果完整保留；本轮没有实际阻力可与它比较，不能判定相容或不相容，也不足以支持“材料太脆”或必须升级软化。直接阻碍是所选牵引候选的预检幅值与静态平衡准入；没有发生可用于归因的正式断裂、非目标损伤或能量失守。

本轮实际为1个独立预检（1132次函数评估，约24.90 s）和1个原生静态核对实例（约6.04 s），**正式原生断裂工况0、动力接受步0**。提供5组审查图和1组明确未接受端点的补充图，不制造断裂演化帧或物理时程。**下一轮唯一优先项**是对同一固定分布牵引试样开展独立的无损伤静态准入：预先选择真正低载的小变形检查与可靠残差恢复方法，核查裂尖/加载端响应及最慢变形时间尺度；不迭代夹具、不改材料，通过后再单独开展自然断裂。本轮到此停止。

'''
links=f'\n本轮 [受控Mode-I分离能审查（34项回答）]({o}/REVIEW_LPS_CONTROLLED_MODE_I_SEPARATION_ENERGY_V1.md) 与 [五组审查图及未接受预检补充图]({root}/results/{c["task"]}/{c["run"]}/OPEN_RESULTS.html) 保存单次预检的冻结门槛、原生静态核对、真实存活拓扑及UNAVAILABLE结论。新增物理结论限于静态准入失败；不改动上一轮标称值及历史计算。\n'
after_s=s.replace(marker,insert+marker,1)+links;after=after_s.encode();sm=difflib.SequenceMatcher(a=s,b=after_s,autojunk=False);assert all(op in ['equal','insert'] for op,*_ in sm.get_opcodes());assert after_s.replace(insert,'',1).removesuffix(links)==s
out=o/'wsl_document';out.mkdir(exist_ok=True);assert not (out/'clot_rupture.before.md').exists()
(out/'clot_rupture.before.md').write_bytes(before);(out/'clot_rupture.after.md').write_bytes(after);diff=''.join(difflib.unified_diff(s.splitlines(True),after_s.splitlines(True),fromfile='clot_rupture.before.md',tofile='clot_rupture.after.md'));(out/'clot_rupture.diff').write_text(diff)
assert p.read_bytes()==before,'Concurrent change before temp write';tmp=p.with_name('.clot_rupture.mode_i_20261012T230100Z.tmp');fd=os.open(tmp,os.O_WRONLY|os.O_CREAT|os.O_EXCL,p.stat().st_mode&0o777)
with os.fdopen(fd,'wb') as f:f.write(after);f.flush();os.fsync(f.fileno())
assert p.read_bytes()==before,'Concurrent change before atomic replace';os.replace(tmp,p);back=p.read_bytes();assert back==after
audit={'status':'PASS_ACTUALLY_WRITTEN_AND_READ_BACK','path':str(p),'before_sha256':sh(before),'after_sha256':sh(after),'readback_sha256':sh(back),'before_bytes':len(before),'after_bytes':len(after),'old_text_unchanged':True,'only_insertions':True,'all_1_to_9_headings_preserved':True,'concurrent_checks':2,'atomic_replace':True,'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':{name:sh((out/name).read_bytes()) for name in ['clot_rupture.before.md','clot_rupture.after.md','clot_rupture.diff']}}
(out/'DOCUMENT_UPDATE_AUDIT.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2));(o/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2));print(json.dumps(audit,ensure_ascii=False,indent=2))
