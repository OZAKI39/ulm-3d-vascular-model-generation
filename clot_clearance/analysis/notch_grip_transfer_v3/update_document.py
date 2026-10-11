"""Minimal guarded update of the actual WSL scientific notebook."""
from pathlib import Path
import hashlib,json,difflib,os,datetime
C=json.loads(Path('/home/lzy/projects/.notch_grip_v3_current.json').read_text());base=Path(C['local']);p=Path('/home/lzy/projects/clot_rupture.md');prov=base/'provenance';sha=lambda b:hashlib.sha256(b).hexdigest();checks=json.loads((base/'IMPLEMENTATION_CHECKS.json').read_text());status=json.loads((base/'FINAL_STATUS.json').read_text());assert checks['status']=='PASS' and status['MATERIAL_GC_STATUS']=='GC_UNAVAILABLE'
old=p.read_bytes();start=json.loads((prov/'WSL_START.json').read_text());assert sha(old)==start['document_sha256'],'Notebook changed since startup: merge latest content before writing';assert (prov/'clot_rupture.before.md').read_bytes()==old
text=old.decode();anchor='**本轮缺口对照**：';assert text.count(anchor)==1;text=text.replace(anchor,'**前轮缺口对照（V2）**：',1)
boundary='\n---\n\n# ⑥ 外加剪切：';assert text.count(boundary)==1
paragraph='''
**本轮夹持传力对照（V3）**：只对最外侧两个端面的 288 个点规定轴向运动，其余 2592 个点自由，全部横向分量按真实受力积分；夹具参考距离从 850 µm 变为 950 µm，仍用 680 µm、24 ms 的原均值开度函数。3204 条双端轴向受控键同处端面，没有运动学必然拉断，实际也没有断开；但首断在 5.2555 ms 的夹具角点，**未保留缺口优先起裂**。峰值 0.29398 kPa，目标带断键占 67.89%、夹具区占 13.33%。到 7.05 ms 按运行前冻结的“单一目标带解释丢失”规则停止：总断键 5598 条，夹具区占 20.72%（旧 M5B 同时刻为 48.21%），非跨原裂面键占 77.60%。中央已有损伤，但 12 个厚度条带均未形成从缺口出发的连续跨面完全失连带，仍为一个 2880 点连通体；**新增裂面面积和表观断裂功/面积都不可认证，GC_UNAVAILABLE 保持**。独立能量余额约为最大外功的 0.0117%，断裂后动能/储能升至 12.23%；负拓扑增量保留，数值阻尼不当作生物黏性。仅运行一个正式 M5C 和一个零步检查，没有重试或材料参数调整；这是实施验证通过的负面科研结果，不是完成材料标定，也不足以支持软化升级。
'''
text=text.replace(boundary,'\n'+paragraph+boundary,1)
prior='下一轮唯一优先项是设计不会必然拉断受规定夹持层、同时保留缺口起裂优势的夹持传力边界，再核查新增单裂面；本轮没有自动开始后续计算。';assert text.count(prior)==1
new=f'本轮 [夹持传力与缺口断裂审查]({base}/REVIEW_LAMMPS_NOTCH_FRACTURE_GRIP_TRANSFER_V3.md) 和 [V3 图与数据索引]({base}/OPEN_RESULTS.html) 记录了端面方案及其负面结果。下一轮唯一优先项是基于已有 M5/M5B/M5C 与 Liu 原文实验，开展只读的断裂观测量可辨识性与本构需求审查，明确哪些实验量可约束 LPS 和软化需求；结束同类夹具试错。本轮到此停止，没有自动开始后续计算。'
text=text.replace(prior,new,1);updated=text.encode();diff=''.join(difflib.unified_diff(old.decode().splitlines(True),text.splitlines(True),fromfile=str(p)+' (before)',tofile=str(p)+' (after)'))
assert sha(p.read_bytes())==sha(old),'Concurrent modification detected immediately before write'
temp=p.with_name('.clot_rupture.v3_'+C['run_id']+'.tmp');assert not temp.exists();temp.write_bytes(updated);temp.chmod(p.stat().st_mode);os.replace(temp,p);assert p.read_bytes()==updated
(prov/'clot_rupture.after.md').write_bytes(updated);(prov/'clot_rupture.changes.diff').write_text(diff)
receipt=dict(state='ACTUALLY_UPDATED_AND_READ_BACK',actual_WSL_path=str(p),before_sha256=sha(old),after_sha256=sha(updated),backup=str(prov/'clot_rupture.before.md'),diff=str(prov/'clot_rupture.changes.diff'),after_copy=str(prov/'clot_rupture.after.md'),concurrent_hash_check='Startup hash matched immediately before atomic replacement; after content read back',scope='Relabel previous V2 notch paragraph; add V3 negative scientific result in section⑤; replace only final next-priority sentence and add V3 links.',utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(prov/'WSL_DOCUMENT_UPDATE.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2));print(json.dumps(receipt,ensure_ascii=False,indent=2))
