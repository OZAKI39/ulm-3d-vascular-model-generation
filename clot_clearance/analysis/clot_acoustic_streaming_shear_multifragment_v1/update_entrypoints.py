"""Prepend the current flow/shear entrypoint; preserve the previous pointer text."""
from common import *
import difflib
assert json.loads((O/'DELIVERY_VALIDATION.json').read_text())['status']=='PASS'
q=json.loads((O/'ANALYSIS.json').read_text());p=R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md';raw=p.read_bytes();old=raw.decode();assert TASK not in old
backup=O/'entrypoints';backup.mkdir(exist_ok=True);(backup/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md').write_bytes(raw)
lines=[f'# 当前声流／剪切开发入口（{RUN}）',f'当前任务：`{TASK}`。两组3.6ms原生小试样已完成，材料与native库保持不变。',
 '**当前主线为流体剪切，不继续压头迭代。** 真实微泡Streaming尚未确定；可计算路径为明确标注的连续反向Stokes剪切模板及独立机械需求牵引。',
 f'- 入口：`analysis/{TASK}/acoustic_streaming_shear.py`、`streaming_interface.py`、`shear_loading.py`、`native.py`、`run_shear_cases.py`。表面牵引、N单位原始ID传递、不可逆损伤与真实存活连接均实际执行。',
 '- 试样1792材料点+336锚点；dx50µm，原95666键/2736预切保持。实际顶面非负三角权重；不新增预切、不改G/K/s0/sc。',
 '- 模式：sin(2πy/700µm)沿x反向剪切，25µm间隙、μ1mPa·s、4mm/s规定上边界；2ms升、2.4ms降、2.8ms停驱、3.6ms结束。B只放大牵引至1000Pa，不代表微泡流体能力。']
for name,x in q.items():lines.append(f'- {name}：峰值{x["peak_shear_Pa"]:.9g}Pa；末态游离{int(x["final"]["detached_nodes"])}点，历史最大块{x["history_max_free_nodes"]}点，历史最多20–200节点块{x["history_max_N20_200"]}，完整目标={x["exploratory_target_met"]}。')
lines += ['- 声学/气泡输入已核对，但缺少当前血栓配准、间隙及近壁一阶/二阶流场。APECSS缓存黏度2mPa·s，与上述独立模板不同；声压、R(t)、质心速度均不能直接当作Streaming剪切。',
 '- V3 3611实际CFD三角面仅作映射审计。旧带符号映射存在明显局部差异；非负核牺牲力矩。不能仅凭全局守恒认证局部传递。未运行新完整V3、FSI或优化trial。',
 '- 参考面冻结，移动/脱离后仍沿原始ID加载；大变形阶段只支持死牵引机械探索。新裂面应力和当前表面重算尚未完成。',
 f'- 先读 `verification/{TASK}/{RUN}/FINAL_STATUS.json`、`ACTIVE_STREAMING_SHEAR_PIPELINE.json`、20问 `REVIEW_CLOT_ACOUSTIC_STREAMING_SHEAR_MULTIFRAGMENT_V1.md`、`DELIVERY_VALIDATION.json`；图表入口 `results/{TASK}/{RUN}/OPEN_RESULTS.html`。',
 f'- 复用库 `{LIB.relative_to(R)}`，SHA `{sha(LIB)}`。本轮没有更改生产本构、MPI或restart；新实验必须新RUN，保留本轮封存输入与结果。',
 '- 后续先补齐微泡近壁流体来源与局部表面耦合条件，再讨论微泡参数反向优化和V3。不得把损伤前兆、单大块或细粉当成多个凝聚性碎块。',
 '## 历史入口原文（保留，不代表当前加载主线）']
new='\n\n'.join(lines)+'\n\n'+old;assert p.read_bytes()==raw;p.write_text(new);assert p.read_text().endswith(old)
(backup/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='before',tofile='after')))
dump(O/'ENTRYPOINT_UPDATE.json',dict(status='PASS',path=str(p.relative_to(R)),before_sha256=hashlib.sha256(raw).hexdigest(),after_sha256=sha(p),history_preserved_exactly=True));print('ENTRYPOINT_UPDATED',sha(p))
