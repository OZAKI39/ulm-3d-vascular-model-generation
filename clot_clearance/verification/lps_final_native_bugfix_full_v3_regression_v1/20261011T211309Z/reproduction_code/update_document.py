from pathlib import Path
import json,hashlib,difflib,re
C=json.loads(Path(__file__).with_name('config.json').read_text());L=Path(C['local']);O=L/'verification'/C['task']/C['run'];F=L/'results'/C['task']/C['run'];doc=Path('/home/lzy/projects/clot_rupture.md');summary=json.loads((O/'SCIENTIFIC_SUMMARY.json').read_text());status=json.loads((O/'FINAL_STATUS.json').read_text())
assert summary['complete'] and status['FULL_V3_REGRESSION']=='FULL_V3_NATIVE_REGRESSION_PASS'
sha=lambda b:hashlib.sha256(b).hexdigest();before=doc.read_bytes();text=before.decode();before_hash=sha(before);expected='5e1db85d6a99285c92565e5a62795d5b22e130ed0e0706810401c84afb9d5eca';assert before_hash==expected,'Concurrent document update: reread before adapting edits'
prov=O/'document_provenance';prov.mkdir(exist_ok=False);(prov/'clot_rupture.before.md').write_bytes(before)
needle='**最终版本未再跑V3，完整动力学认证尚未完成。**';assert text.count(needle)==1;text=text.replace(needle,'**截至该轮结束，最终版本未再跑V3，完整动力学认证尚未完成；后续完整回归见下段。**')
insert='''**本轮最终双修复完整V3回归（20261011T211309Z）**：直接复制原C1接受输入，G、K、s00、200 Pa脉冲、0.5 µs步长、1.0Δx接触及原背景阻力均不变。实际进程映射确认加载上轮最终双修复诊断库（SHA `00b246ebe422…22ae9`），唯一正式运行完成4660步、2.33 ms，两MPI均正常退出。全程4661个状态逐ID检查owned和ghost的固定参考m，最大相对误差1.12×10⁻¹⁵；实际非恒等排序发生在1018、2004、3200、4023步，共34次原子ownership转移。1039条跨rank断键均核对两端owned记录、原生事件及最终存活键，没有漏端点或重复记录。生产库、env.sh和2703个历史/生产保护文件SHA未变，**此结论仅认证本次诊断回归，不是自动部署**。

首个可对齐的力学差异在455步（0.2275 ms）：坐标、半步速度和拓扑仍逐位相同，但最大力差11.650 µN。该步发生真实迁移、没有排序。按旧源码不随ID搬运m的行为，从实际本地顺序重建出386和3480号点的m错配，重建旧力与历史保存力的相对差仅4.90×10⁻¹⁶；这是源码重建支持的归因，旧运行没有直接保存m。455步末半步积分已使接受速度分叉；456步出现位置及半步速度差，462步首次断键事件分叉。旧记录没有覆盖的步仍不补造因果证据。

本体及界面首次断键仍为241/621步；首次真实脱离从旧593步变为611步（0.3055 ms），当步启用背景阻力。关键完整精度状态的原生力与“断前theta、断后存活键＋原生接触”预测最大相对差2.65×10⁻¹⁴。**显式theta时间层限制仍在**，原生PE计量限制也未消除；材料Gc没有标定。粒子接触和早期阻力得到覆盖，管壁罚接触、V3 restart读回、原子数组增长及后续5秒输运本轮未触发或未执行，不记为全覆盖通过。

'''
needle='# ⑥ 外加剪切：模拟微泡对血块表面的“刮”';assert text.count(needle)==1;text=text.replace(needle,insert+needle)
needle='最早的一次脱落发生在 0\\.2965 ms（节点 267）。';assert text.count(needle)==1;text=text.replace(needle,'历史接受C1最早的一次脱落发生在 0\\.2965 ms（节点 267）；本轮最终双修复回归为0.3055 ms，两版各自按真实拓扑事件启用阻力。')
needle='新版末态碎块速度接近血流';assert text.count(needle)==1;text=text.replace(needle,'历史early-flow版末态碎块速度接近血流')
needle='计算时间变成约 3 倍的历史记录保留，原因是每一步都要调用一次 Python。';assert text.count(needle)==1;text=text.replace(needle,needle+'\n\n**本轮逐步实测**：只读观察器在原生final_integrate之后记录每步接受速度，最大3.913394 m/s（759步、0.3795 ms、364号点）；最大步间位移1.95708 µm。全体材料点动能最大4.37535×10⁻⁸ J，游离材料点动能最大2.06164×10⁻⁸ J。原生求解循环187.013 s，含启动189.286 s；两rank进程树采样峰值约1.913 GB。新值是每步求解器记录，旧3.91997 m/s仍只是历史已保存时刻的最大值，两者不混用。')
needle='碎块大小用等效球直径';assert text.count(needle)==1
newstats='''**本轮修复后结果与历史并列**（均为2.33 ms，ATTACHMENT_AWARE_V2）：

| 指标 | 历史接受C1 | 最终双修复诊断回归 |
|---|---:|---:|
| 游离材料点 / 游离组件 | 217 / 208 | 199 / 190 |
| 仍附着材料点 / 质量比例 | 1710 / 88.739% | 1728 / 89.673% |
| 单点 / 双点 / 三点 / 四点游离组件 | 201 / 6 / 0 / 1 | 184 / 4 / 1 / 1 |
| 1–2点碎屑占游离质量 | 98.157% | 96.482% |
| 末态最大真实游离块 | 4点，98.47 µm | 4点，98.47 µm |
| 全过程最大真实游离块 | 13点 | 20点，168.39 µm |
| ≥200 / ≥300 µm真实游离组件 | 0 / 0 | 0 / 0 |
| 累计bulk / interface断键 | 11771 / 718 | 10808 / 643 |

两版均守恒1927个材料点。修复确实改变轨迹和断键，但修复后仍以未解析细碎屑为主，**不足以认定旧细碎屑主要由该bug造成**；20点瞬时组件增大也不是实验认证。下面217点、208组件的叙述继续保留为历史结果。

'''
text=text.replace(needle,newstats+needle)
oldlast='此前静态加载接口建议暂缓；**下一轮唯一优先项**是使用最终“参考m身份＋owned端点检查”诊断库，对原C1工况做一次带逐步ID/迁移/力日志的完整V3验证。本轮不再追加运行，不自动调材料或启动CFD/微泡/AI任务。';assert text.count(oldlast)==1
review=O/'REVIEW_LPS_FINAL_NATIVE_BUGFIX_FULL_V3_REGRESSION_V1.md';html=F/'OPEN_RESULTS.html'
text=text.replace(oldlast,f'该轮提出的完整V3验证已由本轮完成，见[完整回归审查与37项回答]({review})和[四组真实图与数据]({html})。**下一轮唯一优先项**：冻结最终双修复库，做受控小试样的断键时间层与时间步敏感性审查，区分正常拓扑能量释放与显式theta滞后的数值影响。本轮到此停止，不自动部署、调材料或启动新计算。')
after=text.encode();assert re.findall(r'^# [①-⑨].*$',before.decode(),re.M)==re.findall(r'^# [①-⑨].*$',text,re.M);assert '217 个点脱落' in text and '208 个碎块里' in text and 'G_{\\mathrm{drive}}' in text
assert sha(doc.read_bytes())==before_hash,'Concurrent edit while preparing document';doc.write_bytes(after);assert doc.read_bytes()==after
(prov/'clot_rupture.after.md').write_bytes(after);(prov/'clot_rupture.diff').write_text(''.join(difflib.unified_diff(before.decode().splitlines(True),text.splitlines(True),str(doc)+' before',str(doc)+' after')))
audit={'status':'PASS_ACTUAL_WSL_UPDATED','actual_path':str(doc),'before_sha256':before_hash,'after_sha256':sha(after),'backup':str(prov/'clot_rupture.before.md'),'diff':str(prov/'clot_rupture.diff'),'readback_equal':True,'concurrent_SHA_check':'Matched immediately before write','numbered_titles_preserved':True,'historical217_208_preserved':True,'material_Gc_uncalibrated_preserved':True}
(O/'DOCUMENT_UPDATE_AUDIT.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False));status['WSL_DOCUMENT_UPDATE']='PASS_ACTUAL_WSL_UPDATED';(O/'FINAL_STATUS.json').write_text(json.dumps(status,indent=2,ensure_ascii=False));print(json.dumps(audit,ensure_ascii=False))
