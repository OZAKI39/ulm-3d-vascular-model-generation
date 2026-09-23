from pathlib import Path
import json,re,random,hashlib
P=Path(__file__).resolve().parents[1];R=P.parent
manifest=json.loads((P/'audit/figure_manifest.json').read_text());md=(P/'THESIS_DEVELOPMENT_SUMMARY.md').read_text()
captions={int(n):s for n,s in re.findall(r'^图 (\d+) (.+)$',md,re.M)}
sections={1:1,2:2,3:3,4:4,5:5,6:6,7:7,8:8,9:9,10:10,11:10,12:12,13:12,14:14,15:15,16:16}
rows=['# 正文主图与建议附录图索引','','正文共 16 张主图：14 张历史原图逐字节复制，2 张本轮整理图。原图数据不变；本轮没有启动新的轨迹计算。英文历史图题和阶段标识保留，正文使用中文学位论文式图注。原图记号与论文换记关系见符号表。','','## 正文主图','', '| 图 / 节 | 图注标题与阅读重点 | 证据类别 | 历史来源 / 整理方式 |','|---|---|---|---|']
selected=set()
for f in manifest:
 n=f['number'];source=f['source'];selected.add(source)
 link=f'[{Path(source).name}](../{source})' if (R/source).is_file() else source
 rows.append(f'| [图 {n}](figures/fig_{n:02d}.png) / 第 {sections[n]} 节 | {captions[n]} | {f["evidence_type"]} | {link}；'+('原样复制' if f['copy_exact'] else '本轮绘制，非新数值实验')+' |')
rows+=['','图 1 另有 [矢量 PDF](figures/fig_01.pdf)，图 16 另有 [矢量 PDF](figures/fig_16.pdf)。这两张单图 PDF 由绘图工具导出，不是全文 LaTeX 编译预览。','', '## 建议附录图','', '以下列出未放入正文的历史静态图，按原阶段目录便于追溯。建议围绕具体审核问题择用，而非全部插入论文；含 scope、environment、audit、limitations 的图更适合技术补充材料。真实 FEM 与 synthetic 的判别须结合对应原 REVIEW，不凭文件名推断生理证据。']
allfig=[]
for root in [R/'particle_3d/reports',R/'particle_3d/outputs']:
 for f in root.glob('*/figures/*.png'):
  rel=str(f.relative_to(R))
  if 'particle8_2' in rel or rel in selected:continue
  allfig.append(rel)
last=''
for rel in sorted(allfig):
 stage=Path(rel).parts[2]
 if stage!=last:rows+=['',f'### {stage}',''];last=stage
 rows.append(f'- [{Path(rel).stem}](../{rel})')
rows+=['','## 动画补充材料','','动画不计入正文 16 张静态图，也不因显示帧数增加轨迹样本数。真实场景、人工控制及真实未通行的混合接纳场景应保持原分类。']
for rel,title in [('particle_3d/outputs/particle8_1/index.html','完整轨迹、累积和出口着色动画总览'),('particle_3d/outputs/particle8_full3d/index.html','此前完整几何与三维旋转回放总览')]:
 rows.append(f'- [{title}](../{rel})')
(P/'THESIS_FIGURE_INDEX.md').write_text('\n'.join(rows)+'\n')
# Evidence status, including historical test records without claiming they ran this round.
stage=[]
for p in sorted(list((R/'particle_3d/reports').glob('particle*/*VALIDATION.json'))+list((R/'particle_3d/outputs').glob('particle*/*VALIDATION.json'))):
 if 'particle8_2' in str(p):continue
 d=json.loads(p.read_text());entry={'source':str(p.relative_to(R)),'status':{k:d[k] for k in ['automated_checks','automated_checks_pass','stage_result','statuses','manual_visual_review','user_acceptance'] if k in d},'historical_tests':{k:d[k] for k in ['tests','particle1_tests','test_results','regressions','automatic_tests','total_tests','total_tests_run_this_stage'] if k in d},'rerun_this_task':False};stage.append(entry)
(P/'audit/stage_validation_index.json').write_text(json.dumps(stage,ensure_ascii=False,indent=2)+'\n')
# Fixed-seed random sampling with saved excerpts; subsequent reading is recorded separately.
blocks=[];section='摘要'
for block in md.split('\n\n'):
 if block.startswith('## '):section=block[3:];continue
 if section.startswith('参考文献'):break
 if not block.startswith(('#','|','![','图 ','$$')) and len(re.findall('[\u3400-\u9fff]',block))>=65:
  blocks.append({'section':section,'text':block})
rng=random.Random(20260921);picks=rng.sample(range(len(blocks)),10);eqpicks=sorted(rng.sample(range(1,21),10))
sample={'seed':20260921,'paragraph_population':len(blocks),'paragraphs':[dict(index=i,**blocks[i]) for i in picks],'formula_numbers':eqpicks}
(P/'audit/readability_sample.json').write_text(json.dumps(sample,ensure_ascii=False,indent=2)+'\n')
print('Appendix figures',len(allfig),'Stage records',len(stage),'Paragraph sample',picks,'Formula sample',eqpicks)
for i in picks:print(i,blocks[i]['section'],blocks[i]['text'])
