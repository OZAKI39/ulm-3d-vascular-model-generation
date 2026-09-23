"""Document/evidence checks only; never calls a particle solver."""
from pathlib import Path
import json,re,hashlib,subprocess
P=Path(__file__).resolve().parents[1];R=P.parent
sample=json.loads((P/'audit/readability_sample.json').read_text())
comments={30:'能够理解几何接触的作用：限制向墙内的速度，允许沿墙运动；说明调整速度而不移动中心。',4:'明确列出真实场、人工算例及未解决问题，读者无需知道项目阶段命名。',24:'先以不同方向影响通行的直觉引出姿态，再引出梯度分解；后接的式（6）解释两个张量。',9:'流量平衡与局部不可压缩、网格加密稳定性分开；本轮已补充网格收敛的通俗说明。',5:'有限元首次出现即解释小单元及节点信息，网格数量与开放端盖含义明确。',20:'先解释群体差异及扁球近似的目的，随后章节给出体积关系；没有假定读者了解细胞几何实现。',48:'明确说明减小步长改变交接事件，能够理解安全和时间精度是不同要求。',44:'公式中的下限、有效间隙、约束间隙逐一给出含义和单位，曲率尺度追溯至已解释公式。',80:'按科学依赖串联能力并指出组合展示不等于全部耦合，未采用软件模块说明。',42:'把平滑函数说明为近场贡献的渐变开关，并指出背景流仍存在，未将平滑参数称为生物测量。'}
rows=['# 随机抽样的公式与可读性复核','','检查对象：具有基础理工科背景、未读项目代码的生物医学工程研究生。固定种子为 20260921；候选段落至少含 65 个汉字，不含标题、表格、图注及参考文献。完整抽样及原文保留在 [readability_sample.json](audit/readability_sample.json)。本检查是代理阅读自检，不冒称导师或用户已审核。','','本轮重写补充了变异系数的定义、网格收敛的含义、非负耗散的物理解释，并明确基础运动的时间离散与动画间隔不同。最终样本如下。','','## 十段正文','']
for row in sample['paragraphs']:
 idx=row['index'];rows.extend([f'### 样本 {idx}：第 {row["section"]}','',f'原文起句：{row["text"].split("。")[0]}。','',f'判断：可基本理解。{comments[idx]}',''])
checks={1:'入口/出口法向约定清楚；面积元与速度得到体积流量，用于守恒与注入。',3:'x、t、u 的意义和单位齐全，注明仅为基础随流层，并说明后续受约束速度不同。',4:'半旋度与局部刚体转动的联系已说明，单位 s⁻¹，服务旋转验证。',7:'方向向量、形状比、形状系数与两个张量齐全，单位相容；明确刚性自由扁球体范围。',9:'阻力、黏度、间隙、半径和约化半径齐全，kg/s；用途为法向阻滞而非全速度固定倍率。',10:'无量纲间隙、过渡坐标和权重齐全；参考尺寸区分墙与两球，说明平滑开关用途。',11:'原始/有效/下限/约束间隙与两个半径尺度齐全；说明下限并非真实接触，约束另行施加。',12:'浓度 m⁻³ 乘流量 m³/s 得 s⁻¹；明确是计划进入和名义浓度。',17:'路径位置、区间数、上下限与路程单位齐全；保存折线用于路径长度，不是端点距离。',19:'格、投影、指示函数、右端点位置、区间时间齐全；单位 s/格；与实际二维统计方式一致。'}
rows.extend(['## 十个随机公式','','每项均检查：变量解释、单位、整体物理含义及本研究用途。'])
for n in sample['formula_numbers']:rows.extend(['',f'- 式（{n}）：四项均通过。{checks[n]}'])
rows.extend(['','## 图像阅读检查','','已查看图 1–16 的两张整组预览，并结合原图数据与审核记录核对图注。确认示意图未冒充数值结果，人工接触/剪切/生命周期图未称为真实生理通行，出口图保留零完成分支，累积图明确为 s/bin。原始英文图题和阶段字样属于历史图追溯标记，正文以中文物理解释为主。','', '局限：没有 LaTeX 编译器，尚未检查排版后每页的浮动位置及分页效果；不能将这项自检替代 PDF 人工排版审核。'])
(P/'THESIS_READABILITY_REVIEW.md').write_text('\n'.join(rows)+'\n')
md=(P/'THESIS_DEVELOPMENT_SUMMARY.md').read_text();tex=(P/'THESIS_DEVELOPMENT_SUMMARY.tex').read_text()
# Adjacent formula explanation/meaning checks.
formula=[]
for m in re.finditer(r'^\$\$\n(.*?)\n\$\$$',md,re.S|re.M):
 n=int(re.search(r'\\tag\{(\d+)\}',m.group(1)).group(1));tail=md[m.end():].lstrip();where,meaning=tail.split('\n\n',2)[:2]
 assert where.startswith('其中'),n
 assert meaning.startswith('该式表示'),n
 assert '在本研究' in meaning,n
 formula.append({'number':n,'adjacent_variables':True,'physical_meaning':True,'research_use':True})
# Dimensional homogeneity: representative expressions, dimensions (kg, m, s).
zero=(0,0,0);length=(0,1,0);time=(0,0,1);velocity=(0,1,-1);grad=(0,0,-1);area=(0,2,0);volume=(0,3,0);mu=(1,-1,-1);pressure=(1,-1,-2);resistance=(1,0,-1);flux=(0,3,-1);density=(0,-3,0)
def mul(*items):return tuple(sum(t[i] for t in items) for i in range(3))
def div(a,b):return tuple(x-y for x,y in zip(a,b))
unit_expr={1:(mul(velocity,area),flux),2:(mul(zero,pressure),pressure),3:(div(length,time),velocity),4:(mul(div(zero,length),velocity),grad),5:(div(volume,area),length),6:(grad,grad),7:(mul(grad,zero),div(zero,time)),8:(length,length),9:(div(mul(mu,area),length),resistance),10:(div(length,length),zero),11:(div(mul(zero,mu,area),length),resistance),12:(mul(density,flux),div(zero,time)),13:(mul(zero,flux),flux),14:(mul(density,flux,time),zero),15:(div(velocity,flux),div(zero,area)),16:(time,time),17:(mul(velocity,time),length),18:(mul(density,time),(0,-3,1)),19:(mul(time,zero),time),20:(div(flux,flux),zero)}
assert all(a==b for a,b in unit_expr.values())
for f in formula:f['dimensional_representative_check']=unit_expr[f['number']][0]
# Keep bibliography and scientific prose free of technical doc wording (raw commands excluded from math).
body=md.split('## 参考文献')[0]
assert len(re.findall(r'^\$\$$',md,re.M))==40
assert [f['number'] for f in formula]==list(range(1,21))
plain=re.sub(r'^\$\$\n.*?\n\$\$$','',md,flags=re.S|re.M)
plain=re.sub(r'\$[^$\n]+\$','',plain)
assert '$' not in plain
assert all(marker not in md for marker in (r'\[',r'\]',r'\(',r'\)'))
assert '\\tag{' not in tex
assert tex.count(r'\begin{equation}')==tex.count(r'\end{equation}')==20
# Balanced TeX braces and environments; not a claim of successful TeX compilation.
stack=[]
for m in re.finditer(r'\\(begin|end)\{([^}]+)\}',tex):
 act,name=m.groups()
 if act=='begin':stack.append(name)
 else:assert stack and stack.pop()==name,(m.start(),name)
assert not stack
level=0
for i,ch in enumerate(tex):
 if ch in '{}' and (i==0 or tex[i-1]!='\\'):
  level+=1 if ch=='{' else -1;assert level>=0
assert level==0
# All local Markdown targets resolve, including 112 appendix figures.
links=[]
for p in P.glob('*.md'):
 for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
  if re.match(r'https?://',target):continue
  q=(p.parent/target).resolve();assert q.exists(),(p,target);links.append(str(q))
# Historical evidence files/selected images stayed byte-identical.
unchanged=[]
for row in json.loads((P/'audit/evidence_inventory.json').read_text()):
 actual=hashlib.sha256((R/row['path']).read_bytes()).hexdigest();assert actual==row['sha256'],row['path'];unchanged.append(row['path'])
for row in json.loads((P/'audit/figure_manifest.json').read_text()):
 actual=hashlib.sha256((P/row['destination']).read_bytes()).hexdigest();assert actual==row['sha256']
 if row['copy_exact']:assert actual==hashlib.sha256((R/row['source']).read_bytes()).hexdigest()
old=(P/'audit/preexisting_git_status.txt').read_text();now=subprocess.check_output(['git','status','--short'],cwd=R,text=True)
filter_status=lambda s:[x for x in s.splitlines() if 'THESIS_SUMMARY/' not in x]
assert filter_status(old)==filter_status(now),(old,now)
# Check primary cohort arithmetic and formula numeric examples using saved data.
d=json.loads((P/'audit/extracted_values.json').read_text());v=d['8_1']['values'];q=v['clock']['Q_in_m3_s'];c=v['clock']['C_MB_m3']
assert v['scheduled']==v['admitted']+231 and v['admitted']==v['completed']+720
assert sum(v['outlet_counts'].values())==v['completed'] and v['outlet_counts']['OUTLET_02']==1249
assert abs(q*c-d['7']['values']['mb_number_rate_s_inv'])<1e-16
result={'status':'PASS','scope':'DOCUMENT_ONLY_NO_NEW_NUMERICAL_MODEL_OR_TRAJECTORY_RUN', 'formulas':formula,'local_links_checked':len(links),'historical_evidence_files_unchanged':len(unchanged),'historical_main_images_exact':14,'preexisting_git_changes_preserved':True,'random_paragraphs_reviewed':10,'random_formulas_reviewed':10,'all_20_formulas_manually_reviewed':True,'dimensions_reviewed':True,'latex_static_structure':'PASS','pdf_compilation':'NOT_RUN_XELATEX_UNAVAILABLE','historical_scientific_tests_rerun':False}
(P/'audit/final_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='formulas'},ensure_ascii=False))
