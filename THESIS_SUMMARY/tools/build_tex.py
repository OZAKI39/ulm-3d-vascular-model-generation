"""Single-source thesis export. Only converts manuscript markup; no science runs."""
from pathlib import Path
import re,json,shutil
P=Path(__file__).resolve().parents[1];source=P/'THESIS_DEVELOPMENT_SUMMARY.md'
md=source.read_text();lines=md.splitlines();emitted=[]
def esc(s):
 return ''.join({'&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}','\\':r'\textbackslash{}'}.get(c,c) for c in s)
def inline(s):
 pattern=r'(\$[^$\n]+\$|\[[^\]]+\]\([^\)]+\)|\*\*.*?\*\*)'
 out=[]
 for part in re.split(pattern,s):
  if part.startswith('$'):out.append(r'\('+part[1:-1]+r'\)')
  elif part.startswith('[') and '](' in part:
   label,url=re.fullmatch(r'\[([^\]]+)\]\(([^\)]+)\)',part).groups();out.append(r'\href{'+url.replace('%',r'\%')+'}{'+esc(label)+'}')
  elif part.startswith('**'):out.append(r'\textbf{'+esc(part[2:-2])+'}')
  else:out.append(esc(part))
 return ''.join(out)
header=r'''\documentclass[UTF8,a4paper,11pt,fontset=fandol]{ctexart}
\usepackage[margin=23mm]{geometry}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx,booktabs,longtable,array}
\usepackage[hidelinks]{hyperref}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\setlength{\parindent}{2em}
\setlength{\parskip}{0.35em}
\setlength{\emergencystretch}{3em}
\graphicspath{{figures/}}
\title{基于冻结三维血流场的微泡输运模型构建及 ULM 风格轨迹模拟}
\author{}
\date{}
\begin{document}
\maketitle
'''
out=[header];i=1;equations=[];figures=[];tables=[]
while i<len(lines):
 s=lines[i]
 if not s.strip():out.append('');i+=1;continue
 if s.startswith('## '):
  title=s[3:];m=re.match(r'\d+ (.+)',title)
  out.append(('\\section{'+inline(m.group(1))+'}') if m else ('\\section*{'+inline(title)+'}'));i+=1;continue
 if s=='$$':
  j=lines.index('$$',i+1);expr='\n'.join(lines[i+1:j]);tag=int(re.search(r'\\tag\{(\d+)\}',expr).group(1));plain=re.sub(r'\\tag\{\d+\}','',expr)
  out.append('\\begin{equation}\n'+plain+'\n\\label{eq:'+str(tag)+'}\n\\end{equation}');equations.append({'number':tag,'md_expression':plain});i=j+1;continue
 if s.startswith('!['):
  alt,path=re.fullmatch(r'!\[([^\]]+)\]\(([^\)]+)\)',s).groups();j=i+1
  while not lines[j].strip():j+=1
  cap=lines[j];n=int(re.match(r'图 (\d+) ',cap).group(1));caption=re.sub(r'^图 \d+ ','',cap)
  out.append('\\begin{figure}[htbp]\n\\centering\n\\includegraphics[width=0.98\\linewidth,height=0.69\\textheight,keepaspectratio]{'+path+'}\n\\caption{'+inline(caption)+'}\n\\label{fig:'+str(n)+'}\n\\end{figure}')
  figures.append({'number':n,'path':path,'caption':caption});i=j+1;continue
 if s.startswith('|'):
  rows=[]
  while i<len(lines) and lines[i].startswith('|'):
   line=lines[i]
   if not re.fullmatch(r'[|:\-\s]+',line):rows.append([x.strip() for x in line.strip('|').split('|')])
   i+=1
  n=len(rows[0]);weights={2:[.45,.45],3:[.25,.20,.45],4:[.32,.18,.18,.22]}.get(n,[.9/n]*n)
  specs=''.join('>{\\raggedright\\arraybackslash}p{'+str(w)+'\\linewidth}' for w in weights)
  out.append('{\\small\n\\setlength{\\tabcolsep}{4pt}\n\\begin{longtable}{'+specs+'}\n\\toprule')
  for ri,row in enumerate(rows):
   out.append(' & '.join(inline(x) for x in row)+r' \\')
   if ri==0:out.append(r'\midrule\endhead')
  out.append('\\bottomrule\n\\end{longtable}\n}');tables.append(rows);continue
 out.append(inline(s));i+=1
out.append(r'\end{document}')
tex='\n'.join(out)+'\n';(P/'THESIS_DEVELOPMENT_SUMMARY.tex').write_text(tex)
# Verify actual TeX equations, not merely the input data.
texeq=re.findall(r'\\begin\{equation\}\n(.*?)\n\\label\{eq:(\d+)\}\n\\end\{equation\}',tex,re.S)
assert [(int(n),expr.strip()) for expr,n in texeq]==[(x['number'],x['md_expression'].strip()) for x in equations]
assert [x['number'] for x in equations]==list(range(1,21))
assert [x['number'] for x in figures]==list(range(1,17))
assert all((P/x['path']).is_file() for x in figures)
# Recover the exact escaped paragraph stream from the TeX and compare every prose line.
prose=[s for s in lines if s.strip() and not s.startswith(('#','|','![','图 '))]
inmath=False;actualprose=[]
for s in lines:
 if s=='$$':inmath=not inmath;continue
 if not inmath and s.strip() and not s.startswith(('#','|','![','图 ')):actualprose.append(s)
assert all(inline(s) in tex for s in actualprose)
assert all(inline(x['caption']) in tex for x in figures)
assert all(inline(cell) in tex for rows in tables for row in rows for cell in row)
# User-requested technical-document wording scan; math command names are excluded.
body=md.split('## 参考文献')[0]
plain=re.sub(r'^\$\$\n.*?\n\$\$$','',body,flags=re.S|re.M);plain=re.sub(r'\$[^$\n]+\$','',plain)
pattern=r'\.py\b|\b(?:class|def|commit|branch|git|pytest|ssh|function|API|module)\b|json\s+field'
hits=[{'term':m.group(),'offset':m.start()} for m in re.finditer(pattern,plain,re.I)]
assert not hits,hits
stats=dict(numbered_sections=len(re.findall(r'^## \d+ ',md,re.M)),abstract_sections=1,bibliography_sections=1,
 chinese_characters_main_including_abstract=len(re.findall(r'[\u3400-\u9fff]',body)),
 nonwhitespace_characters_main_including_formulas=len(re.sub(r'\s','',body)),
 numbered_display_formulas=len(equations),main_figures=len(figures),tables=len(tables),
 technical_wording_check='PASS',technical_wording_hits=hits,
 md_tex_prose_parity='PASS',md_tex_formula_parity='PASS',md_tex_caption_table_parity='PASS',
 xelatex_available=bool(shutil.which('xelatex')),pdf_compiled=False,
 pdf_note='xelatex not installed; no LaTeX installation attempted. Static source checks do not establish compilation success.',
 count_definition='中文总字数按摘要及第1–18节的汉字字符计数，包含图注表格，排除参考文献和配套说明。')
(P/'audit/document_checks.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2)+'\n')
(P/'audit/formula_manifest.json').write_text(json.dumps(equations,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(stats,ensure_ascii=False))
