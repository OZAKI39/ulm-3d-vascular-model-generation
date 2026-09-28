"""Verify unchanged source bytes and internally consistent audit deliverables."""
from pathlib import Path
from html.parser import HTMLParser
import ast,csv,datetime,hashlib,json,platform
from PIL import Image
import numpy,scipy,pyvista,matplotlib

OUT=Path(__file__).resolve().parents[1];RUN=OUT.parent;PROJECT=RUN.parents[2]
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb') as fh:return hashlib.file_digest(fh,'sha256').hexdigest()
lock=read(OUT/'data/source_hashes.json')
for path,value in lock.items():assert sha(Path(path))==value,path
production=read(RUN/'data/final_summary.json')
for relative,value in production['identity']['protected_sources'].items():assert sha(PROJECT/relative)==value,relative
for relative,value in production['identity']['acceleration_sources'].items():assert sha(RUN/'scripts'/relative)==value,relative
summary=read(OUT/'data/analysis_summary.json')
assert summary['all_pass'] and summary['total']==1500
assert sum(summary['status_counts'].values())==1500
assert sum(summary['outlet_counts'].values())==1315
assert sum(r['count'] for r in summary['terminal_regions'])==185
assert read(OUT/'data/apparent_pauses_summary.json')['all_pass']
assert read(OUT/'data/three_forward_exceptions.json')['all_pass']
csv_counts={}
for name,n in [('all_tracks.csv',1500),('stopped_tracks.csv',185),('terminal_static_checks.csv',185),('apparent_pauses.csv',1500),('three_forward_exceptions.csv',3)]:
    count=sum(1 for _ in csv.DictReader((OUT/'data'/name).open()));assert count==n
    csv_counts[name]=count
figures={}
for name in ['01_stop_locations','02_stop_mechanism']:
    with Image.open(OUT/'figures'/f'{name}.png') as im:
        im.load();dpi=im.info.get('dpi');assert dpi and min(dpi)>=299
        assert im.width>=3000 and im.height>=2000
        figures[name]={'pixels':list(im.size),'dpi':list(dpi)}
    pdf=(OUT/'figures'/f'{name}.pdf').read_bytes();assert pdf.startswith(b'%PDF-') and b'%%EOF' in pdf[-32:]
for p in (OUT/'scripts').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
report=(OUT/'MICROBUBBLE_STOP_ANALYSIS_ZH.md').read_text();assert len(report)>5000
class LinkParser(HTMLParser):
    def __init__(self):super().__init__();self.targets=[]
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ['href','src']:self.targets.append(value)
parser=LinkParser();parser.feed((OUT/'OPEN_RESULTS.html').read_text())
for target in parser.targets:
    if target in ['data/final_verification.json','SHA256SUMS.txt']:continue
    assert (OUT/target).is_file(),target
result={'all_checks_pass':True,'scope':'Existing-data consistency and static model replay; not physiological trapping validation',
        'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'all_recorded_source_hashes_unchanged':True,'recorded_source_count':len(lock),
        'protected_source_count':len(production['identity']['protected_sources']),
        'track_files_verified_in_analysis':summary['track_files_hash_verified'],
        'csv_rows':csv_counts,'figure_files':figures,'html_links_checked':len(parser.targets),
        'figure_visual_review':'Reviewed both PNGs; moved figure 1 colorbar outside plot to prevent overlap',
        'new_CFD_runs':0,'new_trajectory_integrations':0,
        'terminal_static_solves':185,'supplementary_static_solves':3,
        'versions':{'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__,
                    'pyvista':pyvista.__version__,'matplotlib':matplotlib.__version__}}
(OUT/'data/final_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt' and 'mplconfig' not in p.parts)
(OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n' for p in files))
print(json.dumps(result,ensure_ascii=False,indent=2))
