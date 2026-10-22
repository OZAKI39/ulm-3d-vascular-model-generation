"""Check newly authored documentation and current figure links, not archived web pages."""
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import urlparse, unquote
from human_vessel_flow.paths import ROOT, write_json

class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):
        self.links.extend(v for k,v in attrs if k in {'src','href'} and v)

def main():
    files=[ROOT/'README.md',*(ROOT/'docs').glob('*.md'),ROOT/'diagnostics/README.md',ROOT/'scripts/README.md',
           ROOT/'results/reviews/migration_final/OPEN_RESULTS.html']
    rows=[]
    for path in files:
        text=path.read_text()
        if path.suffix=='.html':parser=Links();parser.feed(text);links=parser.links
        else:links=re.findall(r'\[[^\]]*\]\(([^)]+)\)',text)
        for value in links:
            link=urlparse(value.strip('<>'))
            if link.scheme or link.netloc or not link.path:continue
            target=(path.parent/unquote(link.path)).resolve()
            rows.append(dict(document=str(path.relative_to(ROOT)),link=value,exists=target.exists()))
    result=dict(status='PASS' if all(r['exists'] for r in rows) else 'FAIL',
                checked_documents=len(files),checked_links=len(rows),links=rows)
    write_json(ROOT/'validation/CURRENT_LINKS.json',result)
    print(result['status'],result['checked_documents'],result['checked_links'])
    for row in rows:
        if not row['exists']:print(row)
    if result['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
