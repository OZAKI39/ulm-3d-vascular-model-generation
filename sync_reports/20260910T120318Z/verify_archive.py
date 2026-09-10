"""Read archived bytes, tables, run seals and links without running any solver.

This checks archive integrity only. No solver imports, scientific tests,
compilation, recalibration or browser execution take place.
"""
from pathlib import Path
import argparse, ast, collections, csv, hashlib, json, math, re, zipfile
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit

ROOT=Path(__file__).resolve().parents[2]
REPORT=Path(__file__).resolve().parent
def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
class Assets(HTMLParser):
    def __init__(self):super().__init__();self.assets=[];self.links=[]
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if tag in ['script','img','iframe','video','audio','source'] and a.get('src'):self.assets.append(a['src'])
        if tag=='link' and a.get('href'):self.assets.append(a['href'])
        if tag=='a' and a.get('href'):self.links.append(a['href'])
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check-sources',action='store_true');args=parser.parse_args()
    manifest=read(REPORT/'SYNC_MANIFEST.json');review=read(REPORT/'SENSITIVE_SCAN_REVIEW.json')
    checked=collections.Counter();empty=[];csv_nonfinite=collections.Counter();csv_blank=collections.Counter();candidates=[];false_matches=[]
    patterns={
        'github_token':re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,255}|github_pat_[A-Za-z0-9_]{60,255})\b'),
        'private_key':re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----'),
        'aws_access_key':re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
        'openai_key':re.compile(r'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{30,}\b'),
        'credential_url':re.compile(r'https?://[^\s/<>:"\']+:[^\s/<>@"\']+@'),
        'bearer_credential':re.compile(r'(?i)authorization["\']?\s*[:=]\s*["\']?bearer\s+[A-Za-z0-9._~-]{20,}'),
        'literal_secret_assignment':re.compile(r'''(?ix)\b(?:github_token|access_token|api_key|secret_key|password)\b ["']?\s*[:=]\s*["']([^"'\s]{12,})["']''')}
    for e in manifest['files']:
        rel=e['target_path'];p=ROOT/rel;assert p.is_file() and not p.is_symlink(),rel
        assert p.stat().st_size==e['size_bytes'] and digest(p)==e['sha256'],rel
        assert bool(p.stat().st_mode&0o111)==e['executable'],rel;checked['file_hashes_and_modes']+=1
        if args.check_sources and e.get('source_path'):
            src=Path(e['source_path']);assert digest(src)==e['sha256'] and bool(src.stat().st_mode&0o111)==e['executable'],rel;checked['source_hashes_and_modes']+=1
        assert not set(Path(rel).parts)&{'.git','.venv','__pycache__','node_modules','.ssh'},rel
        assert p.suffix not in ('.so','.a','.o','.exe','.lock','.pyc') and not p.name.startswith('.env'),rel
        b=p.read_bytes();assert not b.startswith((b'version https://git-lfs.github.com/spec/v1',b'\x7fELF')),rel
        if not b:empty.append(rel)
        if p.suffix=='.json':json.loads(b);checked['json']+=1
        if p.suffix=='.jsonl':
            for line in b.splitlines():
                if line.strip():json.loads(line);checked['jsonl_records']+=1
        if p.suffix=='.py':ast.parse(b.decode(),filename=rel);checked['python_syntax']+=1
        if p.suffix=='.csv':
            with p.open(newline='') as f:
                rows=csv.reader(f,strict=True);header=next(rows,None)
                if header is not None:
                    assert header and len(header)==len(set(header)),rel
                    for line,row in enumerate(rows,start=2):
                        if not row:checked['blank_csv_rows']+=1;continue
                        assert len(row)==len(header),(rel,line)
                        csv_nonfinite[rel]+=sum(v.strip().lower() in ('nan','-nan','inf','-inf','infinity','-infinity') for v in row)
                        csv_blank[rel]+=sum(v=='' for v in row);checked['csv_data_rows']+=1
            checked['csv']+=1
        if p.suffix=='.npz':
            with zipfile.ZipFile(p) as z:
                assert z.testzip() is None,rel;assert all(n.endswith('.npy') for n in z.namelist()),rel
            checked['npz_crc']+=1
        if b'\0' not in b:
            text=b.decode('utf-8',errors='replace')
            for kind,pattern in patterns.items():
                for match in pattern.finditer(text):
                    row=dict(path=rel,line=text.count('\n',0,match.start())+1,kind=kind)
                    expression=hashlib.sha256(match.group(1).encode()).hexdigest() if kind=='literal_secret_assignment' else None
                    if (rel==review['path'] and e['sha256']==review['file_sha256'] and expression in review['nonliteral_expression_sha256']):false_matches.append(row)
                    else:candidates.append(row)
            checked['sensitive_text_files']+=1
    if candidates:
        print(json.dumps(dict(status='SENSITIVE_REVIEW_REQUIRED',candidates=candidates),ensure_ascii=False,indent=2));raise SystemExit(2)
    def link_exists(page,link):
        url=urlsplit(link.strip('<>'))
        if url.scheme or link.startswith('#'):return
        assert (page.parent/unquote(url.path)).exists(),(str(page),link)
        checked['relative_links']+=1
    for page in [ROOT/'BENCHMARK_REVIEW.md',REPORT/'README.md']:
        for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',page.read_text()):link_exists(page,link)
    final=ROOT/manifest['selected_result']['repository_directory'];page=final/'single_rbc_review.html';html=page.read_text();results=read(final/'results.json')
    embedded=html.split('<script id="rbc-audit-data" type="application/json">',1)[1].split('</script>',1)[0]
    assert json.loads(embedded)==results
    browser=read(final/'browser_delivery_verified/browser_check.json')
    assert browser['html_sha256']==digest(page) and browser['status']=='PASS' and len(browser['checks'])==19 and all(browser['checks'].values())
    assets=Assets();assets.feed(html);assert not assets.assets,assets.assets
    for link in assets.links:link_exists(page,link)
    assert results['workflow']=='PARTIAL' and results['model_comparability']=='PARTIAL' and results['benchmark_screen']=='FAILED'
    assert results['qualified_speedup'] is None and results['human_review']=='PENDING' and results['formal_cold_completed']==2
    data=ROOT/'mirheo_starter/data/single_rbc_benchmark/rbc_shear_20260910';runs=ROOT/'mirheo_starter/runs/single_rbc_benchmark/rbc_shear_20260910'
    delivery=read(data/'final_delivery.json');assert delivery['html_sha256']==digest(page) and delivery['results_sha256']==digest(final/'results.json')
    mapping={e['source_path']:e['target_path'] for e in manifest['files'] if e.get('source_path')}
    for src,sha in delivery['source_sha256'].items():assert digest(ROOT/mapping[src])==sha,src
    records=[]
    for seal in runs.rglob('output_sha256.json'):
        for name,value in read(seal).items():assert digest(seal.parent/name)==value,(str(seal),name)
        checked['sealed_runs']+=1
    for p in runs.rglob('execution.json'):records.append((p,read(p)))
    assert checked['sealed_runs']==34 and len(records)==34
    ledger=read(runs/'solver/budget_ledger.json');assert len(ledger['attempts'])==26
    charged=sum(a['charged_s'] for a in ledger['attempts']);execution_sum=sum(r['elapsed_monotonic_s'] for p,r in records if '/solver/' in str(p))
    assert math.isclose(charged,execution_sum,abs_tol=1e-9) and math.isclose(charged,results['solver_charged_s'],abs_tol=1e-9)
    for r in results['runs']:
        p=ROOT/r['directory'].removeprefix('/home/lzy/projects/')
        assert read(p/'execution.json')==r['execution'],r['task']
    assert read(data/'cpu_tests_final.json')['passed']==26
    print(json.dumps(dict(status='PASS_ARCHIVE_INTEGRITY_ONLY',checks=dict(checked),empty_files_preserved=empty,nonfinite_csv_cells_preserved={k:v for k,v in csv_nonfinite.items() if v},blank_csv_cells_preserved={k:v for k,v in csv_blank.items() if v},reviewed_plotly_nonsecret_matches=false_matches,real_sensitive_candidates=[],original_html_and_browser_identity_match=True,embedded_html_equals_results=True,ledger_equals_execution_and_results=True,scientific_workflow_preserved='PARTIAL',benchmark_screen_preserved='FAILED',qualified_speedup=None,human_review='PENDING',scope='Read/hash/syntax/reference verification only. No new simulation, scientific tests, compilation, recalibration or browser run.'),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
