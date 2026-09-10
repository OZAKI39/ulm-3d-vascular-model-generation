"""Verify archived bytes and references without executing scientific code.

Run from any directory with Python 3. This reads data and parses Python syntax;
it does not import solver modules, compile native code or launch a browser.
"""
from pathlib import Path
import argparse, ast, collections, csv, hashlib, json, re, zipfile
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
        if p.suffix=='.py':ast.parse(b.decode(),filename=rel);checked['python_syntax']+=1
        if p.suffix=='.csv':
            with p.open(newline='') as f:
                rows=csv.reader(f,strict=True);header=next(rows,None)
                assert header and len(header)==len(set(header)),rel
                for line,row in enumerate(rows,start=2):
                    if not row:checked['blank_csv_rows']+=1;continue
                    assert len(row)==len(header),(rel,line)
                    csv_nonfinite[rel]+=sum(v.strip().lower() in ('nan','-nan','inf','-inf','infinity','-infinity') for v in row)
                    csv_blank[rel]+=sum(v=='' for v in row);checked['csv_data_rows']+=1
            checked['csv']+=1
        if p.suffix=='.npz':
            with zipfile.ZipFile(p) as z:assert z.testzip() is None and all(n.endswith('.npy') for n in z.namelist()),rel
            checked['npz_crc']+=1
        if b'\0' not in b:
            text=b.decode('utf-8',errors='replace')
            for kind,pattern in patterns.items():
                for match in pattern.finditer(text):
                    row=dict(path=rel,line=text.count('\n',0,match.start())+1,kind=kind)
                    expression=hashlib.sha256(match.group(1).encode()).hexdigest() if kind=='literal_secret_assignment' else None
                    allowance=review['reviewed_files'].get(rel)
                    if allowance and e['sha256']==allowance['file_sha256'] and expression in allowance['nonliteral_expression_sha256']:false_matches.append(row)
                    else:candidates.append(row)
            checked['sensitive_text_files']+=1
    if candidates:
        print(json.dumps(dict(status='SENSITIVE_REVIEW_REQUIRED',candidates=candidates),ensure_ascii=False,indent=2));raise SystemExit(2)
    def link_exists(page,link):
        url=urlsplit(link.strip('<>'))
        if url.scheme or link.startswith('#'):return
        p=(page.parent/unquote(url.path)).resolve()
        assert p.is_relative_to(ROOT) and p.exists(),(str(page),link)
        checked['relative_links']+=1
    for page in [ROOT/'BENCHMARK_REVIEW.md',REPORT/'README.md',REPORT/'NATIVE_SOURCE_ARCHIVE.md']:
        for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',page.read_text()):link_exists(page,link)
    data=ROOT/manifest['selected_result']['data_directory'];out=ROOT/manifest['selected_result']['output_directory']
    results=read(data/'comparison_results.json');receipt=read(data/'delivery_receipt.json')
    html_path=out/'comparison_review.html';html=html_path.read_text()
    embedded=html.split('<script id="repair-results" type="application/json">',1)[1].split('</script>',1)[0]
    assert json.loads(embedded)==results
    oldout=ROOT/'mirheo_starter/test_code/outputs/single_rbc_benchmark/rbc_shear_20260910'
    legacy=html.split('<script id="rbc-audit-data" type="application/json">',1)[1].split('</script>',1)[0]
    assert json.loads(legacy)==read(oldout/'results.json')
    assert read(data/'review_latest.json')['legacy_html_sha256']==digest(oldout/'single_rbc_review.html')
    assert 'id="repair-legacy"' in html and 'id="repair-current"' in html
    browser=read(out/'browser_cpu_delivery/browser_check.json')
    assert browser['html_sha256']==receipt['html_sha256']==digest(html_path)
    assert browser['status']=='PASS' and len(browser['checks'])==21 and all(browser['checks'].values())
    assert browser['checker_sha256']==digest(ROOT/'mirheo_starter/test_code/check_single_rbc_repair_browser.cjs')
    assets=Assets();assets.feed(html);assert not assets.assets,assets.assets
    for link in assets.links:link_exists(html_path,link)
    assert results['status']=='CPU_REPAIR_PREPARED_AWAITING_NEW_AUTHORIZATION'
    assert results['root_cause_status']=='SUPPORTED_NOT_CONFIRMED' and results['runtime_fix_status']=='NOT_TESTED'
    assert results['material_match']=='NOT_MATCHED' and results['qualified_speedup'] is None and results['human_review']=='PENDING'
    assert results['new_solver_s']==dict(gpu=0,cpu=0) and results['new_compile_s']==0 and not results['new_raw_solver_records']
    assert all(v is None for v in results['new_end_to_end_s'].values())
    assert read(data/'cpu_tests.json')['tests']==16 and read(data/'cpu_tests.json')['status']=='PASS'
    assert read(data/'cpu_tests.json')['log_sha256']==digest(data/'cpu_tests.log')
    request=read(data/'authorization_request.json');assert request['approved'] is False and not (data/'authorization.json').exists()
    plan=read(data/'frozen_benchmark_plan.json');assert request['frozen_plan_sha256']==digest(data/'frozen_benchmark_plan.json')
    mapping={e['source_path']:e for e in manifest['files'] if e.get('source_path')}
    excluded={e['source_path']:e for e in read(REPORT/'EXCLUDED_FILES.json')['files']}
    for name,value in plan['executable_artifact_sha256'].items():
        if name in mapping:assert digest(ROOT/mapping[name]['target_path'])==value,name
        else:assert name in excluded and excluded[name]['sha256']==value and '/.venv/' in name,name
        checked['frozen_artifacts']+=1
    for name,value in read(data/'root_cause_evidence.json')['source_sha256'].items():
        assert digest(ROOT/mapping[name]['target_path'])==value,name;checked['original_native_audit_sources']+=1
    for run in read(data/'failure_timeline.json')['runs']:
        for name,value in run['source_sha256'].items():
            source=str(Path(run['source_directory'])/name)
            assert digest(ROOT/mapping[source]['target_path'])==value,source;checked['timeline_hashed_inputs']+=1
        for sample in run['probe_summary']['per_frame']:
            assert sample['probe_file'] in mapping;checked['timeline_probe_references']+=1
        assert run['error']['exact_native_step'] is None and not run['common_target_reached']
    for e in read(REPORT/'RESULT_INDEX.json')['conclusion_evidence']:
        p=ROOT/e['json_file'];value=read(p)
        for segment in e['json_pointer'].lstrip('/').split('/'):
            segment=segment.replace('~1','/').replace('~0','~')
            value=value[int(segment)] if isinstance(value,list) else value[segment]
        assert value==e['value'],e;checked['conclusion_json_pointers']+=1
        for p in e.get('supporting_files',[]):assert (ROOT/p).is_file(),p
    allowed_empty={e['target_path'] for e in manifest['files'] if e.get('source_path') and (e['target_path'].endswith('.log') or e['target_path'].endswith('__init__.py'))}
    assert set(empty)<=allowed_empty,empty
    print(json.dumps(dict(status='PASS_ARCHIVE_INTEGRITY_ONLY',checks=dict(checked),empty_files_preserved=empty,nonfinite_csv_cells_preserved={k:v for k,v in csv_nonfinite.items() if v},blank_csv_cells_preserved={k:v for k,v in csv_blank.items() if v},reviewed_plotly_nonsecret_matches=false_matches,real_sensitive_candidates=[],original_html_and_browser_identity_match=True,embedded_repair_and_legacy_data_match=True,qualified_speedup=None,human_review='PENDING',scope='Only read/hash/syntax/reference checks; existing scientific and browser results are not rerun or promoted.'),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
