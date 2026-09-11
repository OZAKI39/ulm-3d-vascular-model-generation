"""Read and hash archived files without importing a solver or running experiments."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import argparse, ast, collections, csv, hashlib, json, re, zipfile

ROOT=Path(__file__).resolve().parents[2]
REPORT=Path(__file__).resolve().parent

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def read(path):return json.loads(Path(path).read_text())

PATTERNS={
    'github_token':rb'\b(?:gh[pousr]_[A-Za-z0-9]{30,255}|github_pat_[A-Za-z0-9_]{60,255})\b',
    'private_key':rb'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----',
    'aws_access_key':rb'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b',
    'service_key':rb'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{30,}\b',
    'credential_url':rb'https?://[^\s/<>:\x22\x27]+:[^\s/<>@\x22\x27]+@',
    'bearer_credential':rb'(?i)authorization[\x22\x27]?\s*[:=]\s*[\x22\x27]?bearer\s+[A-Za-z0-9._~-]{20,}',
    'literal_secret_assignment':rb'''(?ix)\b(?:github_token|access_token|api_key|secret_key|password|authorization)\b ["']?\s*[:=]\s*["']([^"'\s]{12,})["']''',
}

class Assets(HTMLParser):
    def __init__(self):super().__init__();self.assets=[]
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if tag in ['script','img','iframe','video','audio','source'] and a.get('src'):self.assets.append(a['src'])
        if tag=='link' and a.get('href'):self.assets.append(a['href'])

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check-sources',action='store_true');parser.add_argument('--output');args=parser.parse_args()
    manifest=read(REPORT/'SYNC_MANIFEST.json');count=collections.Counter();empty=[];csv_empty=[];issues=[];sensitive=[];nonfinite=[]
    allowances=read(REPORT/'SENSITIVE_SCAN_REVIEW.json').get('reviewed_matches',[]) if (REPORT/'SENSITIVE_SCAN_REVIEW.json').exists() else []
    for e in manifest['files']:
        relative=e['target_path'];p=ROOT/relative
        assert p.is_file() and not p.is_symlink(),relative
        assert p.stat().st_size==e['size_bytes'] and digest(p)==e['sha256'],relative
        assert bool(p.stat().st_mode&0o111)==e['executable'],relative
        assert not set(Path(relative).parts)&{'.git','.venv','__pycache__','node_modules','.ssh'},relative
        count['file_hashes_and_modes']+=1
        if args.check_sources and e.get('source_path'):
            s=Path(e['source_path']);assert digest(s)==e['sha256'] and s.stat().st_mtime_ns==e['source_mtime_ns'],relative
            assert bool(s.stat().st_mode&0o111)==e['executable'],relative;count['unchanged_source_files']+=1
        b=p.read_bytes();assert not b.startswith((b'version https://git-lfs.github.com/spec/v1',b'\x7fELF')),relative
        if not b:empty.append(relative)
        try:
            if p.suffix=='.json':json.loads(b);count['json']+=1
            if p.suffix=='.jsonl':
                for line in b.splitlines():
                    if line.strip():json.loads(line);count['jsonl_records']+=1
            if p.suffix=='.py':ast.parse(b.decode('utf-8'),filename=relative);count['python_ast']+=1
            if p.suffix=='.csv':
                with p.open(newline='') as f:
                    rows=csv.reader(f,strict=True);header=next(rows,None);n=0;nf=0
                    if header is None:csv_empty.append(relative)
                    else:
                        for line,row in enumerate(rows,start=2):
                            if not row:continue
                            assert len(row)==len(header),(relative,line,'column count')
                            nf+=sum(v.strip().lower() in ['nan','inf','-inf','-nan'] for v in row);n+=1
                        if not n:csv_empty.append(relative)
                    if nf:nonfinite.append(dict(path=relative,nonfinite_values=nf))
                    count['csv_data_rows']+=n;count['csv']+=1
            if p.suffix=='.npz':
                with zipfile.ZipFile(p) as archive:assert archive.testzip() is None
                count['npz_crc']+=1
            if p.suffix in ['.h5','.hdf5']:
                assert b.startswith(b'\x89HDF\r\n\x1a\n'),relative;count['hdf5_signatures']+=1
            if p.suffix=='.xmf':
                import xml.etree.ElementTree as ET
                for item in ET.fromstring(b).iter('DataItem'):
                    if item.attrib.get('Format')=='HDF':
                        name=(item.text or '').strip().split(':',1)[0];q=(p.parent/name).resolve()
                        assert q.is_relative_to(ROOT) and q.is_file(),(relative,name);count['xmf_hdf5_links']+=1
        except Exception as error:issues.append(dict(path=relative,error=type(error).__name__+': '+str(error)))
        if b'\0' not in b:
            for kind,expression in PATTERNS.items():
                for m in re.finditer(expression,b):
                    hit=dict(path=relative,line=b[:m.start()].count(b'\n')+1,kind=kind,match_sha256=hashlib.sha256(m.group(0)).hexdigest(),file_sha256=e['sha256'])
                    if not any(all(a.get(k)==hit[k] for k in hit) for a in allowances):sensitive.append(hit)
            count['sensitive_text_files']+=1
    for p in [ROOT/'BENCHMARK_REVIEW.md',REPORT/'README.md',REPORT/'NATIVE_SOURCE_ARCHIVE.md']:
        for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            url=urlsplit(link.strip('<>'))
            if url.scheme or link.startswith('#'):continue
            q=(p.parent/unquote(url.path)).resolve();assert q.is_relative_to(ROOT) and q.exists(),(str(p),link);count['relative_links']+=1
    raw=ROOT/'cloud_results/cloud-rbc-full-20260911T000225Z';review=raw.with_name(raw.name+'-review')
    source_manifest=read(raw/'results_manifest.json')
    for row in source_manifest['files']:
        p=raw/row['path'];assert digest(p)==row['sha256'] and p.stat().st_size==row['size'],row['path'];count['cloud_raw_manifest_files']+=1
    assert read(raw/'RESULTS_READY')['manifest_sha256']==digest(raw/'results_manifest.json')
    derived=read(review/'review_manifest.json')
    for row in derived['files']:
        assert digest(review/row['path'])==row['sha256'],row['path'];count['cloud_review_manifest_files']+=1
    receipt=read(ROOT/'cloud_compute/rbc_full/execution_delivery_receipt.json')
    browser=read(review/'browser_final/browser_check.json');page=review/'rbc_full_review.html'
    assert browser['status']=='PASS' and not browser['fixture_only'] and browser['html_sha256']==digest(page)==receipt['html_sha256']
    assert receipt['CLOUD_RBC_RUN_COMPLETE']=='NOT_COMPLETE' and receipt['CLOUD_RBC_NUMERICAL_SCREEN']=='FAILED' and receipt['RESULTS_RETURN_VERIFIED']=='PASS'
    assert receipt['successful_returned_prep_steps'] is None and receipt['successful_returned_shear_steps']==0 and receipt['strain']==0 and receipt['qualified_speedup'] is None
    assert not (raw/'simulation/completion.json').exists()
    assets=Assets();assets.feed(page.read_text());assert not assets.assets,assets.assets
    for inp in read(raw/'provenance.json')['python_inputs']:
        relative=Path(inp['local']).relative_to('/home/lzy/projects');assert digest(ROOT/relative)==inp['sha256'],str(relative);count['frozen_python_inputs']+=1
    result=dict(status='PASS' if not sensitive and not issues else 'REVIEW_REQUIRED',sync_type='archive and read checks only; not new scientific validation',
        counts=dict(count),zero_byte_files=empty,header_only_or_empty_csv=csv_empty,nonfinite_csv_values=nonfinite,readability_issues=issues,sensitive_candidates=sensitive,
        source_check_enabled=args.check_sources,browser_and_science_tests_rerun=False,scientific_results_unchanged=True)
    if args.output:Path(args.output).write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['zero_byte_files','header_only_or_empty_csv','nonfinite_csv_values']},ensure_ascii=False,indent=2))
    return 0 if result['status']=='PASS' else 2

if __name__=='__main__':raise SystemExit(main())
