"""Read-only archive checks: hashes, modes, syntax, CSV/JSON and local links.

No solver imports, builds, browser tests, scientific tests or simulation runs.
"""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
from html.parser import HTMLParser
import json
import math
import re
from urllib.parse import unquote, urlsplit


ROOT=Path(__file__).resolve().parents[2]
REPORT=Path(__file__).resolve().parent


def digest(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    def invalid(value):raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(path.read_text(encoding='utf-8'),parse_constant=invalid)


class Assets(HTMLParser):
    def __init__(self):super().__init__();self.resources=[]
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if tag in ['script','img','iframe','source','video','audio'] and attrs.get('src'):
            self.resources.append(attrs['src'])
        if tag=='link' and attrs.get('href'):
            self.resources.append(attrs['href'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-sources',action='store_true',help='Also verify original WSL source bytes and executable bits; only on the source machine.')
    args=parser.parse_args()
    manifest=read_json(REPORT/'SYNC_MANIFEST.json')
    checks={'file_hashes':0,'executable_modes':0,'source_hashes':0,'json_readbacks':0,'csv_readbacks':0,
            'csv_data_rows':0,'python_syntax':0,'relative_entry_links':0,'empty_files':[],
            'preserved_nonfinite_csv_cells':[],'preserved_empty_csv_cells':[],'sensitive_candidates':[],'reviewed_nonsecret_matches':[]}
    reviewed=read_json(REPORT/'SENSITIVE_SCAN_REVIEW.json')
    patterns={
        'github_token':re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,255}|github_pat_[A-Za-z0-9_]{60,255})\b'),
        'private_key':re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----'),
        'aws_access_key':re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
        'openai_key':re.compile(r'\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{30,}\b'),
        'credential_url':re.compile(r'https?://[^\s/<>:"\']+:[^\s/<>@"\']+@'),
        'bearer_credential':re.compile(r'(?i)authorization["\']?\s*[:=]\s*["\']?bearer\s+[A-Za-z0-9._~-]{20,}'),
        'literal_secret_assignment':re.compile(r'''(?ix)\b(?:github_token|access_token|api_key|secret_key|password)\b ["']?\s*[:=]\s*["']([^"'\s]{12,})["']'''),
    }
    files=manifest['files']
    for entry in files:
        rel=entry['target_path'];path=ROOT/rel
        assert path.is_file() and not path.is_symlink(),rel
        assert path.stat().st_size==entry['size_bytes'] and digest(path)==entry['sha256'],rel
        checks['file_hashes']+=1
        assert bool(path.stat().st_mode&0o111)==entry['executable'],rel
        checks['executable_modes']+=1
        if args.check_sources and entry.get('source_path'):
            source=Path(entry['source_path'])
            assert source.is_file() and digest(source)==entry['sha256'],str(source)
            assert bool(source.stat().st_mode&0o111)==entry['executable'],str(source)
            checks['source_hashes']+=1
        data=path.read_bytes()
        assert not data.startswith(b'version https://git-lfs.github.com/spec/v1'),rel
        assert not ({'.git','.venv','__pycache__','node_modules','.ssh'} & set(path.relative_to(ROOT).parts)),rel
        assert not path.name.startswith('.env'),rel
        assert not path.name.endswith('.lock'),rel
        assert path.suffix not in ['.so','.a','.o','.exe','.zip'],rel
        if not data:checks['empty_files'].append(rel)
        if path.suffix=='.json':read_json(path);checks['json_readbacks']+=1
        if path.suffix=='.py':ast.parse(data.decode('utf-8'),filename=rel);checks['python_syntax']+=1
        if path.suffix=='.csv':
            with path.open(newline='') as stream:
                rows=csv.reader(stream,strict=True);header=next(rows)
                assert header and len(header)==len(set(header)),rel
                for number,row in enumerate(rows,start=2):
                    assert len(row)==len(header),(rel,number)
                    for key,value in zip(header,row):
                        if value=='':checks['preserved_empty_csv_cells'].append({'file':rel,'row':number,'column':key})
                        if value.strip().lower() in ['nan','-nan','inf','-inf','infinity','-infinity']:
                            checks['preserved_nonfinite_csv_cells'].append({'file':rel,'row':number,'column':key,'value':value})
                    checks['csv_data_rows']+=1
            checks['csv_readbacks']+=1
        if b'\0' not in data:
            text=data.decode('utf-8',errors='replace')
            for kind,pattern in patterns.items():
                for match in pattern.finditer(text):
                    # Do not print matched content, even if it later proves to be a placeholder.
                    if (kind=='literal_secret_assignment' and rel==reviewed['path']
                        and entry['sha256']==reviewed['file_sha256']
                        and hashlib.sha256(match.group(1).encode()).hexdigest() in reviewed['nonliteral_expression_sha256']):
                        checks['reviewed_nonsecret_matches'].append({'path':rel,'line':text.count('\n',0,match.start())+1,
                            'classification':'JAVASCRIPT_URL_CONCATENATION_NOT_A_LITERAL_CREDENTIAL'})
                        continue
                    checks['sensitive_candidates'].append({'path':rel,'line':text.count('\n',0,match.start())+1,'kind':kind})
    for page in [ROOT/'BENCHMARK_REVIEW.md',REPORT/'README.md']:
        for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',page.read_text()):
            link=link.strip('<>');url=urlsplit(link)
            if url.scheme or link.startswith('#'):continue
            target=page.parent/unquote(url.path)
            assert target.exists(),(str(page),link)
            checks['relative_entry_links']+=1
    final=ROOT/manifest['selected_result']['repository_directory']
    results=read_json(final/'results.json');html=(final/'benchmark_review.html').read_text()
    embedded=html.split('<script id="benchmark-audit-data" type="application/json">',1)[1].split('</script>',1)[0]
    assert json.loads(embedded)==results
    browser=read_json(final/'browser_local_profile/browser_check.json')
    assert browser['html_sha256']==digest(final/'benchmark_review.html') and browser['status']=='PASS'
    assets=Assets();assets.feed(html)
    assert not assets.resources,assets.resources
    assert results['comparison']['status']=='PARTIAL' and results['human_review']=='PENDING'
    assert len(results['comparison']['results'])==21
    assert all(x['qualified_speedup'] is None for x in results['comparison']['qualified_speedups'])
    runroot=ROOT/'mirheo_starter/runs/solver_benchmark/solver_benchmark_20260909'
    seals=0
    for path in runroot.glob('*/*/output_sha256.json'):
        for name,value in read_json(path).items():assert digest(path.parent/name)==value,(str(path),name)
        seals+=1
    assert seals==21
    for row in results['comparison']['results']:
        directory=runroot/('cpu' if row['backend']=='HemoCell' else 'gpu')/row['task_id']
        assert (directory/'task.json').exists() and (directory/'execution.json').exists()
        execution=read_json(directory/'execution.json')
        if row['status']=='COMPLETED':
            completion=read_json(directory/'completion.json')
            assert completion['actual_steps']==row['actual_steps']
            assert completion['actual_time_si']==row['actual_time_si']
            assert execution['exit_code']==0 and not execution['timeout']
            assert (directory/'profiles.csv').stat().st_size>0 and (directory/'timings.csv').stat().st_size>0
        else:assert row['task_id']=='sdpd_main' and row['actual_steps']==0
    if checks['sensitive_candidates']:
        print(json.dumps({'status':'SENSITIVE_CANDIDATE_REVIEW_REQUIRED','candidates':checks['sensitive_candidates']},ensure_ascii=False,indent=2))
        raise SystemExit(2)
    print(json.dumps({'status':'PASS_ARCHIVE_INTEGRITY_ONLY','checks':checks,'sealed_run_count':seals,
        'embedded_html_equal':True,'original_browser_evidence_matches':True,'html_external_asset_count':len(assets.resources),
        'scientific_status_preserved':'PARTIAL','human_review_preserved':'PENDING',
        'scope':'Syntax/read/hash/reference checks only. No scientific validation, solver import, simulation, compilation or browser rerun.'},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
