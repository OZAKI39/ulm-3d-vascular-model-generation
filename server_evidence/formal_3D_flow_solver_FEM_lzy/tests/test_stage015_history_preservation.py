from fem3d.audit import sha256
from stage015_helpers import ROOT,REPORT,read


def test_all_frozen_history_files_and_existing_source_bytes_unchanged():
    baseline=read(REPORT/'history_baseline.json')
    now={}
    for category in ('reports','outputs','logs','inputs'):
        for stage in ('stage00','stage01','stage02'):
            for p in (ROOT/category/stage).rglob('*'):
                if p.is_file(): now[str(p.relative_to(ROOT))]={'size':p.stat().st_size,'sha256':sha256(p)}
    assert now==baseline['files']
    for p,digest in baseline['existing_source_sha256'].items(): assert sha256(ROOT/p)==digest
