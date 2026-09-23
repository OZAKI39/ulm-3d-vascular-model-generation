from stage016_helpers import *
import pytest
from fem3d.audit import sha256

def test_all_frozen_stage_file_sets_sizes_hashes_and_old_core_unchanged():
 baseline=read(REPORT/'history_baseline.json');current={}
 for category in ('reports','outputs','logs','inputs'):
  for stage in ('stage00','stage01','stage01_5','stage02'):
   for p in (ROOT/category/stage).rglob('*'):
    if p.is_file():current[str(p.relative_to(ROOT))]={'sha256':sha256(p),'size':p.stat().st_size}
 assert current==baseline['files']
 for p,h in baseline['existing_source_sha256'].items():assert sha256(ROOT/p)==h
 assert read(ROOT/'reports/stage01_5/candidate_selection.json')['status']=='FAIL'
