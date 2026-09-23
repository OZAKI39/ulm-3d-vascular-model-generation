from stage017_helpers import *
import pytest
from fem3d.audit import sha256

def test_all_frozen_historical_files_and_core_unchanged():
 frozen=read(REPORT/'history_baseline.json');actual={}
 for category in ('reports','outputs','logs','inputs'):
  for stage in ('stage00','stage01','stage01_5','stage01_6','stage02'):
   for p in (ROOT/category/stage).rglob('*'):
    if p.is_file():actual[str(p.relative_to(ROOT))]={'sha256':sha256(p),'size':p.stat().st_size}
 assert actual==frozen['files']
 for path,h in frozen['existing_source_sha256'].items():assert sha256(ROOT/path)==h
 for stage in ('stage01_5','stage01_6'):assert read(ROOT/'reports'/stage/'candidate_selection.json')['status']=='FAIL'
