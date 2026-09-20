import subprocess
import sys
from fem_sync_support import ROOT


def test_fresh_process_reads_every_frozen_artifact_from_checkout(tmp_path):
    result=subprocess.run([sys.executable,'-I','-B',str(ROOT/'scripts/fem_freeze_sync/validate_frozen.py')],
        cwd=tmp_path,capture_output=True,text=True,timeout=90)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'PASS: fresh process' in result.stdout
