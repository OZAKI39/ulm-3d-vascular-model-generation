from common import *
import subprocess
import datetime
import time

if __name__=='__main__':
    start=time.monotonic();command=[sys.executable,'-m','pytest','-q','-p','no:cacheprovider','test_field_review.py']
    result=subprocess.run(command,cwd=HERE,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    put('TEST_RESULTS.json',dict(command=command,exit_code=result.returncode,stdout=result.stdout,elapsed_s=time.monotonic()-start,
        completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='New diagnostics only, not CFD accuracy certification'))
    print(result.stdout);raise SystemExit(result.returncode)
