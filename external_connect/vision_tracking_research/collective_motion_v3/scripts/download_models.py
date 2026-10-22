"""Bounded official weight download, source and SHA provenance; no redistribution."""
from pathlib import Path
import json,hashlib,subprocess,urllib.request,time,gdown
ROOT=Path('/workspace/f10w_vision_20261017');rows=[]
for name,file,url,file_id in [
 ('LRAFT_ST_DATASET1','lraft_student_dataset1.tar','https://drive.google.com/file/d/18tTvU8WjquJLY84GHrxrgeHQXpfZ5UQh/view','18tTvU8WjquJLY84GHrxrgeHQXpfZ5UQh'),
 ('CoTracker3_online','cotracker3_scaled_online.pth','https://huggingface.co/facebook/cotracker3/resolve/main/scaled_online.pth',None)]:
 p=ROOT/'models'/file;t=time.perf_counter()
 try:
  if not p.exists():
   if file_id:gdown.download(id=file_id,output=str(p),quiet=False,resume=True,timeout=60)
   else:
    with urllib.request.urlopen(url,timeout=60) as f,p.open('wb') as out:
     while True:
      b=f.read(2**20)
      if not b:break
      out.write(b)
  with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
  rows.append(dict(name=name,file=str(p),url=url,sha256=h,bytes=p.stat().st_size,wall_s=time.perf_counter()-t,status='DOWNLOADED_NOT_YET_VALIDATED'))
 except Exception as e:rows.append(dict(name=name,url=url,error=repr(e),status='DOWNLOAD_FAILED'))
 (ROOT/'setup/model_download.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(rows[-1],flush=True)
