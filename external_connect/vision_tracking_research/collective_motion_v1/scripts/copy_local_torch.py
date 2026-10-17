"""Create independent dependency copies from installed distributions; source read-only.
Avoid redownloading several GB of the same CUDA wheels. Pip RECORD provenance saved.
"""
import sys,shutil,importlib.metadata as md,json,hashlib
from pathlib import Path
from packaging.requirements import Requirement
src=Path(__file__).resolve().parents[2]/'smart_online_tracking_v2/.venv_annotation/lib/python3.12/site-packages'
root=Path(__file__).resolve().parents[1];dst=root/'.venv/lib/python3.12/site-packages'
available={d.metadata['Name'].lower().replace('_','-'):d for d in md.distributions(path=[str(src)])}
installed={d.metadata['Name'].lower().replace('_','-') for d in md.distributions(path=[str(dst)])}
queue=['torch','torchvision','huggingface-hub','safetensors'];seen=set();copied=[]
while queue:
    name=queue.pop().lower().replace('_','-')
    if name in seen or name in installed:continue
    seen.add(name);d=available[name]
    for requirement in d.requires or []:
        r=Requirement(requirement)
        if r.marker and not r.marker.evaluate():continue
        queue.append(r.name)
    # cuda-toolkit extras carry the actual shared-library distributions.
    if name=='cuda-toolkit':queue.extend(k for k in available if k.startswith('nvidia-'))
    count=0
    for f in d.files or []:
        if '..' in f.parts:continue # console entry points not needed for inference
        a=src/f;b=dst/f
        if not a.is_file():continue
        b.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(a,b);count+=1
    record=d.read_text('RECORD') or ''
    copied.append({'name':name,'version':d.version,'copied_files':count,'source_RECORD_sha256':hashlib.sha256(record.encode()).hexdigest()})
    print(name,d.version,count,flush=True)
(root/'evidence/torch_dependency_copy.json').write_text(json.dumps({'source':str(src),'destination':str(dst),'mode':'independent copies, not symlinks; original unmodified','distributions':copied},indent=2))
