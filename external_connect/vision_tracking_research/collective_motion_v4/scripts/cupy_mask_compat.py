"""Minimal upstream ABI fix: cumsum index dtype must match CUDA INT pointer."""
import sys,subprocess,difflib,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from xcommon import *
p=THIRD/'openpiv-python-gpu/openpiv_gpu/gpu_process.py';old=p.read_text();needle='k_offset = cp.cumsum(mask)';replacement='k_offset = cp.cumsum(mask, dtype=self.dtype_i)'
if needle in old:
 assert old.count(needle)==1
 before=sha(p);new=old.replace(needle,replacement);p.write_text(new)
 dump(OUT/'evidence/cupy_mask_compat.json',dict(upstream_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=p.parent,text=True).strip(),before_sha256=before,after_sha256=sha(p),reason='Upstream float32 CUDA kernel uses int32 pointer, whereas CuPy bool cumsum defaults to int64; mismatched indexing caused illegal memory access. No PIV numerical algorithm changed.',diff=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='upstream/gpu_process.py',tofile='compatible/gpu_process.py'))))
else:assert replacement in old,'UNEXPECTED_SOURCE'
print('CuPy mask index ABI verified',sha(p))
