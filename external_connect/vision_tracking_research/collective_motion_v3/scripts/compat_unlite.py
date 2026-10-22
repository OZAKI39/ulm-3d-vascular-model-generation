"""Minimal CuPy API compatibility; original CUDA kernels/model unchanged."""
import sys,difflib,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
p=THIRD/'UnLiteFlowNet-PIV/src/model/correlation.py';before=p.read_text()
if '@cupy.util.memoize' not in before:raise RuntimeError('Already patched or unexpected upstream version')
new=before.replace('@cupy.util.memoize(for_each_device=True)','@cupy.memoize(for_each_device=True)').replace('cupy.cuda.compile_with_cache(strKernel).get_function(strFunction)',"cupy.RawModule(code=strKernel, backend='nvcc').get_function(strFunction)")
assert new!=before
e=OUT/'evidence/unlite_compat';e.mkdir(parents=True,exist_ok=True)
dump(e/'before.json',dict(path=str(p),sha256=sha(p),reason='CuPy 13 removed util.memoize/compile_with_cache; NVCC emits native CUDA binary compatible with host driver'))
(e/'change.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),new.splitlines(True),fromfile='official/correlation.py',tofile='compat/correlation.py')))
p.write_text(new);dump(e/'after.json',dict(sha256=sha(p),model_parameters_changed=False,cuda_kernel_source_changed=False))
