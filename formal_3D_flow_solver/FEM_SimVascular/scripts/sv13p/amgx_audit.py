"""Audit exact PETSc package support before an isolated compatibility build."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3p'
assert (R/'P4_SMOKE_acceptance.json').exists()
help=ROOT/'logs/sv1_3p/remote/hypre_configure_help.log';package=ROOT/'external/petsc325/petsc-3.25.5/config/BuildSystem/config/packages/amgx.py'
assert 'download-amgx' in help.read_text();source=package.read_text();assert 'amgx-2.4.0.tar.gz' in source and "maxCxxVersion    = 'c++17'" in source
out=dict(status='AUDITED',result='ISOLATED_BUILD_REQUIRED',PETSc_version='3.25.5',CUDA_version='13.2',package_version='2.4.0',source_URL='https://web.cels.anl.gov/projects/petsc/download/externalpackages/amgx-2.4.0.tar.gz',system_stack_modified=False,requires32bitint=True,precision='double',CXX='C++17',GPU_arch='89',evidence=[dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [help,package]],explanation='当前 PETSc 的官方包接口支持 AmgX 2.4.0、双精度、32 位索引与 C++17。是否兼容固定 CUDA 13.2 由独立目录中的实际编译决定；不升级、降级或覆盖已有环境。')
(R/'amgx_feasibility.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
