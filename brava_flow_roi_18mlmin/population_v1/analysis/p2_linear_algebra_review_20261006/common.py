from pathlib import Path
import json,hashlib,sys,os,importlib.util
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]
TARGET=HERE.parent/'target_bc_pipe_review_20261006'
PREVIOUS=HERE.parent/'common_curved_geometry_review_20261006'
IDENTITY=HERE.parent/'stabilization_local_conservation_review_20261006'
RETAINED=BASE/'gpu_solver_fix/source'
SOURCE=RETAINED/'Code/Source/solver'
P1=BASE/'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu'
def read(p):return json.loads(Path(p).read_text())
def put(n,v):
 p=HERE/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def managed():
 s=importlib.util.spec_from_file_location('managed',BASE/'production_workflow/scripts/managed_job.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def native_env():
 e=dict(os.environ);old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');e['LD_LIBRARY_PATH']=str(Path(old['PETSc']['path']).parent)+':/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/lib:'+e.get('LD_LIBRARY_PATH','');return e
