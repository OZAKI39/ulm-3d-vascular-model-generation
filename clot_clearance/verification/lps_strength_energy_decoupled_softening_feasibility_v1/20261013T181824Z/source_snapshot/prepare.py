from common import *
import os,sys,subprocess,shutil,time,resource,scipy
start=time.perf_counter();O.mkdir(parents=True,exist_ok=False);F.mkdir(parents=True,exist_ok=False);(O/"provenance").mkdir()
def cmd(args):
 z=subprocess.run(args,cwd=R,text=True,capture_output=True);return dict(returncode=z.returncode,stdout=z.stdout,stderr=z.stderr)
prior=json.loads((R/"verification/lps_mode_i_low_load_admission_and_solver_recovery_v1/20261013T122000Z/SOURCE_AND_LIBRARY_IDENTITY.json").read_text())
paths=list(prior["files"])+[str(R/"analysis/timelayer_v1/model.py"),str(R/"analysis/native_state_v1/audits.py")]
files={s:dict(sha256=sha(s),bytes=Path(s).stat().st_size) for s in paths}
assert files[str(OLD/"build_owner_fix/liblammps.so")]["sha256"]==EXPECTED
assert files[str(R/"tools/lammps-native-y-wall-shared-v1-build/liblammps.so")]["sha256"]=="d1b2f53d6b942d1ca9277a8b8cd9c195ce75ea6f616681ec21c311e4faf3b41a"
identity=dict(files=files,pwd=str(R),git_status=cmd(["git","status","--short"]),git_branch=cmd(["git","branch","--show-current"]),git_head=cmd(["git","rev-parse","HEAD"]),python=sys.version,executable=sys.executable,numpy=np.__version__,scipy=scipy.__version__,affinity=sorted(os.sched_getaffinity(0)),cgroup={str(p):p.read_text().strip() for p in [Path("/sys/fs/cgroup/cpu,cpuacct/cpu.cfs_quota_us"),Path("/sys/fs/cgroup/cpu,cpuacct/cpu.cfs_period_us"),Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")]},memory=cmd(["free","-b"]),disk=cmd(["df","-B1","/workspace"]),jobs=cmd(["ps","-eo","pid,comm,pcpu,pmem,etime","--sort=-pcpu"]),source_upstream="8de817dd79bfe4525d5d39246a212d833e6dee07",native_cpp_changed=False)
dump(O/"SOURCE_AND_LIBRARY_IDENTITY.json",identity)
for name in ["reference.npz","body.data","T0.npz","postbreak_same_coordinates.npz","report.json"]:shutil.copyfile(OLD/"final_fracture"/name,O/"provenance"/name)
protected={}
for base in [R/"verification",R/"results"]:
 for p in base.rglob("*"):
  if p.is_file() and C["task"] not in p.parts:protected[str(p.relative_to(R))]=dict(sha256=sha(p),bytes=p.stat().st_size)
dump(O/"HISTORICAL_PROTECTION_BEFORE.json",protected)
z=np.load(O/"provenance/reference.npz")
assert len(z["X"])==341 and len(z["bonds"])==12047 and np.all(z["V"]==DX**3)
dump(O/"SPECIMEN_IDENTITY.json",dict(nodes=341,bonds=12047,reference_sha256=sha(O/"provenance/reference.npz"),source=str(OLD/"final_fracture"),V_m3=float(z["V"][0]),coordinate_min_m=z["X"].min(0),coordinate_max_m=z["X"].max(0),id_order="original immutable 1-based atoms, original zero-based edge endpoint order",m_reference_sha256=hashlib.sha256(z["m"].tobytes()).hexdigest()))
dump(O/"PREPARATION_RESOURCES.json",dict(wall_s=time.perf_counter()-start,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,protected_files=len(protected)))
print("PREPARED",len(protected),flush=True)
