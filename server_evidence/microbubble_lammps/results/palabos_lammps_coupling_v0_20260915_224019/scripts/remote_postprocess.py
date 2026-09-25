from pathlib import Path
import subprocess,json,hashlib,shutil,sqlite3
r=Path('/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019');w=Path('/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019')
state=json.loads((r/'EXECUTION_STATE.json').read_text());assert state['status']=='PASS',state
p=r/'cases/kokkos_compatibility';trace=p/'kokkos_host_fix_trace.nsys-rep';assert trace.exists()
result=subprocess.run(['nsys','export','--type=sqlite','--output='+str(p/'kokkos_host_fix_trace.sqlite'),str(trace)],text=True,capture_output=True,timeout=90);(p/'NSYS_EXPORT.log').write_text(result.stdout+result.stderr);assert result.returncode==0
con=sqlite3.connect(p/'kokkos_host_fix_trace.sqlite')
tables=[x[0] for x in con.execute('select name from sqlite_master where type="table"')]
summary={'scope':'entire compatibility invocation, including initialization; no performance claim','sqlite_tables':tables}
for table in ['CUPTI_ACTIVITY_KIND_RUNTIME','CUPTI_ACTIVITY_KIND_MEMCPY','CUPTI_ACTIVITY_KIND_KERNEL']:
 if table not in tables:continue
 print(table,list(con.execute('pragma table_info('+table+')')))
 if table=='CUPTI_ACTIVITY_KIND_RUNTIME':
  rows=list(con.execute('select s.value,count(*),sum(a.end-a.start) from '+table+' a join StringIds s on a.nameId=s.id group by s.value order by sum(a.end-a.start) desc'))
  summary['CUDA_runtime_API']=[{'name':n,'calls':c,'total_seconds':t/1e9} for n,c,t in rows]
 elif table=='CUPTI_ACTIVITY_KIND_MEMCPY':
  rows=list(con.execute('select copyKind,count(*),sum(bytes),sum(end-start) from '+table+' group by copyKind'))
  summary['CUDA_memcpy']=[{'CUPTI_copyKind':k,'calls':c,'bytes':b,'seconds':t/1e9} for k,c,b,t in rows]
 else:
  row=con.execute('select count(*),sum(end-start) from '+table).fetchone();summary['CUDA_kernel_calls']=row[0];summary['CUDA_kernel_summed_seconds']=row[1]/1e9
summary['GPU_PERFORMANCE_READY']='NO';summary['COUPLING_GPU_ACCELERATION']='PENDING'
(p/'HOST_DEVICE_SYNC_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')
(r/'binaries').mkdir(exist_ok=True)
for name in ['libfrozen_flow_coupling.so','flow_audit_cli','coupled_lmp_cpu','coupled_lmp_gpu']:shutil.copy2(w/'build'/name,r/'binaries'/name)
(r/'binaries/README.md').write_text('Built on this Vast Ubuntu 24.04 instance against unchanged validated LAMMPS static libraries. SHA256 identity supplied. These executables retain build-tree runtime library paths and are evidence, not a portable binary distribution. Rebuild from CMakeLists.txt with ENGINE_ROOT for reuse; or provide compatible libfrozen_flow_coupling.so, HDF5/OpenMPI/CUDA runtime via loader paths.\n')
print('FILE_COUNT',sum(f.is_file() for f in r.rglob('*')));print('SIZE_GIB',sum(f.stat().st_size for f in r.rglob('*') if f.is_file())/2**30)
print(json.dumps({k:v for k,v in summary.items() if k not in ['sqlite_tables','CUDA_runtime_API']},indent=2));print('SYNC', [a for a in summary.get('CUDA_runtime_API',[]) if any(x in a['name'] for x in ['Synchronize','Memcpy'])])
