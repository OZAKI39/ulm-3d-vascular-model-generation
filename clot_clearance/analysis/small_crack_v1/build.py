from common import *
import shlex,subprocess,time,difflib,resource
assert not LIB.exists()
oldsrc=NATIVE_PREV/'source/src/fix_move.cpp';src=O/'source_identity/fix_move.cpp';old=oldsrc.read_text();start=old.index('  // Isolated event integrator:',old.index('void FixMove::reset_dt()'))
new=old[:start]+'''  // Isolated mixed-DOF extension: each Cartesian component must be either
  // a pair of equal-style cumulative-time expressions or two NULLs.
  // NULL components retain the existing native Verlet kick/drift/kick.
  // Displacement-only, velocity-only and atom-style expressions are rejected.
  auto admitted = [&](char *position, char *velocity, const int &ps, const int &vs) {
    return (!position && !velocity) ||
           (position && velocity && ps == EQUAL && vs == EQUAL);
  };
  if (mstyle != VARIABLE || !admitted(xvarstr,vxvarstr,xvarstyle,vxvarstyle) ||
      !admitted(yvarstr,vyvarstr,yvarstyle,vyvarstyle) ||
      !admitted(zvarstr,vzvarstr,zvarstyle,vzvarstyle))
    error->all(FLERR, Error::NOLASTLINE,
               "Variable dt requires paired equal variables or paired NULLs in isolated fix move");
  dt = dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
'''
src.write_text(new);(O/'FIX_MOVE_MIXED_DOF.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='validated/fix_move.cpp',tofile='isolated_small_crack/fix_move.cpp')))
oldbuild=NATIVE_PREV/'build';flags={}
for line in (oldbuild/'CMakeFiles/lammps.dir/flags.make').read_text().splitlines():
 if line.startswith('CXX_'):
  key,value=line.split('=',1);flags[key.strip()]=shlex.split(value)
obj=O/'build/fix_move.cpp.o';compile=['/usr/bin/c++',*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],'-c',str(src),'-o',str(obj)]
link=shlex.split((oldbuild/'CMakeFiles/lammps.dir/link.txt').read_text());link[link.index('-o')+1]=str(O/'build/liblammps.so.0');objects=[]
for i,tok in enumerate(link):
 if tok.startswith('@'):
  objs=shlex.split((oldbuild/tok[1:]).read_text());assert sum(z.endswith('/fix_move.cpp.o') for z in objs)==1
  objects=[str(obj) if z.endswith('/fix_move.cpp.o') else str((oldbuild/z).resolve()) for z in objs]
  rsp=O/'build/objects.rsp';rsp.write_text(' '.join(shlex.quote(z) for z in objects));link[i]='@'+str(rsp)
start=time.perf_counter();cpu=resource.getrusage(resource.RUSAGE_CHILDREN)
with (O/'build/BUILD.log').open('w') as log:
 subprocess.run(compile,check=True,stdout=log,stderr=subprocess.STDOUT,cwd=oldbuild)
 subprocess.run(link,check=True,stdout=log,stderr=subprocess.STDOUT,cwd=oldbuild)
LIB.symlink_to('liblammps.so.0');end=resource.getrusage(resource.RUSAGE_CHILDREN)
dump(O/'BUILD_AUDIT.json',dict(status='PASS',compiled_translation_units=['fix_move.cpp'],reused_immutable_object_count=len(objects)-1,compile_command=compile,link_command=link,commands_cwd=str(oldbuild),old_objects_only_read=True,wall_s=time.perf_counter()-start,CPU_s=end.ru_utime+end.ru_stime-cpu.ru_utime-cpu.ru_stime,peak_RSS_bytes=end.ru_maxrss*1024,object_sha256={p:sha(p) for p in objects}))
for p,h in json.loads((O/'PROTECTED_BASELINES.json').read_text()).items():assert sha(p)==h,p
dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',dict(status='PASS',upstream='8de817dd79bfe4525d5d39246a212d833e6dee07',previous_library=str(OLDLIB),previous_library_sha256=sha(OLDLIB),actual_isolated_library=str(LIB),library_sha256=sha(LIB),changed_CPP=['fix_move.cpp reset_dt admission only'],material_damage_contact_Verlet_sources_unchanged=True,source_sha256={p.name:sha(p) for p in (O/'source_identity').iterdir()},production_environment_unchanged=True,build_inputs='BUILD_AUDIT.json'))
print('BUILD_PASS',sha(LIB),time.perf_counter()-start)
