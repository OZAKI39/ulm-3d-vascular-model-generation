from common import *
import subprocess,difflib
assert sha(PREVIOUS/'runtime/patched.o')==read(PREVIOUS/'COMMON_GEOMETRY_PATCH_IDENTITY.json')['objects']['runtime/patched.o']['object_sha256']
old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');rt=HERE/'runtime';rt.mkdir(exist_ok=True)
# Compile only observation plumbing; scientific object is exactly the validated patch.
original=(SOURCE/'petsc_impl.cpp').read_text();s=original.replace('PetscLogStage stages[6];','PetscLogStage stages[6];\n#include "algebra_observe.h"')
s=s.replace('    /* Create vector data structures */','    algebra_audit::snapshot("before_MatSetUp");\n    /* Create vector data structures */')
s=s.replace('    /* Create vector structure for RCS preconditioner */','    algebra_audit::snapshot("after_MatSetUp",psol[cEq].A,psol[cEq].b);\n    /* Create vector structure for RCS preconditioner */')
needle='  PetscCallAbort(MPI_COMM_WORLD, PetscOptionsGetBool(nullptr, nullptr,\n      "-sv_pc_reuse_within_timestep"'
at=s.index(needle)
s=s[:at]+'''  if(algebra_audit::mode("graph")) {
    std::vector<double> w(4*com_mod.tnNo,1.0), v(4*com_mod.tnNo,0.0);
    petsc_create_linearsystem(4,eq_num,com_mod.nEq,w.data(),v.data());
    algebra_audit::overlaps(psol[eq_num].A); algebra_audit::finish();
  }
'''+s[at:]
needle='    PetscLogStagePush(stages[4]);'
s=s.replace(needle,'''    static int audit_calls=0;
    ++audit_calls;
    algebra_audit::snapshot("before_solve",psol[cEq].A,psol[cEq].b);
    if(audit_calls==1){
      if(!algebra_audit::mode("single")) algebra_audit::dump_system(psol[cEq].A,psol[cEq].b);
      KSPMonitorSet(psol[cEq].ksp,algebra_audit::monitor,nullptr,nullptr);
    }
'''+needle)
needle='    /* Rescale solution for RCS preconditioner */'
s=s.replace(needle,'''    algebra_audit::snapshot("after_solve",psol[cEq].A,psol[cEq].b);
    if(algebra_audit::mode("single")) {
      PetscInt its;KSPConvergedReason why;PetscReal rn;
      KSPGetIterationNumber(psol[cEq].ksp,&its);KSPGetConvergedReason(psol[cEq].ksp,&why);KSPGetResidualNorm(psol[cEq].ksp,&rn);
      if(algebra_audit::rank()==0){std::ofstream out("single_result.json");out<<std::setprecision(17)<<"{\\\"iterations\\\":"<<its<<",\\\"reason\\\":"<<why<<",\\\"residual\\\":"<<rn<<"}\\n";}
      algebra_audit::finish();
    }
'''+needle)
s=s.replace('  if (reuse_within_timestep_) {','  if(algebra_audit::mode("single")) algebra_audit::numeric(psol[com_mod.cEq].A,psol[com_mod.cEq].b,com_mod);\n  if (reuse_within_timestep_) {')
# Validation-only backend split. Preserve native numbering, A, b and RA policy.
needle='    KSPSetOperators(psol[cEq].ksp, psol[cEq].A, psol[cEq].A);'
assert s.count(needle)==1
s=s.replace(needle,'''    static std::vector<Mat> audit_cpu_pmat;
    if(audit_cpu_pmat.size()<=static_cast<size_t>(cEq)) audit_cpu_pmat.resize(cEq+1,nullptr);
    Mat& cpu_pmat=audit_cpu_pmat[cEq];
    if(!cpu_pmat) {
      PetscCallAbort(MPI_COMM_WORLD,MatConvert(psol[cEq].A,MATMPIAIJ,MAT_INITIAL_MATRIX,&cpu_pmat));
      PetscCallAbort(MPI_COMM_WORLD,MatSetVecType(cpu_pmat,VECMPI));
    } else PetscCallAbort(MPI_COMM_WORLD,MatCopy(psol[cEq].A,cpu_pmat,SAME_NONZERO_PATTERN));
    PetscBool identical=PETSC_FALSE;
    PetscCallAbort(MPI_COMM_WORLD,MatEqual(psol[cEq].A,cpu_pmat,&identical));
    if(!identical) SETERRABORT(MPI_COMM_WORLD,PETSC_ERR_PLIB,"CPU Pmat differs from native CUDA Amat");
    PetscCallAbort(MPI_COMM_WORLD,PetscPrintf(MPI_COMM_WORLD,"CPU_PMAT_IDENTICAL 1\\n"));
    algebra_audit::snapshot("CPU_Pmat_before_setup",cpu_pmat);
    PetscCallAbort(MPI_COMM_WORLD,KSPSetOperators(psol[cEq].ksp,psol[cEq].A,cpu_pmat));
''')
needle='psol[cEq].b, *reuse_policy, timestep, cEq));'
assert s.count(needle)==1
s=s.replace(needle,'psol[cEq].b, *reuse_policy, timestep, cEq, cpu_pmat));')
ra=(SOURCE/'sv13q_reuse.h').read_text()
ra2=ra.replace('PetscInt step, PetscInt equation) {','PetscInt step, PetscInt equation, Mat preconditioner = nullptr) {').replace('KSPSetOperators(ksp,matrix,matrix)','KSPSetOperators(ksp,matrix,preconditioner ? preconditioner : matrix)')
assert ra2!=ra
(rt/'sv13q_reuse.h').write_text(ra2)
(HERE/'RA_BACKEND_ONLY.diff').write_text(''.join(difflib.unified_diff(ra.splitlines(True),ra2.splitlines(True),fromfile='native/sv13q_reuse.h',tofile='validation/sv13q_reuse.h')))
(rt/'petsc_impl.cpp').write_text(s)
wrapper=(SOURCE/'PetscLinearAlgebra.cpp').read_text();(rt/'PetscLinearAlgebra.cpp').write_text(wrapper)
main=(SOURCE/'main.cpp').read_text();needle='    // Create LinearAlgebra objects for each equation.'
main=main.replace(needle,'''    if(getenv("BRAVA_ALGEBRA_MODE") && std::string(getenv("BRAVA_ALGEBRA_MODE"))=="initial") {
      auto& cm=simulation->com_mod;auto& st=simulation->get_integrator().get_solutions();
      auto& y=st.old.get_velocity();auto& a=st.old.get_acceleration();
      std::ofstream out("initial_rank_"+std::to_string(cm.cm.idcm())+".csv");out<<std::setprecision(17);
      for(int i=0;i<cm.tnNo;i++){out<<cm.ltg(i);for(int d=0;d<3;d++)out<<","<<cm.x(d,i);for(int d=0;d<4;d++)out<<","<<y(d,i);for(int d=0;d<4;d++)out<<","<<a(d,i);out<<"\\n";}
      out.close();MPI_Barrier(MPI_COMM_WORLD);MPI_Finalize();return 0;
    }
'''+needle)
(rt/'main.cpp').write_text(main)
meta=dict(compile_commands=[],link_commands=[],objects={},production_modified=False)
cmd0=[a.replace(old['source_checkout'],str(RETAINED)) for a in old['compile_command']]+['-I'+str(SOURCE),'-I'+str(HERE),'-I'+str(BASE/'gpu_solver_fix/support')]
for f in ['PetscLinearAlgebra.cpp','main.cpp']:
 cmd=list(cmd0);cmd[cmd.index('-c')+1]=str(rt/f);cmd[cmd.index('-o')+1]=str(rt/(f+'.o'));meta['compile_commands'].append(cmd);put('BUILD_IDENTITY.json',meta);subprocess.run(cmd,cwd=old['link_cwd'],check=True);meta['objects'][f]=dict(source_sha256=sha(rt/f),object_sha256=sha(rt/(f+'.o')))
link=list(old['link_command'])
for i,a in enumerate(link):
 if a.endswith('/PetscLinearAlgebra.cpp.o') or a.endswith('/main.cpp.o'):link[i]=str(rt/Path(a).name)
 if a.endswith('/fluid.cpp.o'):link[i]=str(PREVIOUS/'runtime/patched.o')
link[link.index('-o')+1]=str(rt/'svmultiphysics_audit_v4');meta['link_commands'].append(link);put('BUILD_IDENTITY.json',meta);subprocess.run(link,cwd=old['link_cwd'],check=True)
meta.update(status='BUILT',binary_sha256=sha(rt/'svmultiphysics_audit_v4'),fluid_object_sha256=sha(PREVIOUS/'runtime/patched.o'),retained_petsc_impl_sha256=sha(SOURCE/'petsc_impl.cpp'))
put('BUILD_IDENTITY.json',meta)
(HERE/'CPU_PC_AND_OBSERVATION.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile='native/petsc_impl.cpp',tofile='validation/petsc_impl.cpp'))+''.join(difflib.unified_diff((SOURCE/'main.cpp').read_text().splitlines(True),main.splitlines(True),fromfile='native/main.cpp',tofile='validation/main.cpp')))
print('Built independent native binary',meta['binary_sha256'],flush=True)
