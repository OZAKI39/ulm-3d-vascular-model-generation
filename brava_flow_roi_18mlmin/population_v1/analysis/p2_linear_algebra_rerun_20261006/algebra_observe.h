#pragma once
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/resource.h>
namespace algebra_audit {
inline int rank(){int r;MPI_Comm_rank(MPI_COMM_WORLD,&r);return r;}
inline bool mode(const char* s){const char* p=getenv("BRAVA_ALGEBRA_MODE");return p&&std::string(p)==s;}
inline void snapshot(const char* stage,Mat A=nullptr,Vec b=nullptr,int overlap=-1){
 MatInfo info{};PetscInt m=0,n=0,lo=0,hi=0;double dn=0,on=0;const char *type="NONE";
 if(A){PetscCallAbort(MPI_COMM_WORLD,MatGetInfo(A,MAT_LOCAL,&info));MatGetSize(A,&m,&n);MatGetOwnershipRange(A,&lo,&hi);MatGetType(A,&type);
  if(std::string(type).find("mpiaij")!=std::string::npos){Mat d,o;const PetscInt* g;MatMPIAIJGetSeqAIJ(A,&d,&o,&g);MatInfo di{},oi{};MatGetInfo(d,MAT_LOCAL,&di);MatGetInfo(o,MAT_LOCAL,&oi);dn=di.nz_used;on=oi.nz_used;}}
 PetscLogDouble rss=0,peak=0,allocated=0;PetscMemoryGetCurrentUsage(&rss);struct rusage ru{};getrusage(RUSAGE_SELF,&ru);peak=1024.0*ru.ru_maxrss;PetscMallocGetCurrentUsage(&allocated);
 PetscInt vb=0;const char* vt="NONE";if(b){VecGetLocalSize(b,&vb);VecGetType(b,&vt);}
 std::ofstream f("algebra_rank_"+std::to_string(rank())+".jsonl",std::ios::app);f<<std::setprecision(17)<<"{\"stage\":\""<<stage<<"\",\"rank\":"<<rank()<<",\"overlap\":"<<overlap<<",\"matrix_type\":\""<<type<<"\",\"global_rows\":"<<m<<",\"owned_rows\":"<<hi-lo<<",\"owned_lo\":"<<lo<<",\"owned_hi\":"<<hi<<",\"local_nnz\":"<<info.nz_used<<",\"allocated_nnz\":"<<info.nz_allocated<<",\"MatGetInfo_memory_bytes\":"<<info.memory<<",\"diag_nnz\":"<<dn<<",\"offdiag_nnz\":"<<on<<",\"RSS_bytes\":"<<rss<<",\"peak_RSS_bytes\":"<<peak<<",\"PetscMalloc_current_bytes\":"<<allocated<<",\"vector_type\":\""<<vt<<"\",\"owned_vector_entries\":"<<vb<<",\"native_owned_nodes\":"<<plhs.mynNo<<",\"native_local_nodes\":"<<plhs.nNo<<",\"ghost_rows\":"<<4*(plhs.nNo-plhs.mynNo)<<"}\n";
}
inline void overlaps(Mat A){
 PetscInt lo,hi;MatGetOwnershipRange(A,&lo,&hi);
 for(int level=0;level<=2;++level){
  snapshot("before_overlap",A,nullptr,level);IS is;ISCreateStride(PETSC_COMM_SELF,hi-lo,lo,1,&is);MatIncreaseOverlap(A,1,&is,level);ISSort(is);
  PetscInt nr;ISGetLocalSize(is,&nr);Mat* sub;PetscCallAbort(MPI_COMM_WORLD,MatCreateSubMatrices(A,1,&is,&is,MAT_INITIAL_MATRIX,&sub));snapshot("ASM_submatrix",sub[0],nullptr,level);
  std::ofstream f("overlap_rank_"+std::to_string(rank())+".csv",std::ios::app);f<<level<<","<<nr<<","<<nr*sizeof(PetscInt)<<"\n";f.close();
  snapshot("after_overlap",A,nullptr,level);MatDestroySubMatrices(1,&sub);ISDestroy(&is);
 }
}
inline void numeric(Mat A,Vec b,ComMod& cm){
 long long actual=0,structural=0,dummy=0;double dmin=1e100,dmax=0,off=0,rhs=0;
 for(int i=0;i<plhs.mynNo;i++) for(int component=0;component<4;component++){
  PetscInt row=4*plhs.ltg[i]+component,nc;const PetscInt* cols;const PetscScalar* values;
  PetscCallAbort(MPI_COMM_WORLD,MatGetRow(A,row,&nc,&cols,&values));structural+=nc;
  bool phantom=component==3 && cm.ltg(plhs.map[i])>=36211;
  for(int j=0;j<nc;j++){
   actual+=values[j]!=0.0;
   if(phantom){if(cols[j]==row){dmin=std::min(dmin,(double)values[j]);dmax=std::max(dmax,(double)values[j]);}else off=std::max(off,std::abs((double)values[j]));}
  }
  MatRestoreRow(A,row,&nc,&cols,&values);
  if(phantom){dummy++;PetscScalar v;VecGetValues(b,1,&row,&v);rhs=std::max(rhs,std::abs((double)v));}
 }
 std::ofstream o("numeric_rank_"+std::to_string(rank())+".json");o<<std::setprecision(17)<<"{\"rank\":"<<rank()<<",\"structural_nnz\":"<<structural<<",\"actual_nonzero_entries\":"<<actual<<",\"dummy_pressure_rows\":"<<dummy<<",\"dummy_diag_min\":"<<dmin<<",\"dummy_diag_max\":"<<dmax<<",\"dummy_row_offdiag_maxabs\":"<<off<<",\"dummy_rhs_maxabs\":"<<rhs<<"}\n";
}
inline void dump_system(Mat A,Vec b){PetscViewer v;PetscViewerBinaryOpen(MPI_COMM_WORLD,"first_system.petsc",FILE_MODE_WRITE,&v);MatView(A,v);VecView(b,v);PetscViewerDestroy(&v);}
inline PetscErrorCode monitor(KSP ksp,PetscInt it,PetscReal norm,void*){
 if(it==0){Mat A;Vec b;KSPGetOperators(ksp,&A,nullptr);KSPGetRhs(ksp,&b);snapshot("after_PC_setup",A,b);
  Mat P;KSPGetOperators(ksp,nullptr,&P);snapshot("after_PC_setup_Pmat",P);
  PC pc;KSPGetPC(ksp,&pc);PetscInt ns;KSP* splits;
  PetscCall(PCFieldSplitGetSubKSP(pc,&ns,&splits));
  for(PetscInt j=0;j<ns;j++){
    Mat sA,sP;KSPGetOperators(splits[j],&sA,&sP);const char* tp;MatGetType(sP,&tp);
    PetscCheck(std::string(tp)=="mpiaij"||std::string(tp)=="seqaij",MPI_COMM_WORLD,PETSC_ERR_PLIB,"PC block is not CPU AIJ");
    snapshot(("CPU_split_"+std::to_string(j)).c_str(),sP);
    PC subpc;KSPGetPC(splits[j],&subpc);PetscInt levels;PetscCall(PCMGGetLevels(subpc,&levels));
    for(PetscInt lev=0;lev<levels;lev++){
      KSP smoother;PetscCall(PCMGGetSmoother(subpc,lev,&smoother));Mat m;KSPGetOperators(smoother,nullptr,&m);MatGetType(m,&tp);
      PetscCheck(std::string(tp)=="mpiaij"||std::string(tp)=="seqaij",MPI_COMM_WORLD,PETSC_ERR_PLIB,"AMG level is not CPU AIJ");
      snapshot(("CPU_AMG_"+std::to_string(j)+"_"+std::to_string(lev)).c_str(),m);
    }
  }
  PetscCall(PetscFree(splits));
}
 return PETSC_SUCCESS;
}
inline void finish(){MPI_Barrier(MPI_COMM_WORLD);PetscFinalize();int f;MPI_Finalized(&f);if(!f)MPI_Finalize();std::exit(0);}
}
