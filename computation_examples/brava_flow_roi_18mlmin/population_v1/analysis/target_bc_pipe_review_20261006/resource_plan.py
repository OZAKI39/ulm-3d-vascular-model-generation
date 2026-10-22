from common import *
import numpy as np,re,importlib.util

def main():
 mesh=np.load(HERE/'runtime/PIPE_CURVED_ELEVATION.npz');C=mesh['tetra10'];nn=len(mesh['points_m']);oldnodes=int(C[:,:4].max())+1
 def graph(c,n):
  # Structural nodal graph used by native 4x4 blocked storage, regardless of
  # whether a particular mixed local coefficient evaluates to zero.
  i,j=np.triu_indices(c.shape[1],1);pairs=np.sort(c[:,np.stack([i,j],axis=1)].reshape(-1,2),axis=1)
  key=pairs[:,0].astype(np.int64)*n+pairs[:,1];return n+2*len(np.unique(key))
 p1blocks=graph(C[:,:4],oldnodes);p2blocks=graph(C,nn);ratio=p2blocks/p1blocks
 base=BASE/'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu';log=(base/'run/solver.log').read_text()
 assert p1blocks*16==7719184
 # Original log rank-zero LU view, not an all-rank maximum or P2 measurement.
 pattern=r'factor fill ratio given 5\., needed ([0-9.]+)';fill=float(re.search(pattern,log).group(1))
 lu_nnz=int(re.search(r'Factored matrix:.*?nonzeros=(\d+)',log,re.S).group(1))
 sub_nnz=int(re.search(r'linear system matrix, which is also used to construct the preconditioner:.*?nonzeros=(\d+)',log,re.S).group(1))
 identity=read(BASE/'production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json');rss=identity['resource_record']['peak_tree_RSS_GiB'];p1ind=4*oldnodes;p2ind=3*nn+oldnodes
 stat=dict(l.split() for l in Path('/sys/fs/cgroup/memory.stat').read_text().splitlines())
 kernel=int(stat.get('kernel',sum(int(stat.get(k,0)) for k in ['kernel_stack','pagetables','percpu','sock','slab'])))
 nonrec=(sum(int(stat.get(k,0)) for k in ['anon','shmem','unevictable','file_dirty','file_writeback'])+kernel)/2**30
 limit=int(Path('/sys/fs/cgroup/memory.max').read_text())/2**30
 segment=log.split('SV13Q_BEGIN logical=1',1)[1].split('  0 KSP ',1)[0]
 overlap_nnzs=[int(x) for x in re.findall(r'type: seqaijcusparse\s+rows=\d+, cols=\d+, bs=4\s+total: nonzeros=(\d+)',segment)]
 assert len(overlap_nnzs)==8,overlap_nnzs
 overlap=sum(overlap_nnzs)/(p1blocks*16)
 models=[dict(name='Empirical peak tree RSS scaled by structural NNZ',GiB=rss*ratio,assumption='Same native blocked storage and ASM/LU route; actual fill/overlap may change; planning estimate, not bound'),dict(name='Empirical peak tree RSS scaled by independent DOFs',GiB=rss*p2ind/p1ind,assumption='Even ignoring redundant edge pressure slots, homogeneous per-unknown memory; planning estimate, not bound'),dict(name='LU value/index arrays only, inherited rank0 fill, no overlap',GiB=p2blocks*16*fill*12/2**30,assumption='8-byte values + 4-byte column indices, fill21.2544 from rank0 P1; deliberately omit overlap, other matrices/vectors; conditional model, not rigorous lower bound')]
 models.append(dict(name='Inherited LU fill and measured ASM overlap duplication',GiB=p2blocks*16*fill*overlap*12/2**30,assumption='Same rank0 fill and sum of eight observed P1 ASM submatrix NNZ/global NNZ; P2 partition/fill unmeasured; arrays-only planning estimate'))
 decision_models=[m for m in models if 'no overlap' not in m['name']]
 result=dict(status='PREDICTED_MEMORY_EXCEEDS_CGROUP' if all(m['GiB']+nonrec>limit for m in decision_models) else 'RESOURCE_ESTIMATE_NEEDS_REVIEW',P1_nodes=oldnodes,P2_velocity_nodes=nn,P1_independent_DOF=p1ind,P2_independent_DOF=p2ind,P2_stored_slots=4*nn,P1_structural_blocks=p1blocks,P2_structural_blocks=p2blocks,P1_scalar_structural_nnz=16*p1blocks,P2_scalar_structural_nnz=16*p2blocks,nnz_ratio=ratio,P1_measured_peak_tree_RSS_GiB=rss,P1_rank0_LU_nnz=lu_nnz,P1_rank0_submatrix_nnz=sub_nnz,P1_rank0_fill=fill,models=models,decision_models=[m['name'] for m in decision_models],no_overlap_model_role='Sensitivity only: explicitly excludes mandatory ASM overlap2, not the unchanged route',P1_ASM_submatrix_nnz=overlap_nnzs,P1_ASM_duplication_ratio=overlap,cgroup_limit_GiB=limit,current_nonreclaimable_upper_GiB=nonrec,remaining_to_cgroup_limit_GiB=limit-nonrec,resource_gate='No OOM trial under user section9.3; no preconditioner/rank/system change authorized by this fixed-route study',CFD_launched=False,estimation_limit='No P2 factorization was measured. Models extrapolate original unchanged eight-rank PETSc CUDA ASM overlap2 CPU LU. They do not prove minimum required RAM.',baseline_log_sha256=sha(base/'run/solver.log'),baseline_identity_sha256=sha(BASE/'production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json'))
 put('P2_RESOURCE_PREFLIGHT.json',result);print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
