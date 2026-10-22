from common import *
import numpy as np
from scipy.sparse import csr_matrix
z=np.load(PRIOR/'COMPACT_MIXED_GRAPH_REFERENCE.npz');ptr=z['node_rowptr'];col=z['node_colind'];n=len(ptr)-1;degree=np.diff(ptr)
walks=int(np.sum(degree.astype(np.int64)**2))
g=csr_matrix((np.ones(len(col),dtype=np.bool_),col,ptr),shape=(n,n));two=g@g
nnz=int(two.nnz)
put('CPU_SCHUR_GRAPH_BUDGET.json',dict(status='MEASURED_BOOLEAN_GRAPH_PRODUCT',source_sha256=sha(PRIOR/'COMPACT_MIXED_GRAPH_REFERENCE.npz'),native_node_nnz=len(col),two_hop_node_nnz=nnz,scalar_multiply_paths_upper=3*walks,schur_structural_nnz_upper=nnz,schur_CSR_payload_upper_bytes=nnz*12+(n+1)*4,full_CSR_payload_bytes=111793312*12+(1037320+1)*4,note='Pressure Schur selfp structural graph is a subset of G squared (3 velocity components repeat same graph). CPU algorithms do not store all scalar paths. Boolean product only; no CFD assembly or solve. This bounds final fine Schur CSR, not total GAMG peak.',planning_host_job_GiB=4.64+3*111793312*12/2**30+6*(nnz*12+(n+1)*4)/2**30+3,planning_extra_copies='3 full matrices + 6 fine Schur payloads + 3GiB vectors/metadata; planning estimate, not rigorous total bound'))
print(read(HERE/'CPU_SCHUR_GRAPH_BUDGET.json'),flush=True)
