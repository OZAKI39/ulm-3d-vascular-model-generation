from common import *
import numpy as np,xml.etree.ElementTree as ET
checks={}
for n,h in read(HERE/'INPUT_SHA256.json').items():checks['frozen:'+n]=sha(BASE/n)==h
ref=np.load(HERE/'COMPACT_MIXED_GRAPH_REFERENCE.npz');rp=ref['node_rowptr'];ci=ref['node_colind'];w=ref['components_per_node'];offset=ref['mixed_offsets'];n=len(w)
checks.update(reference_offsets=np.array_equal(np.diff(offset),w),reference_rows=int(offset[-1])==814201,reference_native_graph=len(ci)*16==111793312,reference_CSR_consistent=len(rp)==n+1 and rp[0]==0 and rp[-1]==len(ci) and bool(np.all(np.diff(rp)>0)),reference_cols_in_range=bool(np.all((ci>=0)&(ci<n))))
# Element-wise mathematical mapping checks independently use frozen connectivity.
z=np.load(TARGET/'runtime/PIPE_CURVED_ELEVATION.npz');C=z['tetra10'];nc=int(C[:,:4].max())+1
checks['all_corner_components4']=bool(np.all(w[:nc]==4));checks['all_midpoint_components3']=bool(np.all(w[nc:]==3));checks['all_elements_34_mathematical_slots']=bool(np.all(w[C].sum(1)==34))
checks['compact_reference_not_production_numbering']=read(HERE/'COMPACT_MIXED_GRAPH_STATS.json')['status']=='REFERENCE_ONLY_NO_NUMBERING_CHANGE'
for name in ['baseline','asm0','schur']:
 r=ET.parse(HERE/f'runtime/p1_{name}/run/solver.xml').getroot();old=ET.parse(P1/'run/solver.xml').getroot()
 checks[f'{name}:equation_byte_equivalent']=ET.tostring(r.find('Add_equation'))==ET.tostring(old.find('Add_equation'))
 for k in ['Time_step_size','Spectral_radius_of_infinite_time_step']:
  checks[f'{name}:{k}']=r.find('GeneralSimulationParameters/'+k).text==old.find('GeneralSimulationParameters/'+k).text
 checks[f'{name}:same_checkpoint']=r.find('GeneralSimulationParameters/Simulation_initialization_file_path').text==str(P1/'run/8-procs/stFile_010.bin')
 checks[f'{name}:exact_inlet']=sha(HERE/f'runtime/p1_{name}/run/inlet_exact_velocity.dat')==sha(P1/'run/inlet_exact_velocity.dat')
 checks[f'{name}:exact_traction']=sha(HERE/f'runtime/p1_{name}/run/outlet_exact_traction.vtp')==sha(P1/'run/outlet_exact_traction.vtp')
graph=read(HERE/'P2_PETSC_GRAPH_MEASUREMENTS.json');checks['native_global_rows']=sum(r['owned_rows'] for r in graph['rank_data'])==4*n;checks['diag_offdiag_sum']=graph['scalar_diag_nnz']+graph['scalar_offdiag_nnz']==111793312
checks['native_MPI_ownership_contiguous']=sorted((r['owned_lo'],r['owned_hi']) for r in graph['rank_data'])==[(a,b) for a,b in zip(np.r_[0,np.cumsum([r['owned_rows'] for r in sorted(graph['rank_data'],key=lambda v:v['rank'])])[:-1]],np.cumsum([r['owned_rows'] for r in sorted(graph['rank_data'],key=lambda v:v['rank'])]))]
checks['validated_fluid_object_unchanged']=sha(PREVIOUS/'runtime/patched.o')==read(PREVIOUS/'COMMON_GEOMETRY_PATCH_IDENTITY.json')['objects']['runtime/patched.o']['object_sha256']
checks['only_two_candidates']=read(HERE/'LINEAR_SOLVER_CANDIDATE_PLAN.json')['candidate_count']==2
oldid=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');checks['production_binary_unchanged']=sha(oldid['actual_binary'])==oldid['actual_binary_sha256'];checks['PETSc_library_unchanged']=sha(oldid['PETSc']['path'])==oldid['PETSc']['sha256']
from gpu_budget import gate
checks['GPU_guard_no_budget_refuses']=not gate(None,15477)['PASS']
budget=dict(status='BOUNDED',MPI_ranks=8,source='native_workspace_query_without_numeric_allocation',sum_rank_simultaneous_workspace_bound_MiB=12000,resident_increment_bound_MiB=3000)
checks['GPU_guard_insufficient_free_refuses']=not gate(budget,15477)['PASS']
checks['GPU_guard_sufficient_proven_budget_passes']=gate(budget,20000)['PASS']
checks['GPU_guard_mismatched_ranks_refuses']=not gate(budget,20000,ranks=4)['PASS']
put('TEST_RESULTS.json',dict(status='PASS' if all(checks.values()) else 'FAIL',checks={k:bool(v) for k,v in checks.items()},passed=sum(bool(v) for v in checks.values()),failed=[k for k,v in checks.items() if not v],scope='Engineering identity/graph/config checks, not CFD field validation'))
assert all(checks.values()),[k for k,v in checks.items() if not v]
print('Engineering checks PASS',len(checks),flush=True)
