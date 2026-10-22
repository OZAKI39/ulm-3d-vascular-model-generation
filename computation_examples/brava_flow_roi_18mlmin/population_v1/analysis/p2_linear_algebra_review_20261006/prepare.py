from common import *
import numpy as np, pyvista as pv,xml.etree.ElementTree as ET,shutil
z=np.load(TARGET/'runtime/PIPE_CURVED_ELEVATION.npz');X=z['points_m'];C=z['tetra10'];B=z['boundary6'];tags=z['facet_tags'];n1=int(C[:,:4].max())+1;n=len(X)
i,j=np.triu_indices(10,1);pairs=np.sort(C[:,np.stack([i,j],axis=1)].reshape(-1,2),axis=1);keys=np.unique(pairs[:,0].astype(np.int64)*n+pairs[:,1]);rows=np.r_[np.arange(n),keys//n,keys%n];cols=np.r_[np.arange(n),keys%n,keys//n];order=np.lexsort((cols,rows));rows=rows[order];cols=cols[order].astype(np.int32);ptr=np.r_[0,np.cumsum(np.bincount(rows,minlength=n))];weights=np.where(np.arange(n)<n1,4,3);mixed_offsets=np.r_[0,np.cumsum(weights)];mixed_nnz=int(np.sum(weights[rows]*weights[cols]));full_nnz=len(cols)*16
np.savez_compressed(HERE/'COMPACT_MIXED_GRAPH_REFERENCE.npz',node_rowptr=ptr,node_colind=cols,mixed_offsets=mixed_offsets,components_per_node=weights.astype(np.int8),native_to_compact_node_offsets=mixed_offsets[:-1],corner_nodes=np.int64(n1))
bcnodes=np.unique(B[np.isin(tags,[1,2])]);velbc=len(bcnodes)*3;edge=n-n1
stats=dict(status='REFERENCE_ONLY_NO_NUMBERING_CHANGE',native_rows=4*n,compact_rows=3*n+n1,native_structural_scalar_nnz=full_nnz,compact_structural_scalar_nnz=mixed_nnz,removed_rows=edge,removed_scalar_nnz=full_nnz-mixed_nnz,removed_scalar_nnz_fraction=1-mixed_nnz/full_nnz,node_graph_blocks=len(cols),compact_node_blocks=len(cols),removed_node_blocks=0,compact_uniform_4x4_blocks=None,compact_variable_block_sizes=[3,4],CSR_value_index_saving_bytes=(full_nnz-mixed_nnz)*12+edge*4,meaning='Exact restriction of native structural graph to mathematical mixed DOFs before Dirichlet elimination. Includes structurally allocated corner-pressure block; not actual numeric NNZ.',representation='CSR nodal adjacency + per-node field offsets losslessly defines every compact scalar row/column; does not allocate all scalar indices')
put('COMPACT_MIXED_GRAPH_STATS.json',stats)
put('P2_DOF_STORAGE_AUDIT.json',dict(status='MIXED_STORAGE_ROOT_CAUSE_ESTABLISHED',nodes_corner=n1,nodes_edge=edge,velocity_corner_DOF=3*n1,velocity_edge_DOF=3*edge,pressure_corner_DOF=n1,edge_pressure_storage_only_slots=edge,total_equation_rows=4*n,mathematical_rows=3*n+n1,active_unconstrained_rows=3*n+n1-velbc,prescribed_velocity_rows=velbc,dummy_pressure_constrained_rows=edge,total_constrained_rows=edge+velbc,element_slots=40,element_mathematical_DOF=34,PETSc_block_size=4,block_nnz=len(cols),scalar_structural_nnz=full_nnz,numeric_nnz='NOT_MEASURED_GRAPH_ONLY',numbering='Mesh natural IDs -> native rank local O1 -> ownership O2 -> plhs.ltg PETSc node IDs -> 4*node+component. Native P2 pressure-edge residual=0 and diagonal=1; mixed volume kernels have only 4 pressure bases.',source=['fs.cpp:394-469','petsc_impl.cpp:petsc_create_lhs,petsc_create_vecmat,petsc_set_mat,petsc_set_bc'],pressure_gauge='Full outlet Cauchy traction; no pressure Dirichlet rows'))
rt=HERE/'runtime';rt.mkdir(exist_ok=True);mesh=rt/'p2/SV_MESH';(mesh/'mesh-surfaces').mkdir(parents=True,exist_ok=True)
v=pv.UnstructuredGrid(np.column_stack([np.full(len(C),10),C]).ravel(),np.full(len(C),24,np.uint8),X);v.point_data['GlobalNodeID']=np.arange(1,n+1,dtype=np.int32);v.cell_data['GlobalElementID']=np.arange(1,len(C)+1,dtype=np.int32);v.save(mesh/'mesh-complete.mesh.vtu')
# Boundary owners identify original cell IDs; no geometry/topology modifications.
facekeys={tuple(sorted(f)):i+1 for i,c in enumerate(C[:,:4]) for f in (c[[0,1,2]],c[[0,1,3]],c[[0,2,3]],c[[1,2,3]])}
plan=read(BASE/'production_workflow/gmsh_reference/GMSH_REFERENCE_MESH_V2_PLAN.json');R=plan['radius_m'];L=plan['length_m'];mu=.00345312;U=.1
for tag,name in [(1,'WALL'),(2,'INLET'),(3,'OUTLET_01')]:
 f=B[tags==tag];ids=np.unique(f);inv=np.full(n,-1);inv[ids]=np.arange(len(ids));surf=pv.PolyData(X[ids],np.column_stack([np.full(len(f),6),inv[f]]).ravel());surf.point_data['GlobalNodeID']=ids.astype(np.int32)+1;surf.cell_data['GlobalElementID']=np.array([facekeys[tuple(sorted(x[:3]))] for x in f],np.int32);surf.save(mesh/'mesh-surfaces'/f'{name}.vtp')
 if tag==2:
  speed=2*U*(1-np.sum(X[ids,:2]**2,axis=1)/R**2)
  inlet='3 2 '+str(len(ids))+'\n0\n10000000000\n'+''.join(f'{a+1}\n0 0 {u:.17g}\n0 0 {u:.17g}\n' for a,u in zip(ids,speed))
 if tag==3:
  trac=np.zeros((len(ids),3));trac[:,:2]=-4*mu*U*X[ids,:2]/R**2;trac[:,2]=8*mu*U*(X[ids,2]-L)/R**2;surf.point_data['Traction']=trac;outlet=surf
base_options=(P1/'run/PETSC_OPTIONS.txt').read_text().strip()
opts={'baseline':base_options,'asm0':base_options.replace('-pc_asm_overlap 2','-pc_asm_overlap 0')}
import shlex
s=shlex.split(base_options);remove={'-pc_type','-pc_asm_overlap','-sub_ksp_type','-sub_pc_type','-sub_pc_factor_mat_solver_type','-pc_asm_sub_mat_type'};t=[];k=0
while k<len(s):
 if s[k] in remove:k+=2
 else:t.append(s[k]);k+=1
opts['schur']=' '.join(t)+' -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type full -pc_fieldsplit_schur_precondition selfp -fieldsplit_0_ksp_type preonly -fieldsplit_0_pc_type gamg -fieldsplit_1_ksp_type preonly -fieldsplit_1_pc_type gamg'
put('LINEAR_SOLVER_CANDIDATE_PLAN.json',dict(status='FROZEN_BEFORE_TRIALS',candidate_count=2,MPI_ranks=8,candidates=[dict(name='asm0',choice='ASM overlap0 + existing PETSc CPU LU',basis='Removes measured overlap duplication, retains robust local direct saddle-point solve; actual P2 factor safety still required',expected_memory='Less than overlap2, quantitative graph-based estimate after C',options=opts['asm0']),dict(name='schur',choice='Existing native PCFIELDSPLIT velocity/pressure IS + full Schur + selfp + GAMG on both blocks',basis='Build has GAMG, no HYPRE/external LU; selfp is PETSc A11-A10 diag(A00)^-1 A01; no invented pressure operator; avoids full LU',expected_memory='Sparse matrix products and multigrid hierarchy, measured rather than assumed',options=opts['schur'])],baseline_options=opts['baseline'],scientific_configuration='UNMODIFIED',P1_protocol=dict(preferred_initial_state='original GPU step10 BIN, with independent global-node state correspondence to step10 VTU before solving; if correspondence fails, no use of this checkpoint',controlled_segment='step10 -> step11, dt .01, same initial Y/A, residual/Jacobian and all tolerances',relative_state_tolerance=1e-8,velocity_absolute_tolerance_m_s=1e-9,pressure_absolute_tolerance_Pa=1e-6,initial_R_J='first native PETSc matrix/RHS byte-identical hashes across all runs',convergence='All KSP positive reason; nonlinear native tolerance1e-10; <=24 nonlinear iterations; candidate wall<=1800s; KSPmax2000 unchanged; total Krylov <=10x baseline',GPU='all Mat/Vec CUDA, CPU PC where implemented'),P2_safety=dict(required_headroom_fraction=.2,accounting='limit minus conservative nonreclaimable upper; raw memory.current and reclaimable cache reported separately',no_unbounded_symbolic_LU=True),compact_prototype_peak_benefit_threshold=.15))
for name in ['p2_graph','p2_system','p1_baseline','p1_asm0','p1_schur']:
 case=rt/name;run=case/'run';run.mkdir(parents=True,exist_ok=True)
 (case/'SV_MESH').symlink_to(mesh if name.startswith('p2') else P1/'SV_MESH',target_is_directory=True)
 tree=ET.parse(P1/'run/solver.xml');root=tree.getroot();gp=root.find('GeneralSimulationParameters')
 gp.find('Number_of_time_steps').text='1' if name.startswith('p2') else '11'
 for key in ['Increment_in_saving_VTK_files','Increment_in_saving_restart_files']:gp.find(key).text='1'
 if name.startswith('p2'):
  ET.SubElement(root.find('Add_equation'),'Use_taylor_hood_type_basis').text='true';(run/'inlet_exact_velocity.dat').write_text(inlet);outlet.save(run/'outlet_exact_traction.vtp')
 else:
  ET.SubElement(gp,'Simulation_initialization_file_path').text=str(P1/'run/8-procs/stFile_010.bin')
  for f in ['inlet_exact_velocity.dat','outlet_exact_traction.vtp']:shutil.copy2(P1/'run'/f,run/f)
 tree.write(run/'solver.xml',encoding='utf-8',xml_declaration=True)
 option=opts[name.removeprefix('p1_')] if name.startswith('p1') else opts['baseline']
 (run/'PETSC_OPTIONS.txt').write_text(option+'\n')
put('PREPARED_INPUT_IDENTITY.json',{str(p.relative_to(HERE)):sha(p) for p in rt.rglob('*') if p.is_file() and p.suffix in {'.vtu','.vtp','.xml','.dat','.txt'}})
print('Prepared frozen graph reference and isolated cases',stats,flush=True)
