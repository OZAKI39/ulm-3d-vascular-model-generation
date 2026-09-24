#!/usr/bin/env python3
"""Read-only matrix rank-defect evidence. Reassemble, never factor or solve."""
import json
import re
import sys
from pathlib import Path
import numpy as np
import ufl
from dolfinx import fem
from dolfinx.fem import petsc
from mpi4py import MPI
from petsc4py import PETSc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
from fem3d.vascular import load_config,load_mesh
from fem3d.spaces import create_spaces
from fem3d.boundary import wall_condition
from fem3d.formulation import block_forms
from fem3d.vascular_rank_diagnostics import constrained_pressure_support
comm=MPI.COMM_WORLD;config,config_hash=load_config(ROOT)
assert comm.size==4
base=ROOT/'outputs/stage03/reference';assert not (base/'checkpoints/primary.npz').exists()
error_codes=re.search(r'INFOG\(1\)=(-?\d+), INFO\(2\)=(\d+)',(base/'logs/solve.stderr.txt').read_text())
assert error_codes
assert json.loads((base/'metadata/resources.json').read_text())['status']=='FAIL'
raw=np.load(ROOT/config['mesh']['source']/'mesh/volume_mesh.npz')
topology=constrained_pressure_support(raw)
domain,tags,meshmeta=load_mesh(ROOT,config,comm);spaces=create_spaces(domain);V,P,R,_=spaces
bc,_=wall_condition(V,tags);jit={'cache_dir':str(ROOT/'outputs/stage03/jit_cache')}
volume=comm.allreduce(fem.assemble_scalar(fem.form(1*ufl.dx(domain=domain),jit_options=jit)),op=MPI.SUM)
p=config['physics'];a,L,scales=block_forms(domain,tags,spaces,p['dynamic_viscosity_pa_s'],p['inlet_volume_flow_m3_s'],config['inlet_geometry']['algebraic_length_scale_m'],volume)
problem=petsc.LinearProblem(a,L,u=[fem.Function(s) for s in spaces[:3]],bcs=[bc],kind='mpi',
    petsc_options_prefix='stage03_rank_diagnostic_',jit_options=jit)
problem.A.zeroEntries();petsc.assemble_matrix(problem.A,problem.a,bcs=[bc]);problem.A.assemble()
A=problem.A;start,end=A.getOwnershipRange();offsets,cols,values=A.getValuesCSR()
# Count actual nonzero entries per row, independently of structural CSR entries.
nonzero=(values!=0).astype(np.int64);cumsum=np.r_[0,np.cumsum(nonzero)];counts=cumsum[offsets[1:]]-cumsum[offsets[:-1]]
zeros=np.flatnonzero(counts==0);nvel=V.dofmap.index_map.size_local*3;npressure=P.dofmap.index_map.size_local
coords=P.tabulate_dof_coordinates();local=[]
for row in zeros:
    j=int(row)-nvel
    assert 0<=j<npressure,'Unexpected non-pressure zero row'
    local.append({'matrix_global_row':int(start+row),'owner_rank':comm.rank,'pressure_local_dof':j,
                  'pressure_coordinates_m':coords[j].tolist(),'stored_pattern_entries':int(offsets[row+1]-offsets[row]),
                  'actual_nonzero_entries':0,'row_l1_norm':0.0})
rows=[r for part in comm.allgather(local) for r in part]
for row in rows:
    vector=A.createVecRight();image=A.createVecLeft();vector.set(0)
    if start<=row['matrix_global_row']<end:vector.setValue(row['matrix_global_row'],1)
    vector.assemblyBegin();vector.assemblyEnd();A.mult(vector,image)
    row['unit_pressure_mode_matrix_product_norm']=float(image.norm())
    assert row['unit_pressure_mode_matrix_product_norm']==0
    vector.destroy();image.destroy()
for row in rows:
    candidates=topology['uncoupled_pressure_vertices'];d=np.linalg.norm(np.array([r['coordinates_m'] for r in candidates])-row['pressure_coordinates_m'],axis=1)
    assert min(d)<=1e-18
    row['physical_source_match']=candidates[int(np.argmin(d))]
assert len(rows)==len(topology['uncoupled_pressure_vertices'])
residual=json.loads((ROOT/'inputs/stage03/residual_source_evidence.json').read_text())
residual_rows=[]
centers=domain.geometry.x[domain.geometry.dofmaps[0]].mean(axis=1);nc=domain.topology.index_map(3).size_local
from scipy.spatial import cKDTree
d,indices=cKDTree(centers[:nc]).query([r['centroid_m'] for r in residual['records']]);parts=comm.allgather((d,indices));owners=np.argmin(np.stack([v[0] for v in parts]),axis=0)
for j,owner in enumerate(owners):
    if owner!=comm.rank:continue
    i=int(indices[j]);assert d[j]<=1e-15
    residual_rows.append({'source_centroid_m':residual['records'][j]['centroid_m'],'located_centroid_m':centers[i].tolist(),
                         'displacement_m':float(d[j]),'owner_rank':comm.rank,'current_local_cell':i,
                         'source_array_index_used_for_location':False,'flow_samples':'UNAVAILABLE_SOLVER_FAILED'})
residual_rows=[r for part in comm.allgather(residual_rows) for r in part]
result={'timestamp':timestamp(),'diagnostic_status':'CONFIRMED','stage_status':'FAIL','reason':'SINGULAR_PRESSURE_SUPPORT',
        'matrix_shape':list(A.getSize()),'zero_rows':rows,'zero_row_count':len(rows),'independent_exact_null_modes':len(rows),
        'topological_evidence':topology,'residual_cell_relocation':residual_rows,
        'mumps_infog_1':int(error_codes[1]),'mumps_info_2':int(error_codes[2]),
        'mumps_interpretation_source':'https://petsc.org/release/src/mat/impls/aij/mpi/mumps/impl/imumps.c.html',
        'config_sha256':config_hash,'mesh_sha256':config['mesh']['volume_mesh_sha256'],
        'core_sha256':{f:sha256(ROOT/'src/fem3d'/f) for f in ('spaces.py','boundary.py','formulation.py','solver.py')},
        'matrix_factorized':False,'ksp_solve_called':False,'matrix_altered_for_repair':False,
        'all_three_cells_relocated':len(residual_rows)==len(residual['records']),
        'solution_available':False,'fem_core_changed':False,'geometry_modified':False}
if comm.rank==0:
    write_json(base/'qc/singularity_diagnosis.json',result);print(json.dumps(result,indent=2))
