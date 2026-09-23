#!/usr/bin/env python3
"""Assembly-only preflight using the unchanged Stage 2 spaces, BC and forms."""
import json
import os
import resource
import sys
import time
from pathlib import Path
import numpy as np
import ufl
from dolfinx import fem
from dolfinx.fem import petsc
from mpi4py import MPI
from petsc4py import PETSc
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
from fem3d.vascular import load_config,load_mesh
from fem3d.spaces import create_spaces
from fem3d.boundary import wall_condition
from fem3d.formulation import block_forms

comm=MPI.COMM_WORLD;config,config_hash=load_config(ROOT)
O=ROOT/'outputs/stage03/preflight';cache=ROOT/'outputs/stage03/jit_cache'
assert not (O/'assembly_verified.json').exists()
assert json.loads((ROOT/'outputs/stage03/mesh/qc.json').read_text())['status']=='PASS'
assert comm.size==config['solver']['mpi_ranks']
domain,tags,mesh_meta=load_mesh(ROOT,config,comm)
start=time.perf_counter();spaces=create_spaces(domain)
sizes=[s.dofmap.index_map.size_global*s.dofmap.index_map_bs for s in spaces[:3]]
bc,bc_meta=wall_condition(spaces[0],tags)
jit={'cache_dir':str(cache)}
volume=comm.allreduce(fem.assemble_scalar(fem.form(1*ufl.dx(domain=domain),jit_options=jit)),op=MPI.SUM)
p=config['physics'];a,L,scales=block_forms(domain,tags,spaces,p['dynamic_viscosity_pa_s'],p['inlet_volume_flow_m3_s'],
                                        config['inlet_geometry']['algebraic_length_scale_m'],volume)
fields=[fem.Function(s) for s in spaces[:3]]
problem=petsc.LinearProblem(a,L,u=fields,bcs=[bc],kind='mpi',petsc_options_prefix='stage03_preflight_',
                           petsc_options={'ksp_type':'preonly','pc_type':'lu','pc_factor_mat_solver_type':'mumps'},jit_options=jit)
compiled=time.perf_counter()
# Public assembly functions reproduce LinearProblem's assembly path, stopping before KSP.solve/setup.
problem.A.zeroEntries();petsc.assemble_matrix(problem.A,problem.a,bcs=[bc]);problem.A.assemble()
with problem.b.localForm() as local:local.set(0)
petsc.assemble_vector(problem.b,problem.L)
bcs1=fem.bcs_by_block(fem.extract_function_spaces(problem.a,1),[bc])
petsc.apply_lifting(problem.b,problem.a,bcs=bcs1)
problem.b.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE)
bcs0=fem.bcs_by_block(fem.extract_function_spaces(problem.L),[bc]);petsc.set_bc(problem.b,bcs0)
comm.barrier();end=time.perf_counter()
info=problem.A.getInfo(PETSc.Mat.InfoType.LOCAL)
nnz=int(comm.allreduce(info['nz_used'],op=MPI.SUM))
local={'rank':comm.rank,'pid':os.getpid(),'peak_rss_kib':int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
       'matrix_local_rows':problem.A.getLocalSize()[0],'real_owned_dofs':spaces[2].dofmap.index_map.size_local}
ranks=comm.gather(local,root=0)
assert problem.A.getSize()==(sum(sizes),sum(sizes)) and sizes[2]==1
assert problem.b.norm()>0 and nnz>0
result={'status':'PASS','timestamp':timestamp(),'mpi_ranks':comm.size,'velocity_dofs':sizes[0],
        'pressure_dofs':sizes[1],'real_dofs':sizes[2],'total_dofs':sum(sizes),'matrix_shape':list(problem.A.getSize()),
        'nnz':nnz,'matrix_type':problem.A.getType(),'rhs_l2_norm_scaled':float(problem.b.norm()),
        'spaces_forms_jit_s':comm.allreduce(compiled-start,op=MPI.MAX),
        'assembly_wall_time_s':comm.allreduce(end-compiled,op=MPI.MAX),
        'preflight_wall_time_s':comm.allreduce(end-start,op=MPI.MAX),
        'boundary_conditions':bc_meta,'mesh':mesh_meta,'algebraic_scaling':scales,
        'matrix_factorized':False,'ksp_solve_called':False,'gpu_used':False,'config_sha256':config_hash,
        'core_sha256':{f:sha256(ROOT/'src/fem3d'/f) for f in ('spaces.py','boundary.py','formulation.py','solver.py')}}
if comm.rank==0:
    result['ranks']=ranks;result['max_rank_peak_rss_kib']=max(r['peak_rss_kib'] for r in ranks)
    result['sum_rank_peak_rss_kib']=sum(r['peak_rss_kib'] for r in ranks)
    result['nnz_reduction']='Explicit LOCAL nz_used followed by MPI SUM; supersedes default-global-sum diagnostic in assembly.json'
    write_json(O/'assembly_verified.json',result);print(json.dumps(result,indent=2),flush=True)
