#!/usr/bin/env python3
"""Record the failed formal attempt without constructing surrogate flow fields."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
from fem3d.vascular import load_config,reynolds
base=ROOT/'outputs/stage03/reference';R=ROOT/'reports/stage03'
read=lambda p:json.loads(p.read_text())
config,h=load_config(ROOT);diag=read(base/'qc/singularity_diagnosis.json');resources=read(base/'metadata/resources.json')
assert diag['stage_status']=='FAIL' and diag['zero_row_count']==2 and not (base/'checkpoints/primary.npz').exists()
for directory in ('solution','qc','metadata','checkpoints','logs'):(base/directory).mkdir(parents=True,exist_ok=True)
run=next(p for p in (ROOT/'logs/stage03').glob('*reference_vascular_mumps_r4/metadata.json'))
metadata=read(run)
samples=[json.loads(line) for line in (base/'logs/resources.jsonl').read_text().splitlines()]
accounting={'timestamp':timestamp(),'valid_peak_rss_kib':resources['rusage_children_peak_rss_kib'],
            'valid_peak_scope':'Linux RUSAGE_CHILDREN ru_maxrss: maximum individual child peak RSS, not simultaneous MPI sum',
            'elapsed_time_s':resources['elapsed_time_s'],
            'minimum_sampled_available_memory_bytes':min(row['effective_available_bytes'] for row in samples),
            'peak_sampled_cgroup_memory_bytes':max(int(row['cgroup_memory_current']) for row in samples),
            'cgroup_scope':'Whole container including pre-existing processes and file cache; not FEM-only RSS',
            'sampled_process_group_rss_valid_for_mpi_total':False,
            'sampling_limitation':'MPICH workers use separate process groups, so the initial process-group sampler omitted worker RSS. Its approximately 5.5 MiB measurement is excluded from reported solver peak. No second PDE factorization was performed to replace it.',
            'raw_resource_record_sha256':sha256(base/'metadata/resources.json')}
write_json(base/'metadata/resource_accounting.json',accounting)
failure={'status':'FAIL','reason':'SINGULAR_PRESSURE_SUPPORT','timestamp':timestamp(),
         'condition_type':config['condition_type'],'experimental':False,'experimental_run':'NOT_EXECUTED',
         'warning':config['warning'],'formal_factorization_attempts':1,'pde_solution_available':False,
         'mumps_infog_1':diag['mumps_infog_1'],'mumps_info_2':diag['mumps_info_2'],
         'converged':False,'factorization_status':'FAILED_NUMERICALLY_SINGULAR',
         'converged_reason':None,'converged_reason_note':'Core raised from KSPSetUp/MatLUFactorNumeric before returning; no converged reason or residual was fabricated.',
         'solution_residual':None,'completed_solve_time_s':None,'attempt_elapsed_s':resources['elapsed_time_s'],
         'config_sha256':h,'mesh_sha256':config['mesh']['volume_mesh_sha256'],
         'core_sha256':diag['core_sha256'],'hostname':metadata['hostname'],'mpi_ranks':4,'gpu_used':False,
         'formal_run_record':str(run.relative_to(ROOT)),'formal_run_record_sha256':sha256(run),
         'formal_run_source_hashes':{name:value for name,value in metadata['input_sha256'].items() if '/src/' in name or '/scripts/' in name},
         'assembly_evidence':'outputs/stage03/preflight/assembly_verified.json',
         'failure_evidence':'outputs/stage03/reference/qc/singularity_diagnosis.json',
         'resource_blocked':False,'pressure_pins_added':False,'stabilization_added':False,
         'geometry_changed':False,'iterative_solver_used':False,'formal_retry_performed':False,
         'roundtrip':'NOT_EXECUTED_NO_SOLUTION','solution_finiteness':'NOT_EVALUATED_NO_SOLUTION',
         'flux_gates':'NOT_EVALUATED_NO_SOLUTION','residual_flow_monitor':'NOT_EVALUATED_NO_SOLUTION',
         'residual_cell_relocation':'PASS','residual_cells_relocated':len(diag['residual_cell_relocation']),
         'reference_reynolds':reynolds(config)}
assert failure['formal_run_source_hashes'],'Original run code hashes must be present'
write_json(base/'metadata/failure.json',failure)
write_json(base/'qc/unavailable.json',{'status':'NOT_EVALUATED_NO_SOLUTION','reason':failure['reason'],
    'unavailable':['velocity','pressure','lambda','flux_integrals','mass_closure','wall_velocity','divergence','derived_fields','local_flow_monitor','solution_roundtrip'],
    'no_target_flow_substitution':True,'no_zero_filled_fields':True})
for name,path in [('failure_summary.json',base/'metadata/failure.json'),('singularity_diagnosis.json',base/'qc/singularity_diagnosis.json'),('resource_accounting.json',base/'metadata/resource_accounting.json')]:
    write_json(R/name,read(path))
print(json.dumps({k:failure[k] for k in ('status','reason','formal_factorization_attempts','reference_reynolds')},indent=2))
