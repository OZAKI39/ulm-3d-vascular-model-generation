#!/usr/bin/env python3
"""One formal solve, or independent checkpoint diagnostics; never refactor on reload."""
import argparse
import json
import os
import platform
import resource
import sys
import time
from pathlib import Path
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from fem3d.audit import sha256,timestamp,write_json
from fem3d.vascular import load_config,load_mesh,save_primary,restore_primary
from fem3d.vascular_diagnostics import diagnose,export_visualization
from fem3d.solver import solve_stokes

parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['solve','finish','reload']);args=parser.parse_args()
comm=MPI.COMM_WORLD;config,config_hash=load_config(ROOT);base=ROOT/'outputs/stage03/reference'
cache=ROOT/'outputs/stage03/jit_cache'
started=time.perf_counter()
if args.mode=='solve':
    assert comm.size==4 and os.environ.get('OMP_NUM_THREADS')=='1'
    assert not (base/'checkpoints/primary.npz').exists()
    preflight=json.loads((ROOT/'outputs/stage03/preflight/assembly_verified.json').read_text())
    assert preflight['status']=='PASS' and preflight['config_sha256']==config_hash
    assert json.loads((base/'metadata/resource_gate.json').read_text())['status']=='PASS'
    domain,tags,mesh_meta=load_mesh(ROOT,config,comm)
    PETSc.Log.begin()
    physics=config['physics']
    print(f'Rank {comm.rank}: entering the unique validated direct solve',flush=True)
    fields=solve_stokes(domain,tags,physics['dynamic_viscosity_pa_s'],physics['inlet_volume_flow_m3_s'],
                        config['inlet_geometry']['algebraic_length_scale_m'],cache)
    checkpoint_hash=save_primary(base,fields)
    events={}
    for name in ('KSPSolve','PCSetUp','MatLUFactorSym','MatLUFactorNum','MatSolve'):
        try:events[name]={k:float(v) for k,v in PETSc.Log.Event(name).getPerfInfo().items()}
        except (AttributeError,KeyError,TypeError,PETSc.Error) as exc:events[name]={'unavailable':str(exc)}
    all_events=comm.gather(events,root=0)
    ranks=comm.gather({'rank':comm.rank,'pid':os.getpid(),'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},root=0)
    if comm.rank==0:
        event_times={name:max(row[name].get('time',0) for row in all_events) for name in events}
        metadata={'status':'PASS','timestamp':timestamp(),'hostname':platform.node(),'pid':os.getpid(),
                  'config_sha256':config_hash,'checkpoint_sha256':checkpoint_hash,
                  'mesh':mesh_meta,'solver':fields['solver'],'factorization_status':'SUCCESS',
                  'source_manifest_sha256':sha256(ROOT/'remote/source_manifest.json'),
                  'code_sha256':{str(p.relative_to(ROOT)):sha256(p) for p in sorted((ROOT/'src/fem3d').glob('*.py'))},
                  'runner_sha256':sha256(Path(__file__)),'reference_pde_solve_count':1,
                  'petsc_events_per_rank':all_events,'petsc_event_max_wall_time_s':event_times,
                  'solve_time_s':event_times['KSPSolve'] or None,'ranks_at_checkpoint':ranks,
                  'assembly_preflight':preflight,'experimental':False,'gpu_used':False,
                  'warning':'NOT_EXPERIMENTAL_PUMP_FLOW','unique_formal_condition':'stage03_reference_vascular'}
        write_json(base/'metadata/solve.json',metadata)
        print('Primary P2/P1 coefficients and Real multiplier saved before diagnostics.',flush=True)
    comm.barrier()
else:
    domain,tags,fields,mesh_meta,distance=restore_primary(ROOT,config,comm)

qc,samples=diagnose(ROOT,config,domain,tags,fields)
if args.mode in ('solve','finish'):
    if comm.rank==0:
        write_json(base/'qc/solution.json',qc)
        write_json(base/'qc/boundary_fluxes.json',qc['flux'])
        write_json(base/'qc/residual_cells.json',qc['residual_cells'])
    export_visualization(base,fields,samples)
else:
    original=json.loads((base/'qc/solution.json').read_text())
    compare={}
    def measure(name,old,new,scale):
        error=float(np.max(np.abs(np.asarray(new)-np.asarray(old))))
        compare[name]={'absolute_error':error,'normalization':float(scale),'relative_error':error/max(float(scale),1e-300),
                       'pass':bool(error<=1e-10*max(float(scale),1e-300))}
    Q=config['physics']['inlet_volume_flow_m3_s']
    measure('inlet',original['flux']['Q_in_m3_s'],qc['flux']['Q_in_m3_s'],Q)
    for outlet in ('outlet_01','outlet_02','outlet_03'):
        measure(outlet,original['flux']['Q_outlets_signed_m3_s'][outlet],qc['flux']['Q_outlets_signed_m3_s'][outlet],Q)
    measure('mass_closure_dimensionless',original['flux']['relative_mass_closure'],qc['flux']['relative_mass_closure'],1)
    measure('lambda',original['lambda_pa'],qc['lambda_pa'],abs(original['lambda_pa']))
    for name in ('velocity_L2_norm_m_pow_2p5_s','pressure_integral_pa_m3'):
        measure(name,original['norms'][name],qc['norms'][name],abs(original['norms'][name]))
    for j,(old,new) in enumerate(zip(original['residual_cells']['cells'],qc['residual_cells']['cells'])):
        assert old['source_centroid_m']==new['source_centroid_m']
        for key,value in old['values'].items():
            measure(f'residual_{j}_{key}',value,new['values'][key],np.linalg.norm(value))
        for key,stats in old['neighbors'].items():
            for kind in ('min','median','max'):
                measure(f'residual_{j}_neighbor_{key}_{kind}',stats[kind],new['neighbors'][key][kind],np.linalg.norm(stats[kind]))
    saved_meta=json.loads((base/'metadata/solve.json').read_text())
    result={'status':'PASS' if qc['status']=='PASS' and all(r['pass'] for r in compare.values()) else 'FAIL',
            'timestamp':timestamp(),'new_process':True,'pid':os.getpid(),'solve_pid':saved_meta['pid'],
            'mpi_ranks':comm.size,'pde_solve_called':False,'checkpoint_sha256':saved_meta['checkpoint_sha256'],
            'exact_owned_coefficients_restored':True,'max_coordinate_match_distance_m':distance,
            'comparisons':compare,'recomputed_qc':qc,'wall_time_s':time.perf_counter()-started}
    if comm.rank==0:write_json(base/'qc/roundtrip.json',result)
    assert result['status']=='PASS'
if comm.rank==0:
    print(json.dumps({'mode':args.mode,'status':qc['status'],'flux':qc['flux'],'lambda_pa':qc['lambda_pa'],
                      'Re':qc['reynolds']['Re'],'wall_time_s':time.perf_counter()-started},indent=2),flush=True)
assert qc['status']=='PASS'
