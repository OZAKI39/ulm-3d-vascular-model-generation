"""Independent cloud A/B checks using the archived DPD parameters and native OpenMPI helper."""
from pathlib import Path
import sys,json,ctypes as C,hashlib,time,csv,gc

def main():
    cfg=json.loads(Path(sys.argv[1]).read_text());mode=sys.argv[2];sys.path.insert(0,cfg['cloud_project_root'])
    from py_scripts.solver_benchmark.native_mpi import NativeMPI
    mpi=NativeMPI();mpi.check(mpi.lib.MPI_Init(None,None));compute=mpi.rank==0
    import mirheo as mir
    import libmirheo
    binary=Path(libmirheo.__file__).resolve();assert str(binary)==cfg['runtime_binary_path']
    assert hashlib.sha256(binary.read_bytes()).hexdigest()==cfg['runtime_binary_hash']
    assert str(Path(mir.__file__).resolve())==cfg['python_package']+'/__init__.py'
    assert str(binary) in Path('/proc/self/maps').read_text(),'BINARY_NOT_ACTUALLY_MAPPED'
    runtime=C.CDLL('libcudart.so.12') if compute else None
    def sync():
        if compute:assert runtime.cudaDeviceSynchronize()==0,'CUDA_SYNCHRONIZATION_FAILED'
    sync();mpi.Barrier()
    if mode=='A':
        if compute:Path('loader.json').write_text(json.dumps(dict(status='CLOUD_LOADER_PASS',mirheo_module=mir.__file__,runtime_binary_path=str(binary),runtime_binary_hash=cfg['runtime_binary_hash'],cuda_sync_exit_code=0,ranks=2,mpi_sum=mpi.allreduce(1,mpi.SUM)),indent=2))
        else:mpi.allreduce(1,mpi.SUM)
        mpi.Barrier();mpi.check(mpi.lib.MPI_Finalize());return
    import numpy as np
    dp=cfg['liquid_parameters'];domain=cfg['liquid_domain'];u=mir.Mirheo((1,1,1),domain,debug_level=1,log_filename='fluid',no_splash=True)
    pv=mir.ParticleVectors.ParticleVector('fluid',mass=dp['mass']);u.registerParticleVector(pv,mir.InitialConditions.Uniform(dp['number_density']))
    if compute:
        n=len(pv.getCoordinates());rng=np.random.default_rng(cfg['seed']);v=rng.normal(0,(dp['kBT']/dp['mass'])**.5,(n,3));v-=v.mean(axis=0);pv.setVelocities(v.tolist())
        initial=np.asarray(pv.getCoordinates());assert n==512
    interaction=mir.Interactions.Pairwise('dpd',dp['rc'],kind='DPD',a=dp['a'],gamma=dp['gamma'],kBT=dp['kBT'],power=dp['power'])
    u.registerInteraction(interaction);u.setInteraction(interaction,pv,pv)
    vv=mir.Integrators.VelocityVerlet('vv');u.registerIntegrator(vv);u.setIntegrator(vv,pv)
    u.registerPlugins(mir.Plugins.createParticleChecker('finite_check',check_every=10))
    sync();mpi.Barrier();start=time.perf_counter();u.run(cfg['liquid_steps'],dt=cfg['dt']);sync();mpi.Barrier();elapsed=mpi.allreduce(time.perf_counter()-start)
    if compute:
        actual_steps=int(u.getState().current_step);actual_time=float(u.getState().current_time)
        pos=np.asarray(pv.getCoordinates());vel=np.asarray(pv.getVelocities());ids=np.asarray(pv.get_indices())
        assert actual_steps==cfg['liquid_steps'] and actual_steps>0
        assert len(pos)==n and len(np.unique(ids))==n and np.isfinite(pos).all() and np.isfinite(vel).all()
        assert not np.array_equal(pos,initial),'NO_PARTICLE_MOTION'
        with Path('particles_final.csv').open('w') as f:
            writer=csv.writer(f);writer.writerow(['id','x','y','z','vx','vy','vz'])
            for i in np.argsort(ids):writer.writerow([int(ids[i]),*pos[i],*vel[i]])
        result=dict(status='NATIVE_LIQUID_PASS',actual_steps=actual_steps,actual_time=actual_time,dt=cfg['dt'],particles=n,finite=True,unique_ids=True,
            kinetic_temperature=float(dp['mass']*np.sum((vel-vel.mean(axis=0))**2)/(3*(n-1))),gpu_synchronized_solver_wall_s=elapsed,
            particle_motion=True,parameters=dp,domain=domain,ranks=2,gpus=1,runtime_binary_hash=cfg['runtime_binary_hash'],scientific_calibration_pass=False)
        Path('completion.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
    del u;gc.collect();mpi.Barrier();mpi.check(mpi.lib.MPI_Finalize())

if __name__=='__main__':main()
