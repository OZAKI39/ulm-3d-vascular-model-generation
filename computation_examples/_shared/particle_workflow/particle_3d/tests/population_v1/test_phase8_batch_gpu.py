from pathlib import Path
import os
import numpy as np
import pytest
torch=pytest.importorskip('torch')
pytestmark=pytest.mark.skipif(not torch.cuda.is_available(),reason='Requires CUDA')
from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
from particle_3d.population_v1.gpu_field import GPUField
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.population_v1.gpu_batch import perfuse_batch_until
from particle_3d.population_v1.source import PerfusionSource,perfuse_until
from particle_3d.population_v1.provenance import canonical,write_json
from test_phase4 import source_env


def test_independent_realizations_preserve_clocks_and_sources_across_batch_sizes():
    env=source_env();gf=GPUField(env.field);gw=GPUWall(env.wall)
    def build(rids):
        engines=[GPUPopulationEngine(env,pairs=True,realization_id=i,gpu_field=gf,gpu_wall=gw) for i in rids]
        sources=[PerfusionSource(env,concentration_m3=2e14,seed=2026092908,realization_id=i) for i in rids]
        return engines,sources
    independent,source=build(range(4))
    for e,s in zip(independent,source):perfuse_until(e,s,.004)
    results=[]
    for batch_size in [2,4]:
        engines,sources=build(range(4));dispatches=[]
        for k in range(0,4,batch_size):dispatches.extend(perfuse_batch_until(engines[k:k+batch_size],sources[k:k+batch_size],[.004]*batch_size)['dispatches'])
        errors=[]
        for ref,got in zip(independent,engines):
            assert canonical(ref.source_ledger)==canonical(got.source_ledger)
            assert ref.pool.n==got.pool.n
            assert [r['physical_time_s'] for r in ref.snapshots]==[r['physical_time_s'] for r in got.snapshots]
            err=float((ref.pool.x[:ref.pool.n]-got.pool.x[:got.pool.n]).abs().max()) if got.pool.n else 0.
            assert err<1e-15
            errors.append(err)
        assert max(r['realizations'] for r in dispatches)>=2
        results.append(dict(batch_size=batch_size,maximum_position_errors_m=errors,
            largest_batched_population=max(r['particles'] for r in dispatches),source_ledgers_identical=True,
            physical_time_partitions_identical=True,accounting=[e.accounting() for e in engines]))
    write_json(Path(os.environ.get('POPULATION_VALIDATION_DIR','/tmp/population_v1_validation'))/'phase8_batch_gpu.json',
        dict(PASS=True,phase=8,realizations=4,results=results,concentration_role='EXPLICIT_SYNTHETIC_TEST_ONLY'))
