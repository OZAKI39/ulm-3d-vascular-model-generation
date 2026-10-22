"""Start only after real-data geometry AND full trajectory/audit parity pass."""
import hashlib,json
import subprocess,sys
from pathlib import Path
import campaign
from native_adapter import accelerated_job,native_closest


def main():
    proof=json.loads((campaign.HERE/'data/native_kernel_parity.json').read_text())
    assert proof['PASS'] and proof['kernel_bitwise_pass']
    exact=json.loads((campaign.HERE/'data/native_exact_bytes_verification.json').read_text())
    assert exact['PASS_FOR_ACTUAL_HANDOFF_DEPENDENCY']
    mesh=json.loads((campaign.HERE/'data/gpu_mesh_validation.json').read_text())
    assert mesh['CPU_CUDA_gradient_parity'] and mesh['all_finite']
    subprocess.run([sys.executable,'-B',str(campaign.HERE/'scripts/validate_query_cache.py')],check=True)
    cache=json.loads((campaign.HERE/'data/query_cache_parity.json').read_text())
    assert cache['PASS']
    from query_cache import install
    old_setup=campaign.setup
    def setup():
        root,env,protected=old_setup()
        install(env)
        return root,env,protected
    campaign.setup=setup
    files=['campaign.py','native_adapter.py','closest_native.cpp','production_entry.py','query_cache.py','validate_query_cache.py']
    identities={name:hashlib.sha256((campaign.HERE/'scripts'/name).read_bytes()).hexdigest() for name in files}
    original=campaign.bind_job
    def bind(env,identity):
        identity.update(acceleration_sources=identities,
            acceleration='COMPILED_TRIANGLE_CLOSEST_AND_EXACT_READ_ONLY_QUERY_CACHE',
            query_cache_proof_sha256=hashlib.sha256((campaign.HERE/'data/query_cache_parity.json').read_bytes()).hexdigest(),
            acceleration_proof_sha256=hashlib.sha256((campaign.HERE/'data/native_kernel_parity.json').read_bytes()).hexdigest())
        # Compile/load in parent before forking; no CUDA context enters workers.
        native_closest(env.wall.triangles[0,0],env.wall.triangles[:1])
        return accelerated_job(original(env,identity))
    campaign.bind_job=bind  # campaign-only orchestration, no scientific module patch
    campaign.run('production')


if __name__=='__main__':main()
