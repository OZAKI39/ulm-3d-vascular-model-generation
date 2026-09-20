from sv13p_support import *
def test_isolated_hypre_gpu_build():
 d=read('remote/petsc_hypre_build')
 assert '/sv1_3p/external/petsc325_cuda13_hypre/' in d['prefix']
 assert d['version']=='3.25.5' and d['package_version']=='3.1.0'
 assert d['package_commit']=='85b779557005b2eb94c231c1b516e988b87f4e53' and len(d['package_archive_sha256'])==64
 if d['status']!='PASS':assert candidate('P2')['result'].startswith('BUILD_');return
 assert d['source_unmodified'] and not d['source_modifications']
 f=d['hypre_features'];assert f['HYPRE_USING_CUDA'] and f['HYPRE_USING_GPU'] and f['HYPRE_USING_CUSPARSE']
 assert not f['HYPRE_USING_UNIFIED_MEMORY'] and not f['HYPRE_BIGINT'] and not f['HYPRE_MIXEDINT']
 accepted('remote/hypre_ilu_standalone');accepted('HYPRE_BASELINE_SMOKE_acceptance')
