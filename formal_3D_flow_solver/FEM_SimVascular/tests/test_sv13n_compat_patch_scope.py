from sv13n_support import *
def test_adopted_patch_scope():
 d=accepted('compatibility_adapter');compatibility_gate((ROOT/'patches/sv1_3n/svmp_petsc325_compat.patch').read_text(),d)
@pytest.mark.parametrize('body',['+++ b/Code/Source/solver/fluid.cpp\n+ double rho=1;','+++ b/Code/Source/solver/petsc_impl.cpp\n+ (void)PetscFinalize();','+++ b/Code/Source/solver/petsc_impl.cpp\n+ KSPSetTolerances(ksp,1e-2,0,0,5);'])
def test_science_or_swallowed_error_rejected(body):
 manifest=dict(repair_iterations=1,patch_sha256=hashlib.sha256(body.encode()).hexdigest(),modified_files=['Code/Source/solver/petsc_impl.cpp'])
 with pytest.raises(GateError):compatibility_gate(body,manifest)
