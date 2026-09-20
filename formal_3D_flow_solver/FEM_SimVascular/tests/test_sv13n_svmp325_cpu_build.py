from sv13n_support import *
def test_cpu_linkage():
 d=accepted('svmp_cpu_build');l=accepted('svmp_cpu_link');linkage_gate(l,d['PETSc_prefix'])
def test_wrong_version_rejected():
 d=read('svmp_cpu_link');d['version']='3.19.6'
 with pytest.raises(GateError):linkage_gate(d,read('svmp_cpu_build')['PETSc_prefix'])
