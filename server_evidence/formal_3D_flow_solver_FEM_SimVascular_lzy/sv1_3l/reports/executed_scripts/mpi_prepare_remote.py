"""Read the exact 4.1.6 configure interface before selecting any binding flags."""
import tarfile
from runner_remote import *
assert load('compiler_environment')['status']=='PASS'
manifest=load('MPI_source_manifest')
archive=BASE.parent/'sv1_3g/external/openmpi-4.1.6.tar.gz'
assert archive.is_file() and digest(archive)==manifest['sha256']
src=BASE/'external/openmpi-4.1.6';assert not src.exists()
src.parent.mkdir(parents=True,exist_ok=True)
with tarfile.open(archive) as t:t.extractall(src.parent,filter='data')
help_record=run(['./configure','--help'],'openmpi_configure_help',cwd=src)
assert okay(help_record)
write('openmpi_source',{'status':'PASS','source':manifest,'archive':str(archive),'source_directory':str(src),'fresh_extraction':True,'configure_help':help_record})
for line in text(help_record).splitlines():
    if any(term in line for term in ('mpi-fortran','Fortran','mpif.h','mpi_f08','use mpi','FC ','gl','opencl','nvml','cuda','hwloc')):print(line)
