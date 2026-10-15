"""Extra readback of original native atom/bond IDs, with no LAMMPS or ideal energy call."""
from common import *
import shutil,time
start=time.perf_counter(); data=np.load(W/'reference_topology.npz'); X=data['X']; bonds=data['bonds']; atoms=[]; pb=[]; section=None
source=W/'preflight/full/body.data'
with source.open() as fp:
    for line in fp:
        line=line.strip()
        if line.startswith('Atoms'):section='atoms';continue
        if line=='Bonds':section='bonds';continue
        if not line:continue
        if section=='atoms':
            x=line.split(); atoms.append([int(x[0]),int(x[1]),int(x[2]),*map(float,x[3:6])])
        elif section=='bonds':pb.append(list(map(int,line.split())))
aa=np.array(atoms); bb=np.array(pb,dtype=np.int64)
assert np.array_equal(aa[:,0],np.arange(len(X))+1) and np.array_equal(aa[:,3:6],X)
assert np.array_equal(bb[:,0],data['bond_reference_id']) and np.array_equal(bb[:,2:4]-1,bonds) and np.all(bb[:,1]==1)
inp=W/'preflight/full/in.lammps'; text=inp.read_text()
vol=float(next(l.split()[-1] for l in text.splitlines() if l.strip().startswith('set group all d_vfrac ')));assert vol==V
coeff=next(l.split()[4:] for l in text.splitlines() if l.strip().startswith('bond_coeff 1 lps '))
# The token list for bond_coeff is checked explicitly below; no assumptions from a rounded report.
tokens=next(l.split() for l in text.splitlines() if l.strip().startswith('bond_coeff 1 lps '))
assert list(map(float,tokens[3:]))==[K,G,DELTA,S00,0.]
shutil.copy2(inp,O/'provenance/wide_full_in.lammps')
dump(O/'NATIVE_INPUT_IDENTITY_AUDIT.json',{'status':'PASS','actual_body_data_path':str(source),'body_data_sha256':sha(source),'in_lammps_sha256':sha(inp),'actual_atom_IDs_read':len(aa),'actual_unique_bond_IDs_read':len(bb),'atom_ID_and_coordinate_bitwise_match':True,'bond_ID_endpoints_and_bulk_type_exact_match':True,'native_vfrac_m3':vol,'native_bulk_coefficients':list(map(float,tokens[3:])),'simulation_box_note':'Native nonperiodic padded box is larger than the material specimen. Reported +/-2.5 mm X and +/-0.5 mm Z are material-cell outer bounds, not padded LAMMPS box walls.','runtime_s':time.perf_counter()-start,'additional_formal_integrals':0,'native_calls':0,'scope':'Supplementary archival input readback; preregistration already established identity by original NPZ SHA, original generator/writer, saved native reference and full independent stencil.'})
print('Actual native data IDs, coordinates, complete bond IDs, volumes and material input PASS')
