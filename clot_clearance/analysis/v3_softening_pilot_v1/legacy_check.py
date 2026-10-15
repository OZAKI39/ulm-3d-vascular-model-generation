"""Two static force evaluations of the same real V3 geometry; no dynamics."""
from common import *
import ctypes as ct
from lammps import lammps
ref=dict(np.load(O/'provenance/reference.npz'));N=len(ref['X']);n=int(ref['clot_count']);X=ref['X']
deck=(O/'provenance/in.lammps').read_text().split('fix tangent_props',1)[0]
deck=deck.replace('read_data body.data',f'read_data {O/"provenance/body.data"}').replace('atom_modify map array','atom_modify map array sort 0 0')
deck=deck.replace('fix si_force all setforce v_si_fx v_si_fy v_si_fz','')
L=lammps(cmdargs=['-log',str(O/'preflight/legacy.log'),'-screen','none']);L.commands_string(deck)
maps=sorted(set(s.split()[-1] for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps' in s));assert all(sha(p)==sha(OLDLIB) for p in maps)
probe=ct.CDLL(str(OLDV3/'observer.so'));probe.observe.argtypes=[ct.c_void_p,ct.c_int];probe.observe.restype=ct.POINTER(ct.c_double)
states={};Xp=X.copy();ids0=np.arange(n);Xp[:n,0]+=1e-8*np.sin(ids0*.13);Xp[:n,1]+=1e-8*np.cos(ids0*.11);Xp[:n,2]+=1e-8*np.sin(ids0*.17)
for label,pos in [('initial',X),('perturbed',Xp)]:
 ids=L.numpy.extract_atom('id')[:N].copy();L.numpy.extract_atom('x')[:N]=pos[ids-1];L.command('run 0');ids=L.numpy.extract_atom('id')[:N].copy();order=np.argsort(ids)
 states[label+'_force']=L.numpy.extract_atom('f')[:N].copy()[order]*V;states[label+'_m']=np.ctypeslib.as_array(probe.observe(L.lmp,0),shape=(N,)).copy()[order]
 assert int(L.get_thermo('bonds'))==len(ref['bonds'])
np.savez_compressed(O/'preflight/legacy_force_parity.npz',X=X,perturbed=Xp,**states)
dump(O/'preflight/LEGACY_STATIC_IDENTITY.json',dict(status='PASS',physical_steps=0,static_states=2,library_maps=maps,library_sha256=sha(OLDLIB),input_body_sha256=sha(O/'provenance/body.data'),no_dynamics=True))
L.close();print('LEGACY_STATIC_COMPLETE')
