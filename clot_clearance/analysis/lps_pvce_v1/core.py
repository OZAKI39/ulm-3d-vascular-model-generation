"""PVCE task adapter: reuse locked energy; native virtual cut via public atom arrays.
No production edits, no material parameter edits, one rank and fixed atom ordering.
"""
from pathlib import Path
import os,math,json,hashlib
import numpy as np
from pd_clot.lammps_fragmentation.macro_lab_v1 import geometry,Energy,K,G,RHO,DELTA,save
from pd_clot.lammps_fragmentation.interface_model import write_typed_data
from pd_clot.lammps_fragmentation.lammps_writer import force_conversion
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'verification/lps_pvce_energy_release_feasibility_v1/20261011T132103Z'
DX=50e-6

def model(tip):
 g,_=geometry(DX,'M4');X=g['X'];b=g['bonds'];a,c=X[b[:,0]],X[b[:,1]]
 cross=(a[:,2]<0)!=(c[:,2]<0)
 frac=np.divide(-a[:,2],c[:,2]-a[:,2],out=np.zeros(len(b)),where=c[:,2]!=a[:,2])
 xx=a[:,0]+frac*(c[:,0]-a[:,0]);cut=cross&(xx<tip-1e-12)
 # Half-open crack x<tip; 1pm tie guard enforces exact lattice-plane convention.
 alive=np.ascontiguousarray(~cut,dtype=np.uint8)
 e=Energy(g,DX,ROOT/'analysis/material_macro_audit_v1/energy.so')
 mask=np.zeros_like(X,dtype=bool);low=np.isclose(X[:,2],X[:,2].min(),rtol=0,atol=1e-12);high=np.isclose(X[:,2],X[:,2].max(),rtol=0,atol=1e-12);mask[low|high,2]=True
 return g,e,alive,mask,low,high,xx,cross

class Native:
 def __init__(self,out,g,alive):
  from lammps import lammps
  out.mkdir(parents=True,exist_ok=True);self.g=g;self.alive=alive;self.n=len(g['X']);self.out=out
  write_typed_data(out/'body.data',g,dict(dx_m=DX,delta_m=DELTA,rho_kg_m3=RHO))
  text=f'''units si
 dimension 3
 boundary m m m
 atom_style bond
 atom_modify map array sort 0 0
 newton on off
 special_bonds lj 0 1 1 coul 1 1 1
 lattice sc {DX:.17g}
 read_data {out/'body.data'} extra/bond/per/atom 320 extra/special/per/atom 360
 neighbor {2*DX:.17g} bin
 neigh_modify every 1 delay 0 check no
 comm_modify cutoff {3*DELTA:.17g}
 fix props all property/atom d_vfrac d_theta d_damage ghost yes
 set group all d_vfrac {DX**3:.17g}
 pair_style bpm/peri contact_distance_factor 1.0
 pair_coeff * * {18*K/(math.pi*DELTA**4):.17g} {DELTA:.17g}
 bond_style bpm/peri
 bond_coeff 1 lps {K:.17g} {G} {DELTA:.17g} .15 0
 bond_coeff 2 lps {K:.17g} {G} {DELTA:.17g} .140624859375 0
 bond_coeff 3 lps {K:.17g} {G} {DELTA:.17g} 1 0
 velocity all set 0 0 0 units box
 {force_conversion()}
 timestep 5e-7
 compute contactPE all pe pair
 compute bondPE all pe bond
 thermo 1
 thermo_style custom step atoms bonds ke c_contactPE c_bondPE
 thermo_modify lost error lost/bond error
 run 0
 '''
  (out/'in.lammps').write_text(text+'\n# Python virtual-cut arrays, then run 1 pre no post no per force evaluation. No time integration fix.\n')
  self.L=lammps(cmdargs=['-log',str(out/'log.lammps'),'-screen','none']);self.L.commands_string(text)
  L=self.L;assert np.array_equal(L.numpy.extract_atom('id')[:self.n],np.arange(1,self.n+1))
  cutpairs=set(map(tuple,(g['bonds'][alive==0]+1).tolist()));count=0
  nb=L.extract_atom('num_bond');bt=L.extract_atom('bond_type');ba=L.extract_atom('bond_atom');ns=L.extract_atom('nspecial');sp=L.extract_atom('special')
  for i in range(self.n):
   remove=set()
   for j in range(nb[i]):
    pair=tuple(sorted((i+1,ba[i][j])))
    if pair in cutpairs:bt[i][j]=0;remove.add(ba[i][j]);count+=1
   # Remove only declared cut partners from the 1-2 special list. Later levels
   # are unused (weights1); keep their storage counts consistent without init.
   if remove:
    original=[sp[i][j] for j in range(ns[i][2])];deleted=[j for j,v in enumerate(original[:ns[i][0]]) if v in remove]
    kept=[v for j,v in enumerate(original) if j not in set(deleted)]
    for j,v in enumerate(kept):sp[i][j]=v
    for j in range(3):ns[i][j]-=len(deleted)
  assert count==2*len(cutpairs);self.cut_entries=count
 def evaluate(self,x):
  L=self.L;n=self.n;L.numpy.extract_atom('x')[:n]=x
  # pre no avoids BondBPMPeri::init_style resetting the common reference m.
  # every1/checkno rebuilds active neighbor/bond lists using zeroed bond types.
  L.command('run 1 pre no post no')
  f=L.numpy.extract_atom('f')[:n].copy();theta=L.numpy.extract_atom('d_theta')[:n].copy()
  bt=L.extract_atom('bond_type');nb=L.extract_atom('num_bond');ba=L.extract_atom('bond_atom');pairs=[]
  for i in range(n):
   for j in range(nb[i]):
    if bt[i][j]>0 and ba[i][j]>i+1:pairs.append((i,ba[i][j]-1))
  expected=self.g['bonds'][self.alive.astype(bool)];actual=np.array(sorted(pairs),int)
  assert np.array_equal(expected,actual),'Unprescribed bond loss or topology mismatch'
  return f,dict(contact_J=float(L.extract_compute('contactPE',0,0))*DX**3,native_bond_PE_times_V_J=float(L.extract_compute('bondPE',0,0))*DX**3,active_bonds=len(actual),natural_broken=0,kinetic_J=float(L.get_thermo('ke'))),theta
 def close(self):self.L.close()
