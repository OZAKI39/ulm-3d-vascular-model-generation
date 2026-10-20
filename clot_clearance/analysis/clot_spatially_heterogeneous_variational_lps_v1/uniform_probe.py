"""Run independently against old/new libraries; compare in checks.py."""
from common import *
from native import native
import sys
name=sys.argv[1];ns,C=native(O/'checks'/name,install_properties=name!='old_uniform');ref=ns.ref;out={}
for it,strain in enumerate([.01,.149,.16]):
 x=ref['X']*(1+strain);ns.static(x);assert ns.lib.lpsvar_action(ns.l.lmp,2,it+1,None,None)==0;s=ns.static()
 for key in ['material_force','scalars','d','kappa']:out[f'{it}_{key}']=s[key]
out['m']=ns.field(4);out['original_bond_ids']=ns.field(9);out['on']=ns.field(13);out['off']=ns.field(14)
np.savez_compressed(O/'checks'/(name+'.npz'),**out);ns.close();print('UNIFORM_PROBE',name,flush=True)
