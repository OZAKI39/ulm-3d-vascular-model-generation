"""Prescribed small-affine response of the actual Y material mask, no fitting."""
from common import *
from native import native
from patch_io import write_patch
import time

def material_reference(full):
 n=int(full['clot_count']);X=full['X'][:n].copy();vol=full['V'][:n].copy();b=full['bonds'][full['bond_types']==1].copy()
 i,j=b.T;r=np.linalg.norm(X[j]-X[i],axis=1);h=np.minimum(1.,(HORIZON+DX/2-r)/DX)
 m=np.bincount(i,weights=h*r*vol[j],minlength=n)+np.bincount(j,weights=h*r*vol[i],minlength=n)
 return dict(X=X,bonds=b,r0=r,taper=h,m=m,V=vol,bond_types=np.ones(len(b),int),initial_cut=np.zeros(len(b),bool),atom_ids=np.arange(1,n+1),original_bond_ids=np.arange(1,len(b)+1),clot_count=n)
def main():
 start=time.perf_counter();results={}
 for label,source in [('regular',Y_OLD),('irregular',O)]:
  ref=material_reference(load_npz(source/'provenance/reference.npz'));case=O/'checks/affine'/label;write_patch(case,ref)
  ns,_=native(case,install_properties=False,surface=False);model=energy_model(ref);zero=ns.static();assert abs(zero['material_force']).max()<1e-14
  center=ref['V']@ref['X']/ref['V'].sum();rows=[]
  for mode in ['hydro','xy','xz','yz']:
   Fm=np.eye(3);eps=1e-4
   if mode=='hydro':Fm*=1+eps
   else:Fm['xyz'.index(mode[0]),'xyz'.index(mode[1])]=eps
   x=(ref['X']-center)@Fm.T+center;s=ns.static(x);z=model.evaluate(x)
   ferr=float(abs(z['force']-s['material_force']).max());uerr=abs(z['U']-s['native_U'])
   assert ferr<1e-12 and uerr<1e-17 and s['Uc']==0
   direction=z['force']/np.linalg.norm(z['force']);dh=1e-9;fd=(model.evaluate(x+dh*direction)['U']-model.evaluate(x-dh*direction)['U'])/(2*dh)
   exact=-np.sum(z['force']*direction);fdrel=float(abs(fd-exact)/abs(exact));assert fdrel<2e-5
   modulus=float(2*z['U']/(ref['V'].sum()*eps**2*(9 if mode=='hydro' else 1)))
   P=-z['force'].T@(ref['X']-center)/ref['V'].sum()
   rows.append(dict(mode=mode,modulus_energy_Pa=modulus,first_piola_Pa=P,U_J=z['U'],native_force_error_N=ferr,native_energy_error_J=uerr,gradient_relative_error=fdrel))
  ns.close();results[label]=dict(nodes=len(ref['X']),bonds=len(ref['bonds']),volume_m3=float(ref['V'].sum()),rows=rows)
 changes={r['mode']:i['modulus_energy_Pa']/r['modulus_energy_Pa']-1 for r,i in zip(results['regular']['rows'],results['irregular']['rows'])}
 passed=max(abs(v) for v in changes.values())<=.15
 dump(O/'ELASTIC_CHECK.json',dict(status='PASS' if passed else 'FAIL',test='Same original 1927-node Y material mask, all coordinates affine-prescribed, anchors omitted from diagnostic coupon; original/new void-filtered bulk graphs and actual Vi; no fitted parameters',patches=results,relative_modulus_changes=changes,relative_change_gate=.15,scope='Relative discretization diagnostic. The asymmetric finite sample is not a macroscopic modulus or isotropy calibration.',wall_s=time.perf_counter()-start))
 assert passed,'Relative elastic discretization gate failed'
 print('Y_V4_ELASTIC_PASS',changes,flush=True)
if __name__=='__main__':main()
