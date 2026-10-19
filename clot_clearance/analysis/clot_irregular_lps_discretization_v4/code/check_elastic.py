"""Small affine patch checks, with frozen tolerances; not a modulus calibration."""
from common import *
from scipy.spatial import cKDTree
from native import native
import time

def patch_reference(full, selected, regular):
 X=(full['regular_X'] if regular else full['X'])[selected].copy()
 vol=np.full(len(X),V) if regular else full['V'][selected].copy()
 b=cKDTree(X).query_pairs(HORIZON,output_type='ndarray');b.sort(1);b=b[np.lexsort((b[:,1],b[:,0]))]
 i,j=b.T;r=np.linalg.norm(X[j]-X[i],axis=1);h=np.minimum(1.,(HORIZON+DX/2-r)/DX)
 m=np.bincount(i,weights=h*r*vol[j],minlength=len(X))+np.bincount(j,weights=h*r*vol[i],minlength=len(X))
 return dict(X=X,bonds=b,r0=r,taper=h,m=m,V=vol,bond_types=np.ones(len(b),int),initial_cut=np.zeros(len(b),bool),atom_ids=np.arange(1,len(X)+1),original_bond_ids=np.arange(1,len(b)+1),clot_count=len(X),source_original_atom_ids=selected+1)

def write_patch(path,ref):
 p=path/'provenance';p.mkdir(parents=True,exist_ok=False);np.savez_compressed(p/'reference.npz',**ref)
 N=len(ref['X']);b=ref['bonds'];NB=len(b)
 with (p/'reference.native').open('w') as f:
  f.write(f'LPSVAR_REFERENCE_V2\n{N} {NB} {V:.17g}\n')
  for i,(m,v) in enumerate(zip(ref['m'],ref['V'])):f.write(f'{i+1} {m:.17g} {v:.17g}\n')
  for k,((i,j),r,h) in enumerate(zip(b,ref['r0'],ref['taper'])):f.write(f'{k+1} {i+1} {j+1} 1 {r:.17g} {h:.17g} 0 0 0\n')
 with (p/'body.data').open('w') as f:
  f.write(f'V4 affine patch, original member IDs in reference.npz\n\n{N} atoms\n{NB} bonds\n\n2 atom types\n3 bond types\n\n')
  for axis in 'xyz':f.write(f'-0.003 0.003 {axis}lo {axis}hi\n')
  f.write(f'\nMasses\n\n1 {RHO*V:.17g}\n2 {RHO*V:.17g}\n\nAtoms # bond\n\n')
  for k,x in enumerate(ref['X']):f.write(f'{k+1} 1 1 '+' '.join(f'{v:.17g}' for v in x)+'\n')
  f.write('\nBonds\n\n')
  for k,(i,j) in enumerate(b):f.write(f'{k+1} 1 {i+1} {j+1}\n')

def main():
 start=time.perf_counter();full=load_npz(O/'provenance/reference.npz');g=json.loads((O/'ADMISSION_GATES.json').read_text())
 X=full['regular_X'];selected=np.flatnonzero((X[:,0]>=50e-6-1e-14)&(X[:,0]<=450e-6+1e-14)&(abs(X[:,1])<=200e-6+1e-14)&(abs(X[:,2])<=200e-6+1e-14))
 assert len(selected)==729
 result=dict(test='Small 9x9x9 affine patches; all coordinates prescribed, no free relaxation or fitted material constants',original_material_IDs=selected+1,strain=1e-4,patches={},gates=g)
 for name in ['regular','irregular']:
  ref=patch_reference(full,selected,name=='regular');case=O/'tests/elastic'/name;write_patch(case,ref);a,_=native(case,install_properties=False,surface=False);model=energy_model(ref)
  zero=a.static();assert abs(zero['force']).max()<g['reference_force_max_N'];rows=[]
  for mode in ['hydro','xy','xz','yz']:
   F=np.eye(3);eps=1e-4
   if mode=='hydro':F*=1+eps
   else:F['xyz'.index(mode[0]),'xyz'.index(mode[1])]=eps
   x=ref['X']@F.T;s=a.static(x);z=model.evaluate(x);ferr=float(abs(s['material_force']-z['force']).max());uerr=abs(s['native_U']-z['U'])
   assert ferr<=g['native_independent_force_max_N'] and uerr<=g['native_independent_energy_max_J']
   assert s['Uc']==0 and abs(s['force']-s['material_force']).max()<1e-14
   P=-z['force'].T@ref['X']/ref['V'].sum()
   modulus_reaction=float(np.trace(P)/(9*eps) if mode=='hydro' else P['xyz'.index(mode[0]),'xyz'.index(mode[1])]/eps)
   modulus_energy=2*z['U']/(ref['V'].sum()*eps**2*(9 if mode=='hydro' else 1))
   direction=z['force']/np.linalg.norm(z['force']);dh=1e-9
   fd=(model.evaluate(x+dh*direction)['U']-model.evaluate(x-dh*direction)['U'])/(2*dh);expected=-np.sum(z['force']*direction);fdrel=abs(fd-expected)/abs(expected)
   assert fdrel<=g['energy_gradient_relative_error_max']
   rows.append(dict(mode=mode,modulus_reaction_Pa=modulus_reaction,modulus_energy_Pa=float(modulus_energy),U_J=z['U'],first_piola_Pa=P,force_max_N=float(abs(z['force']).max()),resultant_N=z['force'].sum(0),native_force_error_N=ferr,native_energy_error_J=uerr,energy_gradient_relative_error=float(fdrel)))
   np.savez_compressed(case/(mode+'.npz'),x=x,force=s['force'],material_force=s['material_force'],native_U=s['native_U'],independent_U=z['U'])
  a.close();result['patches'][name]=dict(nodes=len(ref['X']),bonds=len(ref['bonds']),volume_m3=float(ref['V'].sum()),reference_force_max_N=float(abs(zero['force']).max()),rows=rows)
 changes={r['mode']:ir['modulus_energy_Pa']/r['modulus_energy_Pa']-1 for r,ir in zip(result['patches']['regular']['rows'],result['patches']['irregular']['rows'])}
 spread={name:float(np.ptp([r['modulus_energy_Pa'] for r in p['rows'][1:]])/np.mean([r['modulus_energy_Pa'] for r in p['rows'][1:]])) for name,p in result['patches'].items()}
 result.update(relative_modulus_changes=changes,shear_directional_spread=spread,wall_s=time.perf_counter()-start)
 result['pass']=max(abs(v) for v in changes.values())<=g['affine_bulk_and_shear_relative_change_max'] and max(spread.values())<=g['shear_directional_spread_max']
 dump(O/'ELASTIC_PATCH_CHECK.json',result);print(json.dumps({k:result[k] for k in ['relative_modulus_changes','shear_directional_spread','pass']},default=cv),flush=True)
 if not result['pass']:raise RuntimeError('Frozen affine elastic admission gate failed; no formal dynamics authorized')

if __name__=='__main__':main()
