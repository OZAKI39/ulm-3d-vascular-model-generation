from common import *
def write_patch(path,ref):
 p=path/'provenance';p.mkdir(parents=True,exist_ok=False);np.savez_compressed(p/'reference.npz',**ref)
 N=len(ref['X']);b=ref['bonds'];NB=len(b)
 with (p/'reference.native').open('w') as f:
  f.write(f'LPSVAR_REFERENCE_V2\n{N} {NB} {V:.17g}\n')
  for i,(m,v) in enumerate(zip(ref['m'],ref['V'])):f.write(f'{i+1} {m:.17g} {v:.17g}\n')
  for k,((i,j),r,h) in enumerate(zip(b,ref['r0'],ref['taper'])):f.write(f'{k+1} {i+1} {j+1} 1 {r:.17g} {h:.17g} 0 0 0\n')
 with (p/'body.data').open('w') as f:
  f.write(f'V4 affine patch, original member IDs in reference.npz\n\n{N} atoms\n{NB} bonds\n\n2 atom types\n3 bond types\n\n')
  for axis in 'xyz':f.write(f'-0.03 0.03 {axis}lo {axis}hi\n')
  f.write(f'\nMasses\n\n1 {RHO*V:.17g}\n2 {RHO*V:.17g}\n\nAtoms # bond\n\n')
  for k,x in enumerate(ref['X']):f.write(f'{k+1} 1 1 '+' '.join(f'{v:.17g}' for v in x)+'\n')
  f.write('\nBonds\n\n')
  for k,(i,j) in enumerate(b):f.write(f'{k+1} 1 {i+1} {j+1}\n')
