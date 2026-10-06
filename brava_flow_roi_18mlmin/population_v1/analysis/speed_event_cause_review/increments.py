"""Separate background translation/curl and account for every accepted increment."""
import json
import numpy as np
import torch
from diagnose import HERE, save, account, OriginalRHS
from particle_3d.population_v1.gpu_resistance import solve_isolated

NAMES=['bulk_translation','bulk_rotation_from_curl','planar_shear','normal_resistance',
       'planar_resistance','kinematic_constraint','solver_residual']


def main():
    rhs=OriginalRHS()
    for pid in [1189,1329,642,822]:
        with np.load(HERE/f'measured_{pid}.npz') as f:z={k:f[k] for k in f.files}
        t=z['time_s'];R=z['R'];U=z['U'];bulk=np.diff(z['b_bulk'],axis=0)
        loads=np.stack([np.c_[bulk[:,:3],np.zeros_like(bulk[:,:3])],
            np.c_[np.zeros_like(bulk[:,:3]),bulk[:,3:]],np.diff(z['b_planar'],axis=0),
            -np.einsum('nij,nj->ni',np.diff(z['R_normal'],axis=0),U[:-1]),
            -np.einsum('nij,nj->ni',np.diff(z['R_planar'],axis=0),U[:-1]),
            np.diff(z['q'],axis=0),np.diff(z['residual'],axis=0)],axis=1)
        tensor=lambda x:torch.as_tensor(x,dtype=torch.float64,device='cuda')
        matrix=tensor(R[1:]);diag=tensor(z['diag'][1:]);free=tensor(z['free'][1:]);n=len(matrix)
        normals=torch.zeros((n,32,3),dtype=torch.float64,device='cuda')
        counts=torch.zeros(n,dtype=torch.int64,device='cuda')
        response=[]
        for k in range(len(NAMES)):
            v,_,audit=solve_isolated(matrix,tensor(loads[:,k]),diag,free,normals,counts)
            assert not bool(audit['fallback'].any())
            response.append(v.cpu().numpy())
        response=np.stack(response,axis=1);error=response.sum(1)-np.diff(U,axis=0)
        np.savez_compressed(HERE/f'increments_{pid}.npz',time_s=t[1:],names=NAMES,
                            generalized_loads=loads,generalized_velocity_contributions=response)
        result=json.loads((HERE/f'result_{pid}.json').read_text())
        # Replace earlier coarser accounting with the split of the same actual equations.
        for index,entry in enumerate(result['accounting']):
            i,j=np.searchsorted(t,[entry['from_state']['time_s'],entry['to_state']['time_s']])
            result['accounting'][index]=dict(label=entry['label'],**account(z,int(i),int(j),t))
        a=result['accounting'][0];i,j=np.searchsorted(t,[a['from_state']['time_s'],a['to_state']['time_s']])
        sums=response[i:j].sum(0)
        result['pathwise_increment_accounting']=dict(reference='R at each later adjacent state',
            names=NAMES,summed_translation_mm_s=(sums[:,:3]*1e3).tolist(),
            summed_rotation_s_inv=sums[:,3:].tolist(),
            maximum_reconstruction_error_translation_m_s=float(np.max(abs(error[:,:3]))),
            maximum_reconstruction_error_rotation_s_inv=float(np.max(abs(error[:,3:]))),
            actual_translation_delta_mm_s=((U[j,:3]-U[i,:3])*1e3).tolist(),
            actual_rotation_delta_s_inv=(U[j,3:]-U[i,3:]).tolist())
        for interface in result['interfaces']['measurements']:
            for probe in interface['approaches']:
                if probe['status']!='MEASURED':continue
                if not (probe['from_state']['wall_weight'] or probe['to_state']['wall_weight']):continue
                positions=[probe['from_state']['x_m'],probe['to_state']['x_m']]
                packet=rhs.batch(positions,result['radius_m'])
                replacement=account(packet,0,1)
                assert abs(np.linalg.norm(packet['U'][1,:3]-packet['U'][0,:3])-probe['delta_RHS_m_s'])<1e-12
                probe.update(replacement)
        save(f'result_{pid}.json',result)
        print('INCREMENTS',pid,'max_translation_error',float(np.max(abs(error[:,:3]))),flush=True)


if __name__=='__main__':main()
