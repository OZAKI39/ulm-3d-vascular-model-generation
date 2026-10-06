"""Exact deployed TET4 metric, evaluated on existing accepted/Gauss samples."""
from common import *
import csv,gc
import pyvista as pv
RHO=1056.;MU=.00345312;DT=.01;CTM=1.;CTC=36.

def metric(x):
    # nn::gnn, after nn solver->basis order {1,2,3,0}: origin=node3.
    J=np.swapaxes(x[...,:3,:]-x[...,3:4,:],-1,-2)
    inv=np.linalg.inv(J)
    return np.swapaxes(inv,-1,-2)@inv

def terms(K,u):
    kT=np.full(len(u),4*(CTM/DT)**2)
    kU=np.einsum('ni,nij,nj->n',u,K,u)
    kS=CTC*(MU/RHO)**2*np.sum(K*K,axis=(1,2)); total=kT+kU+kS
    tauM=1/(RHO*np.sqrt(total));tauC=1/(tauM*np.trace(K,axis1=1,axis2=2))
    steadyM=1/(RHO*np.sqrt(kU+kS));steadyC=1/(steadyM*np.trace(K,axis1=1,axis2=2))
    return np.c_[kT,kU,kS,tauM,tauC,kT/total,kU/total,kS/total,steadyM,steadyC]
COLS=['kT_s_inv2','kU_s_inv2','kS_s_inv2','tauM_m3_s_kg','tauC_Pa_s','fT','fU','fS','hypothetical_omega0_tauM','hypothetical_omega0_tauC']

def main():
    protocol=read(PREVIOUS/'EVENT_FIELD_PROTOCOL.json');inputs=read(HERE/'INPUT_SHA256.json')
    blocks=[];summaries={};orientation={}
    def reg(p):inputs[str(p)]=sha(p)
    def add(name,pid,kind,x,ids,u,phase,arc,weight):
        valid=ids>=0;v=np.full((len(x),10),np.nan);K=np.full((len(x),3,3),np.nan)
        K[valid]=metric(points[tet[ids[valid]]]);v[valid]=terms(K[valid],u[valid])
        blocks.append(dict(field=np.full(len(x),name),pid=np.full(len(x),pid),kind=np.full(len(x),kind),
             x_m=x,cell=ids,u_m_s=u,phase=phase,arc_m=arc,arc_weight_m=weight,Kxi_m_inv2=K,values=v,valid=valid))
        return v
    for name in NAMES:
        mp,fp=paths(name);reg(mp);reg(fp);grid=pv.read(mp);flow=pv.read(fp)
        points=np.asarray(grid.points,float);tet=grid.cells.reshape(-1,5)[:,1:].copy();u_nodes=np.asarray(flow.point_data['Velocity'],float)
        s=np.einsum('ij,ij->i',points[tet[:,3]]-points[tet[:,2]],np.cross(points[tet[:,1]]-points[tet[:,0]],points[tet[:,2]]-points[tet[:,1]]))
        assert not np.any(s==0);positive=s>0;tet[positive,:2]=tet[positive,:2][:,::-1]
        orientation[name]=dict(input_positive=int(positive.sum()),input_negative=int((~positive).sum()),reader_swap='0<->1 when positive; origin3 unchanged',metric_origin=3)
        summaries[name]={}
        for pid in PIDS:
            e=protocol['events'][str(pid)];path=np.load(PREVIOUS/f'path_{pid}.npz');reg(PREVIOUS/f'path_{pid}.npz')
            raw=np.load(CAUSE/f'original_{pid}.npz');reg(CAUSE/f'original_{pid}.npz')
            for kind in ['accepted','path']+(['rise'] if pid==822 else []):
                p=PREVIOUS/f'{name}_{kind}_{pid}.npz';reg(p);z=np.load(p);ids=z['cell'];u=z['flow'][:,:3]
                if kind=='accepted':
                    x=raw['x_m'];arc=path['original_arc_m'];phase=np.clip(np.searchsorted(e['phase_boundaries_s'],raw['physical_time_s'],side='right')-1,0,2);w=np.zeros(len(x))
                else:
                    c=np.load(PREVIOUS/f'{kind}_{pid}.npz');reg(PREVIOUS/f'{kind}_{pid}.npz');x=c['x_m'];arc=c['arc_m'];phase=c['phase'] if kind=='path' else np.full(len(x),-1);w=c['arc_weight_m'] if kind=='path' else np.zeros(len(x))
                v=add(name,pid,kind,x,ids,u,phase,arc,w)
                if kind=='accepted':
                    reps={}
                    for tag in ['A','B']:
                        ix=np.flatnonzero(np.all(x==e[tag]['x_m'],axis=1));assert len(ix)==1;ii=int(ix[0]);reps[tag]=dict(index=ii,cell=int(ids[ii]),x_m=x[ii].tolist(),u_m_s=u[ii].tolist(),**dict(zip(COLS,v[ii].tolist())))
                    summaries[name][str(pid)]=dict(representatives=reps,accepted_count=len(x),invalid=int((ids<0).sum()))
                if kind=='path':
                    # Keep all original order1/2/4 samples; statistics use order2 only, no repeated-weight bias.
                    summaries[name][str(pid)]['phase_distributions']={}
                    for ph in range(3):
                        sel=(phase==ph)&(path['quadrature_order']==2)&(ids>=0);vv=v[sel]
                        if len(vv):summaries[name][str(pid)]['phase_distributions'][e['phase_names'][ph]]=dict(N=len(vv),fT_min_median_max=np.quantile(vv[:,5],[0,.5,1]).tolist(),max_tauM_relative_change=float(np.max(vv[:,8]/vv[:,3]-1)),fT_arc_weighted=float(np.average(vv[:,5],weights=w[sel])))
            print('TERMS',name,pid,summaries[name][str(pid)]['representatives']['A']['fT'],flush=True)
        if name=='M0':
            # Four symmetric TET4 quadrature nodes plus centroid; no new mesh.
            a=(5+3*np.sqrt(5))/20;b=(5-np.sqrt(5))/20
            bary=np.full((4,4),b);np.fill_diagonal(bary,a);bary=np.vstack([np.full(4,.25),bary])
            for cid in [388769,138584,402304,400899,59816,212669]:
                x=bary@points[tet[cid]];u=bary@u_nodes[tet[cid]]
                add(name,1189 if cid==388769 else 822 if cid in [402304,400899,59816,212669] else -1,f'cell_{cid}_centroid_gauss4',x,np.full(5,cid),u,np.full(5,-1),np.full(5,np.nan),np.zeros(5))
        del grid,flow;gc.collect()
    out={k:np.concatenate([b[k] for b in blocks]) for k in blocks[0]}
    np.savez_compressed(HERE/'STABILIZATION_TERMS.npz',columns=COLS,**out)
    with (HERE/'STABILIZATION_TERMS.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['field','event','sample_kind','row','cell','phase','x_m','y_m','z_m','arc_m','valid']+COLS)
        for i in range(len(out['cell'])):w.writerow([out['field'][i],out['pid'][i],out['kind'][i],i,out['cell'][i],out['phase'][i],*out['x_m'][i],out['arc_m'][i],out['valid'][i],*out['values'][i]])
    put('STABILIZATION_TERMS.json',dict(status='OFFLINE_EXISTING_FIELD_EVALUATION',rows=len(out['cell']),invalid=int((~out['valid']).sum()),dt_s=DT,rho_kg_m3=RHO,mu_Pa_s=MU,ctM=CTM,ctC=CTC,
       Kxi='inverse(J).T @ inverse(J); J=[x0-x3,x1-x3,x2-x3] AFTER solver read_msh orientation',orientation=orientation,
       field_semantics='Frozen accepted steady velocity, not saved generalized-alpha intermediate quadrature payload. Formula exact; stage velocity reconstructed from steady output.',
       hypothetical_status='OFFLINE_HYPOTHETICAL_TIME_CONSISTENT_EVALUATION',hypothetical='omega=0 steady limit, same ctC=36 and Kxi; NOT actual CFD or full Jia code replication',
       actual_global_paper_omega=None,actual_global_paper_omega_reason='Global accepted partial-time acceleration payload not used; convective Gu is NOT paper numerator',
       distributions='All raw rows retained; phase distributions use original common path Gauss order2; min/median/max unweighted, arc mean weighted',fields=summaries,
       overall_fT_max=float(np.nanmax(out['values'][:,5])),overall_tauM_change_max=float(np.nanmax(out['values'][:,8]/out['values'][:,3]-1))))
    put('INPUT_SHA256.json',inputs)
    print('COMPLETE',len(out['cell']),np.nanmax(out['values'][:,5]),flush=True)
if __name__=='__main__':main()
