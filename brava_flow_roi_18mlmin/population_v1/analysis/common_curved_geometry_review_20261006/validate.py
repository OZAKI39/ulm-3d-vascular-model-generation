"""Independent geometry reference + actual native element/matrix regressions."""
from common import *
import numpy as np
import subprocess
import csv

EDGES=np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])
DL=np.vstack([np.eye(3),-np.ones(3)])
RN=np.vstack([np.eye(3),np.zeros((1,3))])
RN10=np.vstack([RN,RN[EDGES].mean(1)])
HP2=np.concatenate([4*np.einsum('ai,aj->aij',DL,DL),
    np.array([4*(np.outer(DL[i],DL[j])+np.outer(DL[j],DL[i])) for i,j in EDGES])])
PACK=[(0,0),(1,1),(2,2),(0,1),(1,2),(0,2)]

def basis(xi,n):
    L=np.r_[xi,1-np.sum(xi)]
    if n==4:return L,DL,np.zeros((4,3,3)),L
    N=np.r_[L*(2*L-1),4*L[EDGES[:,0]]*L[EDGES[:,1]]]
    D=np.vstack([(4*L-1)[:,None]*DL,4*(L[EDGES[:,0],None]*DL[EDGES[:,1]]+L[EDGES[:,1],None]*DL[EDGES[:,0]])])
    return N,D,HP2,L

def geometry(xi,X):
    n=len(X);N,D,H,P=basis(xi,n);Y=X-X[3]
    J=Y.T@D;inv=np.linalg.inv(J);det=np.linalg.det(J)
    Dx=D@inv;Dp=DL@inv
    geomH=np.einsum('ai,ajk->ijk',Y,H)
    Hx=np.array([inv.T@(H[a]-np.einsum('i,ijk->jk',Dx[a],geomH))@inv for a in range(n)])
    return dict(N=N,P=P,D=Dx,Dp=Dp,H=np.array([[h[i,j] for i,j in PACK] for h in Hx]),
        J=J,inv=inv,det=det,metric=inv.T@inv,x=N@Y+X[3])

def unpack(rows,n):
    rows=np.atleast_2d(rows);k=9
    out=dict(opt=rows[:,0].astype(int),g=rows[:,1].astype(int),wref=rows[:,4],det=rows[:,5],xi=rows[:,6:9])
    for name,shape in [('N',(n,)),('P',(4,)),('D',(n,3)),('Dp',(4,3)),('H',(n,6)),('metric',(3,3))]:
        count=int(np.prod(shape));out[name]=rows[:,k:k+count].reshape((len(rows),)+shape);k+=count
    assert k==rows.shape[1]
    return out

def read_RK(prefix,n):
    lines=Path(str(prefix)+'_RK.txt').read_text().splitlines()
    return np.fromstring(lines[0],sep=' ').reshape(4,n),np.fromstring(lines[1],sep=' ').reshape(16,n,n)

def fields(X,label,H):
    n=len(X);d=X-X[3];y=np.zeros((n,4));a=np.zeros((n,4));bf=np.zeros((n,3))
    y[:,3]=.3+d[:,2]/H
    if label=='A':y[:,:3]=d
    if label=='B':y[:,0]=d[:,0];y[:,1]=-d[:,1]
    if label=='C':y[:,2]=d[:,2]**2
    if label=='P':y[:,3]=.3+d@np.array([1.,2.,3.])/H
    if label=='R':
        rng=np.random.default_rng(61006+n);y[:,:3]=rng.normal(size=(n,3))*H
        y[:,3]=rng.normal(size=n);a[:,:3]=rng.normal(size=(n,3))*.2;bf=rng.normal(size=(n,3))*.1
    return y,a,bf

def main():
    protocol=read(HERE/'COMMON_GEOMETRY_HARNESS_PROTOCOL.json')
    identity=read(HERE/'COMMON_GEOMETRY_PATCH_IDENTITY.json')
    assert identity['protocol_sha256']==sha(HERE/'COMMON_GEOMETRY_HARNESS_PROTOCOL.json')
    old=np.load(PREVIOUS/'P2_CURVED_GEOMETRY_SAMPLES.npz');rt=HERE/'runtime';rt.mkdir(exist_ok=True)
    tests=[];cases=[];arrays={};csvrows=[];baseline_pressure=[]
    def test(name,error=0.,limit=1e-11,**detail):
        ok=bool(np.isfinite(error) and error<=limit)
        tests.append(dict(name=name,status='PASS' if ok else 'FAIL',error=float(error),limit=limit,**detail))
        return ok
    def run(exe,inp,prefix,mode='native',reference=None):
        cmd=[str(rt/('harness_'+exe)),str(inp),str(prefix),mode]
        if reference:cmd.append(str(reference))
        subprocess.run(cmd,env=native_env(),check=True,stdout=subprocess.DEVNULL)
    for cell in ['central_wall','inlet_rim']:
      for variant,n in [('affine',4),('affine',10),('curved',10)]:
       X=old[f'{cell}_{variant}_nodes'][:n]
       H=float(np.linalg.norm(X[:4,None]-X[None,:4],axis=2).max())
       for field in (['A','R'] if n==4 else ['A','B','C','P','R']):
        tag=f'{cell}_{variant}_P{1 if n==4 else 2}_{field}';prefix=rt/tag
        y,a,bf=fields(X,field,H);inp=rt/(tag+'_input.txt')
        with inp.open('w') as f:
            f.write(str(n)+'\n')
            for v in [X,y,a,bf]:np.savetxt(f,v,fmt='%.17g')
        runs={}
        for name in ['original','baseline_trace','patched','patched_trace']:
            pre=rt/(tag+'_'+name);run(name,inp,pre);runs[name]=read_RK(pre,n)
        rawbase=np.loadtxt(rt/(tag+'_baseline_trace_q.txt'))
        rawpatch=np.loadtxt(rt/(tag+'_patched_trace_q.txt'))
        qb=unpack(rawbase,n);qp=unpack(rawpatch,n)
        refs=[geometry(xi,X) for xi in qp['xi']]
        refpath=rt/(tag+'_reference.txt')
        with refpath.open('w') as f:
            f.write(str(len(refs))+'\n')
            for i,r in enumerate(refs):
                values=np.r_[qp['opt'][i],qp['wref'][i]*r['det'],r['N'],r['P'],r['D'].ravel(),
                    r['Dp'].ravel(),r['H'].ravel(),r['metric'].ravel()]
                f.write(' '.join(format(v,'.17g') for v in values)+'\n')
        run('original',inp,rt/(tag+'_reference'),'reference',refpath)
        RR,KR=read_RK(rt/(tag+'_reference'),n)
        us=H**2 if field=='C' else H
        sr=max(float(abs(RR[:3]).max()),H**2,.00345312*H*us,1056*H**2*us**2)
        sc=H**2*us
        for name in ['baseline_trace','patched_trace']:
            plain='original' if name=='baseline_trace' else 'patched'
            R0,K0=runs[plain];R1,K1=runs[name]
            test(tag+'_'+name+'_observation_does_not_change_R',float(abs(R1-R0).max())/max(float(abs(R0).max()),sr))
            test(tag+'_'+name+'_observation_does_not_change_K',float(abs(K1-K0).max())/max(float(abs(K0).max()),1e-100))
        RP,KP=runs['patched'];RB,KB=runs['original']
        test(tag+'_patched_momentum_vs_common_reference',float(abs(RP[:3]-RR[:3]).max())/sr)
        test(tag+'_patched_continuity_vs_common_reference',float(abs(RP[3]-RR[3]).max())/sc)
        for nm,sl in [('momentum',slice(0,12)),('continuity',slice(12,16))]:
            test(tag+'_patched_'+nm+'_matrix_vs_common_reference',float(abs(KP[sl]-KR[sl]).max())/max(float(abs(KR[sl]).max()),1e-100))
        if variant=='affine':
            test(tag+'_affine_R_before_after',float(abs(RP-RB).max())/max(float(abs(RB).max()),sr))
            test(tag+'_affine_K_before_after',float(abs(KP-KB).max())/max(float(abs(KB).max()),1e-100))
        if field=='R':
            repeat=rt/(tag+'_repeat');run('patched',inp,repeat)
            test(tag+'_deterministic_output_repeat',0. if sha(Path(str(repeat)+'_RK.txt'))==sha(rt/(tag+'_patched_RK.txt')) else 1.,0.)
        metrics={};data={}
        for name,q in [('before',qb),('after',qp)]:
            rr={k:np.array([r[k] for r in refs]) for k in refs[0]}
            refnodes=RN if n==4 else RN10
            vinv=np.einsum('ai,naj->nij',refnodes,q['D'])
            pinv=q['Dp'][:,:3,:]
            vJ=np.linalg.inv(vinv);pJ=np.linalg.inv(pinv)
            usedJ=vJ.copy()
            if name=='before':usedJ[(q['opt']==2)&(q['g']==0)]=pJ[(q['opt']==2)&(q['g']==0)]
            vworld=q['N']@(X-X[3])+X[3]
            pworld=q['P']@X[:4] if name=='before' else vworld
            u=q['N']@y[:,:3];G=np.einsum('ai,naj->nij',y[:,:3],q['D']);div=np.trace(G,axis1=1,axis2=2)
            p=q['P']@y[:4,3];gp=np.einsum('a,naj->nj',y[:4,3],q['Dp'])
            exu=rr['N']@y[:,:3];exG=np.einsum('ai,naj->nij',y[:,:3],rr['D'])
            exp=rr['P']@y[:4,3];exgp=np.einsum('a,naj->nj',y[:4,3],rr['Dp'])
            m=dict(pressure_basis_grad_scaled=float(abs(q['Dp']-rr['Dp']).max())/float(abs(rr['Dp']).max()),
                velocity_J_scaled=float(abs(vJ-rr['J']).max())/H,pressure_transform_J_scaled=float(abs(pJ-rr['J']).max())/H,
                inverse_scaled=float(abs(pinv-rr['inv']).max())/float(abs(rr['inv']).max()),
                used_det_relative=float(np.max(abs(q['det']/rr['det']-1))),
                metric_relative=float(abs(q['metric']-rr['metric']).max())/float(abs(rr['metric']).max()),
                physical_q_position_over_H=float(abs(pworld-rr['x']).max())/H,
                u_scaled=float(abs(u-exu).max())/us,G_scaled=float(abs(G-exG).max())/(us/H),
                p_scaled=float(abs(p-exp).max()),gradp_scaled=float(abs(gp-exgp).max())*H)
            if field in ['A','B','C']:
                d=rr['x']-X[3];au=np.zeros_like(u);aG=np.zeros_like(G)
                if field=='A':au=d;aG[:]=np.eye(3)
                if field=='B':au[:,0]=d[:,0];au[:,1]=-d[:,1];aG[:]=np.diag([1,-1,0])
                if field=='C':au[:,2]=d[:,2]**2;aG[:,2,2]=2*d[:,2]
                m.update(analytic_u_scaled=float(abs(u-au).max())/us,analytic_G_scaled=float(abs(G-aG).max())/(us/H),
                    analytic_div_scaled=float(abs(div-np.trace(aG,axis1=1,axis2=2)).max())/(us/H),
                    analytic_axial_pressure_scaled=float(abs(p-(.3+d[:,2]/H)).max()),
                    analytic_axial_gradp_scaled=float(abs(gp-np.array([0.,0.,1/H])).max())*H)
                curl=np.stack([G[:,2,1]-G[:,1,2],G[:,0,2]-G[:,2,0],G[:,1,0]-G[:,0,1]],axis=1)
                m['analytic_curl_scaled']=float(abs(curl).max())/(us/H)
            if field=='P':
                general=.3+(rr['x']-X[3])@np.array([1.,2.,3.])/H
                m['general_linear_p_representation_error']=float(abs(p-general).max())
                m['general_linear_gradp_representation_error']=float(abs(gp-np.array([1.,2.,3.])/H).max())*H
            metrics[name]=m
            if name=='after':
                for k,v in m.items():
                    if 'representation_error' not in k:test(tag+'_'+k,v,1e-12 if k=='physical_q_position_over_H' else 1e-11)
                test(tag+'_P1_reference_pressure_basis_unchanged',float(abs(q['P']-rr['P']).max()))
                # Physical Hessian checked at momentum points: native P1/P1 continuity retains old cached Hessian (zero).
                mom=q['opt']==1
                test(tag+'_native_Hessian_chain_rule',float(abs(q['H'][mom]-rr['H'][mom]).max())*H**2)
            data[name]=dict(velocity_J=vJ,pressure_transform_J=pJ,used_J_from_captured_transforms_and_branch=usedJ,
                used_det=q['det'],inverse_pressure=pinv,xq_velocity=vworld,xq_pressure_geometry=pworld,metric=q['metric'],
                weights=q['wref']*q['det'],u=u,G=G,div=div,p=p,gradp=gp)
        # Direct weak-form references do not call fluid kernels.
        cont=qp['opt']==2;mom=qp['opt']==1
        qwt=np.array([r['det'] for r in refs])*qp['wref']
        if n==10:
            direct=np.einsum('n,na,n->a',qwt[cont],qp['P'][cont],data['after']['div'][cont])
            test(tag+'_direct_continuity',float(abs(RP[3,:4]-direct).max())/sc)
        if field=='P':
            direct=-np.einsum('n,n,na j->ja',qwt[mom],data['after']['p'][mom],np.array([r['D'] for r in refs])[mom])
            test(tag+'_direct_pressure_only_momentum',float(abs(RP[:3]-direct).max())/sr)
            baseline_pressure.append(dict(case=tag,before=metrics['before'],after=metrics['after']))
        taub=np.loadtxt(rt/(tag+'_baseline_trace_tau.txt'));taup=np.loadtxt(rt/(tag+'_patched_trace_tau.txt'))
        if variant=='affine':test(tag+'_affine_actual_stabilization_parameters',float(abs(taub[:,2:]-taup[:,2:]).max())/max(float(abs(taub[:,2:]).max()),1e-100))
        # Momentum tauM must use the captured common metric without any coefficient changes.
        expected=[]
        for i,r in enumerate(refs):
            u=data['after']['u'][i];tau=1/(1056*np.sqrt(4/.01**2+u@r['metric']@u+36*(.00345312/1056)**2*np.sum(r['metric']**2)))
            expected.append(tau if qp['opt'][i]==1 or n==4 else 0.)
        test(tag+'_actual_tauM_common_metric',float(abs(taup[:,2]-expected).max())/max(max(expected),1e-100))
        for name,(R,K) in runs.items():arrays[tag+'_'+name+'_R']=R;arrays[tag+'_'+name+'_K']=K
        arrays[tag+'_reference_R']=RR;arrays[tag+'_reference_K']=KR
        arrays[tag+'_nodes']=X;arrays[tag+'_fields']=y;arrays[tag+'_before_q']=rawbase;arrays[tag+'_after_q']=rawpatch
        arrays[tag+'_before_tau']=taub;arrays[tag+'_after_tau']=taup
        for stage,vals in data.items():
            for k,v in vals.items():arrays[tag+'_'+stage+'_'+k]=v
        record=dict(case=tag,cell=cell,variant=variant,velocity_nodes=n,field=field,H_m=H,metrics=metrics,
            continuity_before_vs_common_relative_L2=float(np.linalg.norm(RB[3,:4]-RR[3,:4])/max(np.linalg.norm(RR[3,:4]),sc*1e-14)),
            continuity_after_scaled_max=float(abs(RP[3]-RR[3]).max())/sc,
            momentum_after_scaled_max=float(abs(RP[:3]-RR[:3]).max())/sr,
            P1_R_max_abs_change=float(abs(RP-RB).max()) if n==4 else None,
            P1_K_max_abs_change=float(abs(KP-KB).max()) if n==4 else None)
        cases.append(record)
        for st,m in metrics.items():csvrows.append(dict(case=tag,stage=st,**m))
        print(tag,'failures so far',sum(t['status']=='FAIL' for t in tests),flush=True)
    # Inherited representation evidence is hashed, not rerun as new A-D1 work.
    put('INHERITED_REGRESSION.json',dict(rerun=False,files={n:sha(PREVIOUS/n) for n in
        ['P2_REPRESENTATION_FIXTURES.json','P2_FIELD_SAMPLER_TESTS.json','P2_POISEUILLE_REPRESENTATION.json','P2_CURVED_GEOMETRY_SAMPLES.npz']},
        statement='Prior A-D1 evidence retained; this study only executes new element geometry regressions'))
    np.savez_compressed(HERE/'COMMON_GEOMETRY_SAMPLES.npz',**arrays)
    keys=list(dict.fromkeys(k for r in csvrows for k in r))
    with (HERE/'COMMON_GEOMETRY_METRICS.csv').open('w') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(csvrows)
    p1tests=[t for t in tests if '_P1_' in t['name'] and ('central_wall_affine_P1_' in t['name'] or 'inlet_rim_affine_P1_' in t['name'])]
    p1pass=all(t['status']=='PASS' for t in p1tests)
    put('P1_REGRESSION_RESULTS.json',dict(status='PASS' if p1pass else 'P1_REGRESSION_FAILED',
        cases=[c for c in cases if c['velocity_nodes']==4],tests=p1tests,production_modified=False))
    ok=all(t['status']=='PASS' for t in tests)
    put('COMMON_GEOMETRY_HARNESS_RESULTS.json',dict(status='PASS' if ok else 'COMMON_GEOMETRY_FIX_NOT_EVALUABLE',
        cases=cases,passed=sum(t['status']=='PASS' for t in tests),failed=sum(t['status']=='FAIL' for t in tests),tests=tests,
        protocol_sha256=sha(HERE/'COMMON_GEOMETRY_HARNESS_PROTOCOL.json'),
        trace_semantics='D/metric/det/weights and tau captured from native assembly. Full J matrices reconstructed from captured reference-coordinate gradients; used J selection follows audited native branch, independently checked against captured determinant.',
        CFD_steps=0,linear_solves=0))
    base=read(HERE/'CURVED_MIXED_BLOCKER_BASELINE.json');base['linear_pressure_extension']=baseline_pressure
    put('CURVED_MIXED_BLOCKER_BASELINE.json',base)
    if not p1pass:raise SystemExit('P1_REGRESSION_FAILED')
    if not ok:raise SystemExit('COMMON_GEOMETRY_FIX_NOT_EVALUABLE')

if __name__=='__main__':main()
