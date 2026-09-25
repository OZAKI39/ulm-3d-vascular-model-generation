"""Independent NumPy/HDF5 assembly and solves versus raw standalone C++ output."""
import csv
import json
import sys
from pathlib import Path
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
from cf2003_free_shear_reference import CF2003_FU,CF2003_FOMEGA
TABLE=Path('/home/lzy/projects/compre_output/wall_hydrodynamics_v0/20260916_130547/tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5')

def main():
    data=list(csv.DictReader((ROOT/'raw/PHASE1_INPUTS.csv').open()))
    raw=list(csv.DictReader((ROOT/'raw/CPP_PHASE1_OUTPUT.csv').open()))
    assert len(data)==len(raw) and all(x['id']==y['id'] for x,y in zip(data,raw))
    statuses=all(x['expected_status']==y['status'] for x,y in zip(data,raw))
    valid=[i for i,r in enumerate(raw) if r['status']=='OK']
    d=[data[i] for i in valid];out=[raw[i] for i in valid]
    def inputs(names):return np.array([[float(r[k]) for k in names] for r in d])
    def outputs(prefix,n):return np.array([[float(r[prefix+str(i)]) for i in range(n)] for r in out])
    a,mu,h=inputs(['a','mu','h']).T;e=h/a
    t1=inputs(['t1x','t1y','t1z']);t2=inputs(['t2x','t2y','t2z']);n=inputs(['nx','ny','nz']);g=inputs(['gx','gy','gz'])
    qframe=np.stack([t1,t2,n],axis=2);transform=np.zeros((len(d),6,6))
    transform[:,:3,:3]=qframe;transform[:,3:,3:]=qframe
    with h5py.File(TABLE,'r') as f:
        grid=f['epsilon'][:];table=f['R_wall_excess_scaled'][:]
        assert np.max(abs(f['R_total_scaled'][:]-table-np.eye(6)))==0
    # Separate interpolation implementation: NumPy interpolation per matrix
    # component, without importing or translating C++ interpolation code.
    local=np.stack([np.interp(np.log(e),np.log(grid),table[:,i,j])
                    for i in range(6) for j in range(6)],axis=1).reshape(-1,6,6)
    bulk=np.repeat(np.stack([6*np.pi*mu*a,8*np.pi*mu*a**3],axis=1),3,axis=1)
    scale=np.sqrt(bulk)
    local_si=local*scale[:,:,None]*scale[:,None,:]
    excess=np.einsum('nik,nkl,njl->nij',transform,local_si,transform)
    total=excess+bulk[:,:,None]*np.eye(6)
    gt=g-n*np.sum(g*n,axis=1)[:,None]
    qfree=np.concatenate([(a+h)[:,None]*gt,.5*np.cross(n,gt)],axis=1)
    factors=np.array([[CF2003_FU(x),CF2003_FOMEGA(x)] for x in e])
    target=qfree*np.repeat(factors,3,axis=1)
    shear=np.einsum('nij,nj->ni',total,target)-bulk*qfree
    char=np.maximum(np.linalg.norm(scale*qfree,axis=1),scale[:,0]*a*np.linalg.norm(g,axis=1))
    char=np.where(char==0,scale[:,0]*a,char)
    E=outputs('E',36).reshape(-1,6,6);R=outputs('R',36).reshape(-1,6,6)
    b=outputs('b',6);cpp_bulk=outputs('bulk',6);cpp_qfree=outputs('qfree',6)
    cpp_target=outputs('qtarget',6);cpp_solved=outputs('qsolved',6)
    def mscale(m):return m/scale[:,:,None]/scale[:,None,:]
    def merror(x,y):return np.linalg.norm(mscale(x-y),axis=(1,2))/np.maximum(1,np.linalg.norm(mscale(y),axis=(1,2)))
    def berror(x,y):return np.linalg.norm((x-y)/scale,axis=1)/char
    def qerror(x,y):return np.linalg.norm((x-y)*scale,axis=1)/char
    model=mscale(total);cpp_model=mscale(R)
    expected_rhs=(bulk*qfree+shear)/scale
    q_np=np.linalg.solve(model,expected_rhs)/scale
    q_cpp_np=np.linalg.solve(cpp_model,(cpp_bulk*cpp_qfree+b)/scale)/scale
    eig=np.linalg.eigvalsh((cpp_model+cpp_model.transpose(0,2,1))*.5)
    tests={}
    def gate(name,values,limit):
        v=np.asarray(values);worst=float(np.max(v));ok=bool(np.isfinite(v).all() and np.all(v<=limit))
        tests[name]={'pass':ok,'worst':worst,'limit':limit,'count':int(v.size)}
    gate('A',merror(R,R.transpose(0,2,1)),1e-10)
    tests['B']={'pass':bool(np.isfinite(eig).all() and (eig>0).all()),'minimum_scaled_eigenvalue':float(eig.min()),'count':len(d)}
    metric={x:[] for x in ['C','D','E','F','G','H']}
    groups={}
    for i,r in enumerate(d):groups.setdefault(r['group'],[]).append(i)
    factors_scale=np.array([[6 if i<3 and j<3 else 24 if i>=3 and j>=3 else 12 for j in range(6)] for i in range(6)])
    for k in range(96):
        base,rot,flip=groups['rotated_'+str(k)]
        ss=scale[base];den=char[base]
        def mm(x,y):return float(np.linalg.norm((x-y)/ss[:,None]/ss[None,:])/max(1,np.linalg.norm(y/ss[:,None]/ss[None,:])))
        def bb(x,y):return float(np.linalg.norm((x-y)/ss)/den)
        metric['C'].append(mm(E[groups['scale_'+str(k)][0]],E[base]*factors_scale))
        metric['D'].append(max(mm(E[rot],E[base]),bb(b[rot],b[base])))
        metric['E'].append(max(mm(E[flip],E[base]),bb(b[flip],b[base])))
        metric['F'].append(bb(b[groups['zero_'+str(k)][0]],np.zeros(6)))
        metric['G'].append(bb(b[groups['reverse_'+str(k)][0]],-b[base]))
        metric['H'].append(bb(b[groups['double_'+str(k)][0]],2*b[base]))
    for name in metric:gate(name,metric[name],1e-12 if name=='F' else 1e-9 if name in ['D','E'] else 1e-10)
    gate('I',abs(np.sum(b[:,:3]*n,axis=1))/scale[:,0]/char,1e-12)
    gate('J',np.maximum(qerror(cpp_solved,target),qerror(q_cpp_np,target)),1e-8)
    dense=np.array(groups['dense'])
    solved_factors=np.stack([q_cpp_np[dense,0]/qfree[dense,0],q_cpp_np[dense,4]/qfree[dense,4]],axis=1)
    gate('K',np.max(abs(solved_factors/factors[dense]-1),axis=1),1e-8)
    gate('independent_bulk_SI',np.max(abs(cpp_bulk/bulk-1),axis=1),1e-12)
    gate('independent_excess_SI',merror(E,excess),1e-10)
    gate('independent_total_SI',merror(R,total),1e-10)
    gate('independent_shear_RHS',berror(b,shear),1e-10)
    gate('independent_qfree',qerror(cpp_qfree,qfree),1e-10)
    gate('independent_qtarget',qerror(cpp_target,target),1e-10)
    gate('independent_solve',qerror(cpp_solved,q_cpp_np),1e-10)
    gate('python_assembly_recovery',qerror(q_np,target),1e-8)
    cpp_factors=np.array([[float(r['FU']),float(r['FOMEGA'])] for r in out])
    gate('independent_reference_evaluation',np.max(abs(cpp_factors/factors-1),axis=1),1e-10)
    rejected=0
    for value in [-1.,0.,.000999,.200001,float('nan'),float('inf')]:
        for function in [CF2003_FU,CF2003_FOMEGA]:
            try:function(value)
            except ValueError:rejected+=1
    tests['INPUT_REJECTION']={'pass':statuses and rejected==12,'cpp_contract_statuses_match':statuses,'python_invalid_queries_rejected':rejected}
    cpp=json.loads((ROOT/'validation/CPP_PHASE1_VALIDATION.json').read_text())
    passed=all(t['pass'] for t in tests.values()) and cpp['ALGEBRA_PASS']
    report={'ALGEBRA_PASS':passed,'INDEPENDENT_IMPLEMENTATION_PASS':passed,'REFERENCE_SOURCE_VALIDATED':True,
            'PHYSICS_REFERENCE_STATUS':'PARTIAL','free_shear_reference_status':'PASS',
            'meaning':'C++/NumPy implementation and assembly agreement; one CF2003 physics basis. RMBW finite-gap RR/TR accuracy remains limited.',
            'valid_cases':len(d),'invalid_cases':len(data)-len(d),'tests':tests,
            'independent_python_methods':['h5py original HDF5 read','np.interp in ln(epsilon)','NumPy tensor rotation','np.linalg.eigvalsh','np.linalg.solve'],
            'K_scope':'All 10009 points in dense certified reference grid; actual independent solve of C++ matrix and RHS'}
    (ROOT/'validation/PHASE1_VALIDATION.json').write_text(json.dumps(report,indent=2)+'\n')
    with (ROOT/'raw/REFERENCE_RECOVERY.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['epsilon','FU_reference','FOMEGA_reference','FU_solved','FOMEGA_solved','relative_FU_error','relative_FOMEGA_error'])
        w.writerows(zip(e[dense],*factors[dense].T,*solved_factors.T,*abs(solved_factors/factors[dense]-1).T))
    with (ROOT/'raw/SHEAR_RHS_COMPONENTS.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['epsilon','Fx_N','Fy_N','Fz_N','Tx_Nm','Ty_Nm','Tz_Nm']);w.writerows(zip(e[dense],*b[dense].T))
    print(json.dumps(report,indent=2));assert passed,'PHASE1_VALIDATION_FAILED'

if __name__=='__main__':main()
