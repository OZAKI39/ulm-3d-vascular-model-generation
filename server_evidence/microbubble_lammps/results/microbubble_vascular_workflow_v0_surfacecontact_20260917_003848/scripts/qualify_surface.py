#!/usr/bin/env python3
"""Independent topology/QP checks and actual C++ constrained trajectories on immutable test STL."""
import sys,json,struct,subprocess,itertools,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import Delaunay,ConvexHull
from surface_reference import SurfaceReference,project_velocity
S=Path(__file__).resolve().parents[1];F=S/'raw/synthetic';F.mkdir(exist_ok=True)
EXE=Path(sys.argv[1]) if len(sys.argv)>1 else Path('/home/lzy/projects/compre_output/work/surfacecontact_20260917_003848/surface_contact_probe')

def mesh(name,v,f):
    v=np.asarray(v,dtype=np.float32).astype(float);f=np.asarray(f,int)
    with (F/(name+'.stl')).open('wb') as o:
        o.write(b'IMMUTABLE SYNTHETIC WALL'.ljust(80,b' '));o.write(struct.pack('<I',len(f)))
        for tri in v[f]:o.write(struct.pack('<12fH',*([0.]*3+tri.ravel().tolist()),0))
    np.savez(F/(name+'.npz'),vertices=v,faces=f)
    return SurfaceReference(v,f)

def plane(name,seed):
    rng=np.random.default_rng(seed);xy=np.vstack([[[a,b] for a in [-5.,0,5.] for b in [-5.,0,5.]],rng.uniform(-4.8,4.8,(40,2))]);f=Delaunay(xy).simplices;v=np.column_stack([xy,np.zeros(len(xy))])*1e-6
    for a in f:
        if np.cross(v[a[1]]-v[a[0]],v[a[2]]-v[a[0]])[2]<0:a[[1,2]]=a[[2,1]]
    return mesh(name,v,f)

def curved(name,sphere=False):
    if sphere:
        phi=np.arange(96)*2*np.pi/96;theta=np.arange(1,48)*np.pi/48;v=np.array([[8e-6*np.sin(t)*np.cos(p),8e-6*np.sin(t)*np.sin(p),8e-6*np.cos(t)] for t in theta for p in phi]+[[0,0,8e-6],[0,0,-8e-6]])
        f=ConvexHull(v).simplices
        for a in f:
            if np.dot(np.cross(v[a[1]]-v[a[0]],v[a[2]]-v[a[0]]),v[a].mean(axis=0))<0:a[[1,2]]=a[[2,1]]
    else:
        v=np.array([[8e-6*np.cos(p),8e-6*np.sin(p),z] for z in [-5e-6,-2.5e-6,0,2.5e-6,5e-6] for p in np.arange(96)*2*np.pi/96]);f=[]
        for j in range(4):
            for k in range(96):
                a=j*96+k;b=j*96+(k+1)%96;c=(j+1)*96+k;d=(j+1)*96+(k+1)%96;f.extend([[a,b,d],[a,d,c]])
    return mesh(name,v,f)

def wedge(alpha):
    h=np.radians(alpha/2);v=[[0,0,-5e-6],[0,0,5e-6],[10e-6*np.cos(h),10e-6*np.sin(h),-5e-6],[10e-6*np.cos(h),10e-6*np.sin(h),5e-6],[10e-6*np.cos(h),-10e-6*np.sin(h),-5e-6],[10e-6*np.cos(h),-10e-6*np.sin(h),5e-6]]
    return mesh('wedge'+str(alpha),v,[[0,2,3],[0,3,1],[0,1,5],[0,5,4]])

class Probe:
    def __init__(self,name,production=False):
        cmd=[str(EXE),str(S/'geometry/WALL_ONLY.stl') if production else str(F/(name+'.stl'))]
        if production:cmd.append(str(S/'geometry/WALL_TRIANGLE_REGIONS.txt'))
        self.p=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
    def ask(self,x,v,factor=1,theta=15,contact=None):
        args=['Q' if contact is None else 'C',*x,*v,factor,theta]
        if contact is not None:args+=list(contact)
        self.p.stdin.write(' '.join(map(str,args))+'\n');self.p.stdin.flush();line=self.p.stdout.readline()
        if not line:raise RuntimeError('C++ probe stopped')
        return json.loads(line)
    def close(self):self.p.stdin.close();assert self.p.wait()==0

# The probe must flush per response for interactive trajectory stepping.
metrics=dict(query_count=0,max_normal_error=0.,max_velocity_error=0.,max_objective_error=0.,max_distance_error=0.,features={},modes={},failures=[])
records=[]
def compare(ref,probe,x,v,factor=1,theta=15,label=''):
    py=ref.query(x,factor,theta);cpp=probe.ask(x,v,factor,theta);q=cpp['query'];metrics['query_count']+=1
    assert py['status']==q['status'],(label,py['status'],q['status'])
    assert [c['id'] for c in py['candidates']]==[c['id'] for c in q['candidates']],label
    assert len(py['clusters'])==len(q['clusters']),label
    de=abs(py['minimum_distance']-q['minimum_distance']);metrics['max_distance_error']=max(de,metrics['max_distance_error']);assert de<2e-19
    for a,b in zip(py['candidates'],q['candidates']):
        assert a['feature']==b['feature'] and a['support']==b['support'],(label,a,b)
        metrics['features'][a['feature']]=metrics['features'].get(a['feature'],0)+1
        for j in b['adjacency']:
            key=tuple(sorted((a['id'],j['id'])));assert abs(ref.dihedral[key]-j['dihedral'])<1e-5
    for a,b in zip(py['clusters'],q['clusters']):
        assert a['support']==b['support'] and a['id']==b['id'] and a['candidate_ids']==b['candidate_ids'],label
        assert np.linalg.norm(np.array(a['surface_normal'])-b['surface_normal'])<2e-10
        assert abs(np.linalg.norm(b['surface_normal'])-1)<1e-12 and np.dot(b['surface_normal'],np.asarray(x)-b['point'])>0
        e=np.linalg.norm(np.array(a['normal'])-b['normal']);metrics['max_normal_error']=max(e,metrics['max_normal_error']);assert e<2e-10
        assert abs(np.linalg.norm(b['normal'])-1)<1e-12 and np.dot(b['normal'],np.asarray(x)-b['point'])>0
    if py['clusters']:
        # Separate geometric rounding at 10-nm seam probes from KKT accuracy.
        # C++ normals have already been independently reconstructed above.
        N=np.array([c['normal'] for c in py['clusters']]);C=np.array([c['normal'] for c in q['clusters']]);obj,used,active,lam=project_velocity(v,N)
        ue=np.linalg.norm(used-cpp['used']);metrics['max_velocity_error']=max(ue,metrics['max_velocity_error']);assert ue<4e-10*max(np.linalg.norm(v),1e-20),(label,ue)
        exact_obj,exact_used,exact_active,_=project_velocity(v,C)
        assert np.linalg.norm(exact_used-cpp['used'])<2e-12*max(np.linalg.norm(v),1e-20)
        oe=abs(exact_obj-cpp['objective']);metrics['max_objective_error']=max(oe,metrics['max_objective_error']);assert oe<2e-12*max(np.dot(v,v),1e-40)
        assert np.min(C@cpp['used'])>=-2e-12*max(np.linalg.norm(v),1e-20)
        py_tight=np.flatnonzero(np.abs(C@exact_used)<2e-12*max(np.linalg.norm(v),1e-20));cpp_tight=np.flatnonzero(np.abs(C@cpp['used'])<2e-12*max(np.linalg.norm(v),1e-20));assert np.array_equal(py_tight,cpp_tight)
    metrics['modes'][q['status']]=metrics['modes'].get(q['status'],0)+1
    records.append(dict(label=label,x=np.asarray(x).tolist(),v=np.asarray(v).tolist(),factor=factor,theta=theta,cpp=cpp))
    return cpp

def trajectory(name,ref,dtmax,factor=1,theta=15,T=.01,plane_case=False):
    probe=Probe(name);radius=1e-6;margin=1e-10
    if plane_case:x=np.array([-1e-6,0.,radius+margin+1e-17])
    else:
        direction=np.array([np.cos(.021),np.sin(.021),0.]);lo=0.;hi=8e-6
        for _ in range(70):
            a=(lo+hi)/2
            if ref.query(direction*a)['minimum_distance']>radius+margin+1e-17:lo=a
            else:hi=a
        x=direction*lo
    t=0.;accepted=0;rejected=0;offsets=[];rows=[];faces=set();mingap=1.
    def raw(pos):
        if plane_case:return np.array([1e-4,0.,-2e-5])
        radial=pos.copy();radial[2]=0;radial/=np.linalg.norm(radial);return 1e-4*np.array([-radial[1],radial[0],0])+2e-5*radial
    while t<T-1e-15:
        dt=min(dtmax,T-t);allow=False
        for attempt in range(40):
            first=probe.ask(x,raw(x),factor,theta,[radius,margin,dt/2,int(allow),*x]);mid=np.array(first['endpoint'])
            second=probe.ask(mid,raw(mid),factor,theta,[radius,margin,dt,int(allow),*x]) if first['safe'] else None
            if second is not None and second['safe']:break
            dt/=2;rejected+=1;allow=True
            if dt<1e-14:raise AssertionError((name,'TIME_REFINEMENT_EXHAUSTED',t,first,second))
        else:raise AssertionError('retry cap')
        for r in [first,second]:
            assert r['swept_gap']>=margin-3e-20;assert len(r['query']['clusters'])==1
            faces.update(c['id'] for c in r['query']['candidates']);mingap=min(mingap,r['swept_gap'])
            if r['correction']:assert r['correction']<=margin/4;offsets.append(r['correction'])
            N=[c['normal'] for c in r['query']['clusters']]
            expected=project_velocity(raw(x if r is first else mid),N)[1] if not r['raw_safe'] else raw(x if r is first else mid)
            assert np.linalg.norm(expected-r['used'])<2e-16
        x=np.array(second['endpoint']);t+=dt;accepted+=1;rows.append([t,*x,dt,len(second['query']['clusters']),second['correction']])
        if accepted>100000:raise AssertionError('trajectory cap')
    probe.close();return dict(name=name,dt=dtmax,factor=factor,theta=theta,time=t,accepted=accepted,rejected=rejected,endpoint=x.tolist(),minimum_gap=mingap,faces_visited=len(faces),max_correction=max(offsets,default=0),total_correction=sum(offsets),trajectory=rows)

def main():
    refs={f'plane{k}':plane(f'plane{k}',k+37) for k in range(3)};refs['cylinder']=curved('cylinder');refs['sphere']=curved('sphere',True)
    for a in [30,60,90,120,150]:refs['wedge'+str(a)]=wedge(a)
    v=np.array([[0,0,0],[5,0,0],[0,5,0],[0,0,5],[5,5,0],[0,5,5],[5,0,5]],float)*1e-6;refs['corner']=mesh('corner',v,[[0,2,4],[0,4,1],[0,1,6],[0,6,3],[0,3,5],[0,5,2]])
    refs['nonmanifold']=mesh('nonmanifold',np.array([[0,0,0],[1,0,0],[0,1,0],[0,-1,0],[0,0,1]])*1e-6,[[0,1,2],[1,0,3],[0,1,4]])
    for name,ref in list(refs.items()):
        data=np.load(F/(name+'.npz'));refs[name+'_reordered']=mesh(name+'_reordered',data['vertices'],data['faces'][np.random.default_rng(42).permutation(len(data['faces']))])
    a=1.0001e-6;tests={}
    for name,ref in refs.items():
        if name.endswith('_reordered'):continue
        probe=Probe(name);reorder=Probe(name+'_reordered');points=[]
        if name.startswith('plane'):
            points=[np.array([x,0,a]) for x in np.linspace(-1e-6,1e-6,31)]+[ref.t[0].mean(axis=0)+[0,0,a],ref.t[0,:2].mean(axis=0)+[0,0,a],ref.t[0,0]+[0,0,a]]
        elif name in ['cylinder','sphere']:
            points=[ref.t[k].mean(axis=0)-a*ref.normal[k] for k in range(0,len(ref.faces),max(1,len(ref.faces)//30))]
            for edge,inc in list(ref.edges.items())[::max(1,len(ref.edges)//20)]:
                if len(inc)==2:
                    n=np.sum(ref.normal[inc],axis=0);n/=np.linalg.norm(n);points.append(ref.vertices[list(edge)].mean(axis=0)+a*n)
            for vertex in range(0,len(ref.vertices),max(1,len(ref.vertices)//20)):
                n=ref.vertices[vertex].copy()
                if name=='cylinder':n[2]=0
                n/=np.linalg.norm(n);points.append(ref.vertices[vertex]+a*n)
        elif name.startswith('wedge'):points=[np.array([a/np.sin(np.radians(int(name[5:])/2)),0.,0.])]
        elif name=='corner':points=[np.array([a,a,a])]
        else:points=[np.array([.5e-6,0.,.2e-6])]
        for x in points:
            for factor,theta in itertools.product([.5,1,2],[10,15,25]):
                for v in [np.array([-1.1, -.7,-.3])*1e-4,np.array([.4,.8,1.2])*1e-4]:
                    out=compare(ref,probe,x,v,factor,theta,name);other=compare(refs[name+'_reordered'],reorder,x,v,factor,theta,name+'_reordered')
                    assert out.get('used')==other.get('used') and out['query']['clusters']==other['query']['clusters'],'REORDER_VARIANCE'
                    n=len(out['query']['clusters'])
                    if name.startswith('plane') or name in ['cylinder','sphere']:assert n==1,(name,n)
                    if name.startswith('wedge'):assert n==2,(name,n)
                    if name=='corner':assert n==3
                    if name=='nonmanifold':assert out['query']['status']=='FAIL_NONMANIFOLD'
        probe.close();reorder.close();tests[name]=dict(points=len(points),status='PASS')
    print('QUERY_CHECKS_PASS',metrics,flush=True)
    trajectories=[]
    for name in ['plane0','plane1','plane2','plane0_reordered','cylinder','cylinder_reordered','sphere','sphere_reordered']:
        for dt in ([1e-4] if name.startswith('plane') or name.endswith('reordered') else [1e-4,5e-5,2.5e-5,1.25e-5]):
            tr=trajectory(name,refs[name],dt,plane_case=name.startswith('plane'));trajectories.append(tr);print('TRAJECTORY',name,dt,tr['accepted'],tr['rejected'],tr['faces_visited'],flush=True)
    for base in ['plane0','cylinder','sphere']:
        a=next(t for t in trajectories if t['name']==base);b=next(t for t in trajectories if t['name']==base+'_reordered');assert a['trajectory']==b['trajectory']
    plane_ends=[t['endpoint'] for t in trajectories if t['name'].startswith('plane')];assert np.max(np.ptp(plane_ends,axis=0))<1e-18
    for name in ['cylinder','sphere']:
        ts=[t for t in trajectories if t['name']==name];fine=np.array(ts[-1]['endpoint']);errors=[float(np.linalg.norm(np.array(t['endpoint'])-fine)) for t in ts[:-1]]
        assert all(b<=a+1e-15 for a,b in zip(errors,errors[1:])),(name,errors)
        for t in ts:assert t['faces_visited']>=2
        tests[name]['endpoint_error_against_finest_dt']=errors
        tests[name]['finest_dt']=ts[-1]['dt']
    sensitivities=[]
    for name in ['plane0','cylinder','sphere']:
        for factor,theta in itertools.product([.5,1,2],[10,15,25]):
            tr=trajectory(name,refs[name],1e-4,factor,theta,T=.002,plane_case=name.startswith('plane'));sensitivities.append(tr)
        same=[t['endpoint'] for t in sensitivities if t['name']==name];assert np.max(np.ptp(same,axis=0))<1e-15
    (S/'validation/SYNTHETIC_SURFACE_QUALIFICATION.json').write_text(json.dumps(dict(status='PASS',metrics=metrics,tests=tests,trajectories=[{k:v for k,v in t.items() if k!='trajectory'} for t in trajectories],sensitivity=[{k:v for k,v in t.items() if k!='trajectory'} for t in sensitivities]),indent=2)+'\n')
    (S/'raw/SYNTHETIC_QUERY_COMPARISON.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in records))
    (S/'raw/SYNTHETIC_TRAJECTORIES.json').write_text(json.dumps(trajectories+sensitivities)+'\n');print('SYNTHETIC_PASS',flush=True)
if __name__=='__main__':
    try:main()
    except Exception as e:
        (S/'validation/SYNTHETIC_SURFACE_FAILURE.json').write_text(json.dumps(dict(error=repr(e),metrics=metrics,last_queries=records[-4:]),indent=2)+'\n');raise
