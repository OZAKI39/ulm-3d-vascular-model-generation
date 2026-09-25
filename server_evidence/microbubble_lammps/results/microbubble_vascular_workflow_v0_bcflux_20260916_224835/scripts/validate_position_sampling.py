from pathlib import Path
import json,csv
import numpy as np
from flow_geometry import Geometry,FrozenSampler,refine,quadrature
S=Path(__file__).resolve().parents[1];g=Geometry(S/'geometry/GEOMETRY_ARRAYS.npz');s=FrozenSampler(S/'fields/FROZEN_FLOW_FIELD_V0.h5');section=dict(np.load(S/'geometry/INJECTION_SECTION.npz'))
tri=section['triangles']
for i in range(4):tri=refine(tri)
x,w=quadrature(tri);status,u=s.query(x);un=-u@section['normal'];radius=9.683592065545495e-7
valid=(status==0)&(un>0);validids=np.flatnonzero(valid)
wall=np.array([g.distance(p)[0] for p in x[validids]]);valid[validids[wall<radius+s.dx]]=False
x=x[valid];w=w[valid];un=un[valid]
with (S/'raw/POSITION_SAMPLES.csv').open() as f:records=list(csv.DictReader(f))
results=[]
for mode in ['VALID_AREA_UNIFORM','VALID_SUPPORT_FLUX_WEIGHTED']:
    rows=[r for r in records if r['mode']==mode];p=np.array([[float(r[k]) for k in ['x_m','y_m','z_m']] for r in rows]);a=np.array([float(r['radius_m']) for r in rows]);st,vel=s.query(p);assert (st==0).all();v=-vel@section['normal'];assert (v>0).all();assert all(g.distance(q)[0]>=rr+s.dx-2e-14 for q,rr in zip(p,a));assert all(r['disclaimer']=='NOT EXPERIMENTAL CONCENTRATION' for r in rows)
    weight=w*(un if mode.endswith('WEIGHTED') else 1);weight/=weight.sum();ks={}
    # Three one-dimensional marginal CDF comparisons with independent numerical
    # quadrature; no renormalization of unknown full-plane flow occurs here.
    for label,pv,qv in [('axis1',(p-section['center'])@section['axis1'],(x-section['center'])@section['axis1']),('axis2',(p-section['center'])@section['axis2'],(x-section['center'])@section['axis2']),('inward_velocity',v,un)]:
        order=np.argsort(qv);qv=qv[order];cdf=np.cumsum(weight[order]);ps=np.sort(pv);idx=np.searchsorted(qv,ps,side='right')-1;F=np.where(idx<0,0,cdf[np.maximum(idx,0)]);e=np.maximum(abs(F-np.arange(len(ps))/len(ps)),abs(F-(np.arange(len(ps))+1)/len(ps)));ks[label]=float(e.max());assert ks[label]<.05,(mode,label,ks)
    results.append(dict(mode=mode,count=len(rows),marginal_KS=ks,gate=.05,mean_positive_velocity_m_s=float(v.mean()),quadrature_expected_mean_m_s=float(weight@un)))
(S/'validation/POSITION_DISTRIBUTION_COMPARISON.json').write_text(json.dumps(dict(status='PASS',scope='conditional center distribution at fixed diagnostic median radius; no active congestion; not full physical inlet sampling',position_status='WORKFLOW_V0_APPROXIMATION',samples=results,disclaimer='NOT EXPERIMENTAL CONCENTRATION'),indent=2)+'\n');print(json.dumps(results,indent=2))
