"""Unmodified function bodies extracted from archived physics.py; no legacy workflow imports."""
from pathlib import Path
import math
import numpy as np

def read_off(path):
    lines=[x.strip() for x in Path(path).read_text().splitlines() if x.strip() and not x.startswith('#')]
    if lines[0]!='OFF':raise ValueError('NOT_OFF')
    nv,nf,_=map(int,lines[1].split());v=np.array([list(map(float,s.split())) for s in lines[2:2+nv]])
    raw=[list(map(int,s.split())) for s in lines[2+nv:2+nv+nf]]
    if any(s[0]!=3 for s in raw):raise ValueError('TRIANGLES_REQUIRED')
    return v,np.array([s[1:] for s in raw],int)

def write_off(path,v,f):
    with Path(path).open('x') as out:
        out.write(f'OFF\n{len(v)} {len(f)} 0\n')
        np.savetxt(out,v,fmt='%.17g');np.savetxt(out,np.column_stack([np.full(len(f),3),f]),fmt='%d')

def order_vertices(ids,positions,n):
    ids=np.asarray(ids,dtype=np.int64)
    if len(ids)!=n or len(np.unique(ids))!=n:raise ValueError('MISSING_OR_DUPLICATE_VERTEX')
    order=np.argsort(ids)
    if not np.array_equal(ids[order]-ids.min(),np.arange(n)):raise ValueError('NONCONTIGUOUS_MEMBRANE_IDS')
    return np.asarray(positions)[order]

def geometry(v,f,degenerate_ratio=1.05):
    v=np.asarray(v,float);f=np.asarray(f,int)
    if not np.isfinite(v).all():raise ValueError('NONFINITE_MESH')
    q=v-v.mean(axis=0);tri=q[f];cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(cross,axis=1).sum()/2
    signed_volume=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2])).sum()/6;volume=abs(signed_volume)
    edges={}
    for face in f:
        for i,j in zip(face,np.roll(face,-1)):edges.setdefault(tuple(sorted((int(i),int(j)))),[]).append(1 if i<j else -1)
    closed=all(len(e)==2 and sum(e)==0 for e in edges.values())
    # Vertex-weighted second moment in the fixed X-Z shear plane; no fitted rotation.
    covariance=q[:,[0,2]].T@q[:,[0,2]]/len(q);w,axes=np.linalg.eigh(covariance);B,L=np.sqrt(np.maximum(w,0))*2
    ratio=L/max(B,1e-30);D=(L-B)/(L+B);theta=math.degrees(math.atan2(axes[1,1],axes[0,1]));theta=(theta+90)%180-90
    return dict(area=float(area),volume=float(volume),signed_volume=float(signed_volume),reduced_volume=float(6*math.sqrt(math.pi)*volume/area**1.5),a=float(math.sqrt(area/(4*math.pi))),
                center=v.mean(axis=0).tolist(),D=float(D),theta_deg=theta if ratio>=degenerate_ratio else None,
                angle_degenerate=bool(ratio<degenerate_ratio),axis_ratio=float(ratio),closed=closed,vertices=len(v),faces=len(f),
                euler_characteristic=len(v)-len(edges)+len(f),minimum_face_area=float(np.linalg.norm(cross,axis=1).min()/2),
                bounds=[v.min(axis=0).tolist(),v.max(axis=0).tolist()],deformation_definition='vertex_covariance_XZ_axes_sqrt_eigenvalues')
