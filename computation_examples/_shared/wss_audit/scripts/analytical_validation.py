"""Manufactured nodal fields through the exact historical production WSS functions.
This does NOT solve Navier-Stokes and is NOT vascular mesh convergence.
"""
from pathlib import Path
import sys, importlib.util, json, csv
import numpy as np
from scipy.spatial import Delaunay
A=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
sys.path.insert(0,str(A/'evidence/source'))
spec=importlib.util.spec_from_file_location('actual_wss',A/'evidence/source/compute_field_diagnostics.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
MU=.00345312

def cylinder(nr,nt,nz,R=4e-6,L=40e-6):
    xy=[[0.,0.]]
    for r in np.linspace(R/nr,R,nr):
        xy.extend(np.column_stack([r*np.cos(np.arange(nt)*2*np.pi/nt),r*np.sin(np.arange(nt)*2*np.pi/nt)]))
    xy=np.asarray(xy);tri=np.sort(Delaunay(xy/R).simplices,axis=1);N=len(xy)
    x=np.vstack([np.column_stack([xy,np.full(N,z)]) for z in np.linspace(0,L,nz+1)])
    ts=[]
    for k in range(nz):
        for a,b,c in tri+k*N:
            ts.extend([[a,b,c,c+N],[a,b,b+N,c+N],[a,a+N,b+N,c+N]])
    t=np.array(ts);sgn=np.linalg.det(x[t[:,1:]]-x[t[:,:1]])
    assert np.all(abs(sgn)>1e-30)
    neg=sgn<0;t[neg,:2]=t[neg,:2][:,::-1]
    faces=t[:,[[0,1,2],[0,1,3],[0,2,3],[1,2,3]]].reshape(-1,3)
    _,first,count=np.unique(p.face_keys(faces),return_index=True,return_counts=True)
    b=faces[first[count==1]];rr=np.linalg.norm(x[b][:,:,:2],axis=2)
    wall=b[np.all(np.isclose(rr,R,atol=1e-14,rtol=0),axis=1)]
    return x,t,wall

def main():
    # Arbitrary affine tensor detects transposition, omitted transpose, factor two and projection errors.
    x=np.array([[0,0,0],[2,0,0],[0,3,0],[0,0,4]],float)*1e-6;t=np.array([[0,1,2,3]])
    G=np.array([[1,2,3],[4,5,6],[7,8,9]],float)*100
    u=x@G.T+np.array([.001,-.002,.003]);actual=p.p1_gradients(x,t,u)[0]
    ns=np.array([[0,0,1],[0,1,0],[1,0,0],[1,2,3]],float);ns/=np.linalg.norm(ns,axis=1)[:,None]
    exact=np.array([(np.eye(3)-np.outer(n,n))@(MU*(G+G.T)@n) for n in ns])
    computed=p.tangential_traction(np.repeat(actual[None],len(ns),axis=0),ns,MU)
    reverse=p.tangential_traction(np.repeat(actual[None],len(ns),axis=0),-ns,MU)
    rotation=np.array([[0,-400,300],[400,0,-200],[-300,200,0]],float)
    rotgrad=p.p1_gradients(x,t,x@rotation.T)
    rottau=p.tangential_traction(np.repeat(rotgrad,len(ns),axis=0),ns,MU)
    # Couette plane z=0: three wall nodes zero, fourth interior node nonzero.
    g=1250.;uc=np.column_stack([g*x[:,2],np.zeros(len(x)),np.zeros(len(x))])
    cg=p.p1_gradients(x,t,uc);wall=np.array([[0,1,2]]);owners=p.boundary_owners(t,wall);_,area,n=p.wall_geometry(x,t,wall,owners)
    ct=p.tangential_traction(cg[owners],n,MU);cw=np.linalg.norm(ct,axis=1)
    # Pressure reference cancels in tangent projection; direct stress test (not a new solver run).
    pressure_errors=[]
    for pr in [0.,295.57231837550194,5000.]:
        total=np.array([(np.eye(3)-np.outer(nn,nn))@((MU*(G+G.T)-pr*np.eye(3))@nn) for nn in ns])
        pressure_errors.append(float(abs(total-exact).max()))
    checks=dict(affine_gradient_max_error_s_inv=float(abs(actual-G).max()),affine_traction_max_error_Pa=float(abs(computed-exact).max()),
        normal_reversal_vector_sum_max_Pa=float(abs(computed+reverse).max()),normal_reversal_magnitude_error_Pa=float(abs(np.linalg.norm(computed,axis=1)-np.linalg.norm(reverse,axis=1)).max()),
        rigid_rotation_max_wss_Pa=float(np.linalg.norm(rottau,axis=1).max()),couette_theory_Pa=MU*g,couette_computed_Pa=float(cw[0]),couette_error_Pa=float(abs(cw[0]-MU*g)),
        pressure_reference_projection_errors_Pa=pressure_errors,meaning='Actual production gradient, normal, owner, traction functions. No PDE solver was run.')
    assert checks['affine_traction_max_error_Pa']<1e-10
    assert checks['rigid_rotation_max_wss_Pa']<1e-10
    assert checks['couette_error_Pa']<1e-10
    (A/'data/analytical_tensor_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    R=4e-6;L=40e-6;umean=.002;theory=4*MU*umean/R;Q=np.pi*R**2*umean;dp=8*MU*umean*L/R**2
    rows=[]
    for nr,nt,nz in [(2,16,4),(4,32,8),(8,64,16)]:
        x,t,wall=cylinder(nr,nt,nz,R,L);u=np.zeros_like(x);u[:,2]=2*umean*(1-(x[:,0]**2+x[:,1]**2)/R**2)
        u[np.isclose(np.linalg.norm(x[:,:2],axis=1),R,rtol=0,atol=1e-14),2]=0.
        own=p.boundary_owners(t,wall);ctr,area,n=p.wall_geometry(x,t,wall,own)
        grad=p.p1_gradients(x,t,u);tr=p.tangential_traction(grad[own],n,MU);w=np.linalg.norm(tr,axis=1)
        v,weight=p.nodal_average(wall,w,area,len(x));surf,ids=p.surface(x,wall)
        surf.cell_data['WSS_raw_Pa']=w;surf.point_data['WSS_display_Pa']=v[ids]
        surf.save(A/f'validation_cases/pipe_nr{nr}_wall.vtp')
        np.savez_compressed(A/f'validation_cases/pipe_nr{nr}_input.npz',points_m=x,tetra=t,wall_triangles=wall,velocity_m_s=u)
        sel=(ctr[:,2]>L*.2)&(ctr[:,2]<L*.8);aa=area[sel];ww=w[sel]
        row=dict(case=f'analytic_pipe_nr{nr}',validation_type='manufactured_nodal_field_not_CFD',radial_intervals=nr,circumferential_intervals=nt,axial_intervals=nz,nodes=len(x),tetra=len(t),sample_wall_facets=sel.sum(),sample_z_min_um=ctr[sel,2].min()*1e6,sample_z_max_um=ctr[sel,2].max()*1e6,R_um=R*1e6,L_um=L*1e6,Umean_m_s=umean,mu_Pa_s=MU,Q_m3_s=Q,delta_p_Pa=dp,theory_WSS_Pa=theory,radial_step_um=R/nr*1e6,computed_area_mean_Pa=np.average(ww,weights=aa),min_Pa=ww.min(),max_Pa=ww.max(),area_relative_L2_error_pct=100*np.sqrt(np.average((ww-theory)**2,weights=aa))/theory,area_mean_relative_error_pct=100*(np.average(ww,weights=aa)/theory-1))
        rows.append(row)
    with (A/'data/poiseuille_validation.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    assert all(rows[i+1]['area_relative_L2_error_pct']<rows[i]['area_relative_L2_error_pct'] for i in range(2))
    print(json.dumps(checks,indent=2));print(json.dumps([{k:(v.item() if isinstance(v,np.generic) else v) for k,v in r.items()}for r in rows],indent=2))

if __name__=='__main__':main()
