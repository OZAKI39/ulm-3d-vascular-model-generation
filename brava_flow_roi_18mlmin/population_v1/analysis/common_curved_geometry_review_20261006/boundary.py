from common import *
from validate import EDGES,RN,basis,geometry
import numpy as np
import subprocess

def main():
    assert read(HERE/'COMMON_GEOMETRY_HARNESS_RESULTS.json')['status']=='PASS'
    assert read(HERE/'P1_REGRESSION_RESULTS.json')['status']=='PASS'
    old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json');rt=HERE/'runtime'
    cmd=[a.replace(old['source_checkout'],str(RETAINED)) for a in old['compile_command']]
    cmd[cmd.index('-c')+1]=str(HERE/'native_boundary.cpp');cmd[cmd.index('-o')+1]=str(rt/'boundary.o');cmd+=['-I'+str(HERE)]
    link=list(old['link_command'])
    for i,a in enumerate(link):
        if a.endswith('/main.cpp.o'):link[i]=str(rt/'boundary.o')
        if a.endswith('/fluid.cpp.o'):link[i]=str(rt/'patched.o')
    link[link.index('-o')+1]=str(rt/'native_boundary')
    build=dict(compile_command=cmd,link_command=link,cwd=old['link_cwd'],source_sha256=sha(HERE/'native_boundary.cpp'),
        boundary_source_unmodified=True,protocol_sha256=sha(HERE/'BOUNDARY_HARNESS_PROTOCOL.json'))
    put('BOUNDARY_BUILD_IDENTITY.json',build)
    subprocess.run(cmd,cwd=old['link_cwd'],check=True);subprocess.run(link,cwd=old['link_cwd'],check=True)
    build['binary_sha256']=sha(rt/'native_boundary');put('BOUNDARY_BUILD_IDENTITY.json',build)
    data=np.load(PREVIOUS/'P2_CURVED_GEOMETRY_SAMPLES.npz')
    frozen=read(BASE/'production_workflow/gmsh_reference/GMSH_REFERENCE_MESH_V2_PLAN.json')
    R,L=frozen['radius_m'],frozen['length_m'];mu=.00345312;U=.1
    fixtures=[];tests=[];arrays={}
    def check(name,value,limit=1e-11):
        tests.append(dict(name=name,error=float(value),limit=limit,status='PASS' if value<=limit else 'FAIL'))
    for cell,kind in [('central_wall','wall'),('inlet_rim','wall'),('inlet_rim','inlet_cap'),('inlet_rim','outlet_rigid_control')]:
      for variant in ['affine','curved']:
        tag=f'{cell}_{kind}_{variant}';X=data[f'{cell}_{variant}_nodes'].copy()
        if kind=='outlet_rigid_control':X[:,1]*=-1;X[:,2]=L-X[:,2]
        H=float(np.linalg.norm(X[:4,None]-X[None,:4],axis=2).max())
        vertices=[0,2,3] if kind=='wall' else [0,1,3]
        ids=vertices+[4+next(i for i,e in enumerate(EDGES) if set(e)==set(pair)) for pair in
            [(vertices[0],vertices[1]),(vertices[1],vertices[2]),(vertices[2],vertices[0])]]
        XF=X[ids]
        if kind=='wall':traction=np.column_stack([.2+XF[:,0]/R,-.1+XF[:,1]/R,.4+XF[:,2]/L])
        else:
            sign=-1 if kind=='inlet_cap' else 1
            traction=sign*np.column_stack([-4*mu*U/R**2*XF[:,0],-4*mu*U/R**2*XF[:,1],np.full(6,-.4)])
        inp=rt/(tag+'_boundary_input.txt');pre=rt/tag
        with inp.open('w') as f:
            np.savetxt(f,X,fmt='%.17g');f.write(' '.join(map(str,ids))+'\n');np.savetxt(f,traction,fmt='%.17g')
        subprocess.run([str(rt/'native_boundary'),str(inp),str(pre)],env=native_env(),check=True)
        raw=np.loadtxt(str(pre)+'_q.txt');w=raw[:,1];qxi=raw[:,2:4];N=raw[:,4:10];D=raw[:,10:22].reshape(-1,6,2)
        xp=raw[:,22:25];nv=raw[:,25:28];known=raw[:,28:31];vxi=raw[:,31:34];pxi=raw[:,34:37];valid=raw[:,37:39].astype(bool)
        lt=np.column_stack([qxi,1-qxi.sum(1)]);dt=np.vstack([np.eye(2),-np.ones(2)]);edges=np.array([[0,1],[1,2],[2,0]])
        Nr=np.column_stack([lt*(2*lt-1),4*lt[:,edges[:,0]]*lt[:,edges[:,1]]])
        Dr=np.concatenate([(4*lt-1)[:,:,None]*dt,4*(lt[:,edges[:,0],None]*dt[edges[:,1]]+lt[:,edges[:,1],None]*dt[edges[:,0]])],axis=1)
        tangent=np.einsum('ai,naj->nij',XF-XF[2],Dr)
        normal=np.cross(tangent[:,:,0],tangent[:,:,1]);inside=[i for i in range(4) if i not in vertices][0]
        normal[np.einsum('ni,i->n',normal,XF[0]-X[inside])<0]*=-1
        jac=np.linalg.norm(normal,axis=1);x_ref=Nr@(XF-XF[2])+XF[2]
        full_ref=np.array([geometry(q,X)['x'] for q in known])
        check(tag+'_TRI6_basis',float(abs(N-Nr).max()))
        check(tag+'_TRI6_derivatives',float(abs(D-Dr).max()))
        check(tag+'_surface_volume_position',float(abs(full_ref-x_ref).max())/H,1e-12)
        check(tag+'_native_surface_position',float(abs(xp-x_ref).max())/H,1e-12)
        check(tag+'_native_surface_normal',float(abs(nv-normal).max())/float(abs(normal).max()))
        reference=-np.einsum('n,na,ni->ia',w*jac,Nr,Nr@traction)
        lines=Path(str(pre)+'_traction_RK.txt').read_text().splitlines();actual=np.fromstring(lines[0],sep=' ').reshape(4,6)
        check(tag+'_native_traction_residual',float(abs(actual[:3]-reference).max())/(H*H*max(float(abs(traction).max()),1.)))
        dirs=np.loadtxt(str(pre)+'_dir.txt').reshape(10,4,2);expected=dirs[:,:,0].copy();expected[ids,:3]=0.
        check(tag+'_six_node_strong_Dirichlet_and_preserved_other_DOF',float(abs(dirs[:,:,1]-expected).max()),0.)
        pressure_common=np.column_stack([known,1-known.sum(1)]);pressure_native=np.column_stack([pxi,1-pxi.sum(1)])
        pressure_error=float(abs(pressure_native-pressure_common).max())
        velocity_ref_error=float(abs(vxi-known).max())
        native_status=read(Path(str(pre)+'_status.json'))
        defect=pressure_error>1e-11 or not valid[:,1].all() or not native_status['native_weak_Dirichlet_completed']
        physical_reconstruction=np.array([geometry(q,X)['x'] for q in vxi])
        if kind!='wall':
            capz=0. if kind=='inlet_cap' else L
            check(tag+'_cap_plane',float(abs(xp[:,2]-capz).max())/H,1e-12)
            sign=-1 if kind=='inlet_cap' else 1
            exact_t=sign*np.column_stack([-4*mu*U/R**2*x_ref[:,0],-4*mu*U/R**2*x_ref[:,1],np.full(len(w),-.4)])
            check(tag+'_full_analytic_Cauchy_at_q',float(abs(Nr@traction-exact_t).max())/max(float(abs(exact_t).max()),1.))
            exact_R=-np.einsum('n,na,ni->ia',w*jac,Nr,exact_t)
            check(tag+'_full_Cauchy_native_integral',float(abs(actual[:3]-exact_R).max())/(H*H*max(float(abs(traction).max()),1.)))
        errpath=Path(str(pre)+'_weak_error.txt')
        fixtures.append(dict(case=tag,source_cell=143654 if cell=='central_wall' else 35265,kind=kind,variant=variant,
            rigid_control=kind=='outlet_rigid_control',native=native_status,face_local_nodes=ids,
            pressure_reference_coordinate_difference_max=pressure_error,
            velocity_reference_coordinate_difference_max=velocity_ref_error,
            native_velocity_inverse_position_error_over_H=float(abs(physical_reconstruction-xp).max())/H,
            pressure_inverse_admissible_points=int(valid[:,1].sum()),total_quadrature_points=len(w),
            same_geometry_defect=bool(defect),native_weak_exception=errpath.read_text() if errpath.exists() else None,
            general_pressure_field_basis_difference_max=float(abs((pressure_native-pressure_common)@np.array([.2,.8,-.3,.1])).max())))
        for k,v in dict(nodes=X,face_ids=np.array(ids),native_rows=raw,reference_normal=normal,reference_x=x_ref,
            common_pressure_basis=pressure_common,native_pressure_basis=pressure_native,traction_nodes=traction,
            reference_traction_R=reference,native_traction_R=actual,strong_Dirichlet=dirs).items():arrays[tag+'_'+k]=v
    blocked=any(f['same_geometry_defect'] for f in fixtures)
    np.savez_compressed(HERE/'BOUNDARY_MIXED_GEOMETRY_SAMPLES.npz',**arrays)
    put('BOUNDARY_MIXED_GEOMETRY_AUDIT.json',dict(status='MIXED_CURVED_BOUNDARY_ASSEMBLY_BLOCKED' if blocked else 'PASS',
        fixtures=fixtures,tests=tests,engineering_failed=sum(t['status']=='FAIL' for t in tests),
        root_cause='set_bc_dir_wl independently calls get_nnx with the four-corner pressure geometry after computing the curved velocity reference coordinate',
        full_Cauchy_scope='Native prescribed TRI6 traction load and planar cap full analytical Poiseuille Cauchy integration tested; weak-Dirichlet pressure path is separate and remains blocked',
        no_boundary_patch=True,CFD_steps=0,linear_solves=0,automatic_followup=False))
    print('Boundary audit:', 'BLOCKED' if blocked else 'PASS', 'engineering failures',sum(t['status']=='FAIL' for t in tests))
    if any(t['status']=='FAIL' for t in tests):raise SystemExit('Boundary fixture engineering check failed')

if __name__=='__main__':main()
