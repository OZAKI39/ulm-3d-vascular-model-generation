"""Common geometry and explicit units. Never imports either solver."""
from pathlib import Path
import math
import numpy as np
import yaml
from py_scripts.fluid_physics.common import PROJECT_ROOT, sha256_file, fingerprint

def load_config(path):
    path=Path(path).resolve(); c=yaml.safe_load(path.read_text())
    if c['schema_version']!=1 or '/' in c['campaign_id'] or c['campaign_id'] in ('.','..'):raise ValueError('BAD_CONFIG')
    if c['geometry']['axes']!={'flow':'x','gradient':'z','periodic':['x','y']}:raise ValueError('AXES_MISMATCH')
    if c['protocol']['cold_repetitions']<2:raise ValueError('TWO_COLD_REPEATS_REQUIRED')
    if c['budget']['authorized_gpu_s']!=0:raise ValueError('CONFIG_IS_NOT_AUTHORIZATION')
    if c['budget']['proposed_total_solver_s']>3600 or c['budget']['proposed_max_task_s']>1800:raise ValueError('NEW_BUDGET_REVIEW_REQUIRED')
    c['_config_path']=str(path);c['_config_sha256']=sha256_file(path);return c

def paths(c):
    name=c['campaign_id'];return (PROJECT_ROOT/'data/single_rbc_benchmark'/name,PROJECT_ROOT/'runs/single_rbc_benchmark'/name,PROJECT_ROOT/'test_code/outputs/single_rbc_benchmark'/name)

def units():
    # A complete bookkeeping embedding required by HemoCell's SI-based API.
    # These constants define numerical units, not a calibrated physical fluid.
    L0=1e-6; E0=4.100531391e-21; M0=1000*L0**3/8; T0=math.sqrt(M0*L0*L0/E0)
    return dict(L0=L0,E0=E0,M0=M0,T0=T0,F0=E0/L0,stress0=E0/L0**3,
                rho0=M0/L0**3,nu0=L0*L0/T0,mu0=M0/(L0*T0),Gs0=E0/L0**2,
                interpretation='API bookkeeping embedding; reported time is common t*, not physiological seconds')

def lattice(c,nu,strict=False):
    d=c['hemocell']['dx']; dt=(c['hemocell']['tau_candidate']-.5)*d*d/(3*nu)
    # Make all strain sampling/endpoints exact integer multiples of the same dt.
    interval=c['protocol']['sample_strain']/c['protocol']['shear_rate']
    dt=interval/math.ceil(interval/dt)
    if strict:dt/=2
    return dict(dx=d,dt=dt,tau=.5+3*nu*dt/d**2,nx=round(c['geometry']['periodic_length']/d),nz=round(c['geometry']['gap']/d)+1)

def groups(rho,nu,shear,a,Gs,kb,mass,volume,area,eta=0):
    if min(rho,nu,a,Gs,volume,area)<=0:raise ValueError('NONPOSITIVE_MATERIAL')
    mu=rho*nu
    return dict(Re=shear*a*a/nu,Ca=mu*shear*a/Gs,lambda_viscosity=1.,B=kb/(Gs*a*a),
                surface_Boussinesq=eta/(mu*a),membrane_to_displaced_fluid_mass=mass/(rho*volume),
                thermal_to_shear_energy=1/(Gs*a*a),a_definition='sqrt(A0/(4*pi))')

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

def unwrap(v,faces,lengths,periodic=(True,True,False)):
    v=np.asarray(v,float);out=np.full_like(v,np.nan);seen=np.zeros(len(v),bool);out[0]=v[0];seen[0]=True
    adjacency=[set() for _ in v]
    for face in faces:
        for i,j in zip(face,np.roll(face,-1)):adjacency[i].add(j);adjacency[j].add(i)
    queue=[0];lengths=np.array(lengths,float);mask=np.array(periodic)
    while queue:
        i=queue.pop()
        for j in adjacency[i]:
            delta=v[j]-v[i];delta[mask]-=lengths[mask]*np.rint(delta[mask]/lengths[mask]);point=out[i]+delta
            if seen[j]:
                if not np.allclose(point,out[j],atol=1e-5):raise ValueError('INCONSISTENT_UNWRAPPING_OR_CELL_SPANS_BOX')
            else:out[j]=point;seen[j]=True;queue.append(j)
    if not seen.all():raise ValueError('DISCONNECTED_CELL')
    return out

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

def sample_steps(c,dt,mode='main'):
    p=c['protocol'];prep=round(p['relaxation_time']/dt);steps=round(p['end_strain']/p['shear_rate']/dt);every=round(p['sample_strain']/p['shear_rate']/dt)
    if not math.isclose(steps*dt*p['shear_rate'],p['end_strain'],abs_tol=1e-8):raise ValueError('STRAIN_ENDPOINT_NOT_RESOLVED')
    return dict(prep_steps=prep,steps=steps,sample_steps=every)

def speedup(cpu_s,gpu_s,*,workflow,model,screen,repeats):
    if workflow!='PASS' or model!='PASS' or screen!='PASS' or repeats<2:return None
    if min(cpu_s,gpu_s)<=0:raise ValueError('INVALID_TIMING')
    return cpu_s/gpu_s

def xml_case(c,nu=1.5,strict=False,cell=True,probe=False,short_steps=None):
    u=units();l=lattice(c,nu,strict);s=sample_steps(c,l['dt'])
    if short_steps is not None:s.update(prep_steps=0,steps=short_steps,sample_steps=max(1,short_steps//4))
    def block(name,d):return '<'+name+'>'+''.join(f'<{k}>{v}</{k}>' for k,v in d.items())+'</'+name+'>'
    return '<?xml version="1.0"?>\n<hemocell>'+block('parameters',dict(outputDirectory='native_output',warmup=0))+block('domain',dict(rhoP=8*u['rho0'],nuP=nu*u['nu0'],dx=l['dx']*u['L0'],dt=l['dt']*u['T0'],kBT=u['E0'],particleEnvelope=24,shearrate=c['protocol']['shear_rate']/u['T0']))+block('benchmark',dict(dxStar=l['dx'],dtStar=l['dt'],length=c['geometry']['periodic_length'],gap=c['geometry']['gap'],shear=c['protocol']['shear_rate'],cell=int(cell),probe=int(probe)))+block('sim',dict(prepSteps=s['prep_steps'],steps=s['steps'],sampleEvery=s['sample_steps'],tmax=s['steps']+s['prep_steps'],tmeas=s['sample_steps'],tcheckpoint=999999999))+'</hemocell>'

def rbc_files(c,d):
    d=Path(d);p={**c['hemocell']['material'],'radius':c['geometry']['radius']*1e-6,'minNumTriangles':c['geometry']['min_triangles'],'Volume':90,'enableInteriorViscosity':0,'viscosityRatio':1}
    (d/'RBC.xml').write_text('<?xml version="1.0"?>\n<hemocell><MaterialModel>'+''.join(f'<{k}>{v}</{k}>' for k,v in p.items())+'</MaterialModel></hemocell>')
    (d/'RBC.pos').write_text('1\n'+' '.join(map(str,c['geometry']['center']+c['geometry']['hemocell_angles_degrees']))+'\n')
