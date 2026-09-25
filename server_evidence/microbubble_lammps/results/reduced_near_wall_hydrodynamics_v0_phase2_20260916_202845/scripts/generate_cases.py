"""Only synthetic analytic planes. Initial coordinate rounding is recorded, never clamped in engine."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
PATHS=json.loads((ROOT/'provenance/STAGE_PATHS.json').read_text())
a=2.**-20
formal=[.001,.002,.005,.01,.02,.05,.08,.1,.2]
cases=[]

def rotation(axis,angle):
    axis=np.asarray(axis,dtype=float);axis/=np.linalg.norm(axis)
    k=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    return np.eye(3)+np.sin(angle)*k+(1-np.cos(angle))*(k@k)

def add(name,e,gamma=100.,Q=None,nr=1,steps=64,matrix=False,expected='COMPLETE',group='formal'):
    Q=np.eye(3) if Q is None else Q
    local=np.array([-2*a,0,(1+e)*a]);nominal=local.copy();fixture_steps=0
    if e in [.001,.2]:
        while (local[2]-a)/a<.001 or (local[2]-a)/a>.2:
            local[2]=np.nextafter(local[2],np.inf if e==.001 else -np.inf);fixture_steps+=1
        assert fixture_steps<=1
    x=Q@local
    eps=float((np.sum(x.astype(np.longdouble)*Q[:,2].astype(np.longdouble))-np.longdouble(a))/a)
    if expected=='COMPLETE':assert .001<=eps<=.2,(name,eps)
    case={'name':name,'epsilon_nominal':e,'initial_epsilon_from_coordinates':eps,'radius':a,'initial_position':x.tolist(),'initial_position_wall_frame':local.tolist(),'nominal_initial_position_wall_frame':nominal.tolist(),'input_center_rounding_ULPs':fixture_steps,'gamma_dot':gamma,'Q':Q.tolist(),'origin':[0.,0.,0.],'mpi_ranks':nr,'steps':steps,'dt':1e-4,'max_time':steps*1e-4,'expected':expected,'group':group,'dump_matrix':matrix}
    target=ROOT/'cases'/name;target.mkdir(exist_ok=True)
    values={'geometry':'ANALYTIC_PLANE','shear_source':'ANALYTIC','resistance_table':PATHS['table_remote'],'data':str(Path(PATHS['remote_stage'])/'cases'/name/'particle.data'),'dt_max':case['dt'],'max_time':case['max_time'],'steps_cap':steps+10,'stride':1,'gamma_dot':gamma,'dump_matrix':'YES' if matrix else 'NO'}
    for label,v in [('origin_',case['origin']),('t1_',Q[:,0]),('t2_',Q[:,1]),('n_',Q[:,2])]:
        for axis,val in zip('xyz',v):values[label+axis]=float(val)
    (target/'case.cfg').write_text(''.join(f'{k} {v:.17g}\n' if isinstance(v,float) else f'{k} {v}\n' for k,v in values.items()))
    # Fixed nonperiodic computational box; MPI 4 decomposes its x direction.
    text='Phase 2 single sphere synthetic data\n\n1 atoms\n1 atom types\n\n'
    for axis in 'xyz':text+=f'{-20*a:.17g} {20*a:.17g} {axis}lo {axis}hi\n'
    text+='\nAtoms # sphere\n\n'+f'1 1 {2*a:.17g} 1000 '+' '.join(f'{v:.17g}' for v in x)+'\n'
    (target/'particle.data').write_text(text);(target/'initial.json').write_text(json.dumps(case,indent=2)+'\n');cases.append(case)

for i,e in enumerate(formal):
    for suffix,gamma in [('zero',0.),('reverse',-100.),('half',50.),('standard',100.),('double',200.)]:
        add(f'gap{i:02d}_{suffix}',e,gamma,steps=1000 if suffix=='standard' else 64)
    add(f'gap{i:02d}_mpi4',e,nr=4,steps=1000,group='mpi')
rotations=[rotation([1,0,0],.7),rotation([0,1,0],.5),rotation([1,2,3],.8)]
for j,Q in enumerate(rotations):
    for i in [2,4,6]:add(f'rot{j}_gap{i:02d}',formal[i],Q=Q,steps=256,group='rotation')
add('rot2_gap04_mpi4',.02,Q=rotations[2],nr=4,steps=256,group='mpi_rotation')
for i,e in enumerate([.0095,.01,.0105,.095,.1,.105,.01-1e-12,.01+1e-12,.1-1e-11,.1+1e-11]):
    add(f'transition{i:02d}',e,steps=2,matrix=True,group='transition')
for name,e,status in [('below_domain',.0009,'OUTSIDE_CERTIFIED_WALL_DOMAIN'),('above_domain',.2001,'OUTSIDE_CERTIFIED_WALL_DOMAIN'),('zero_gap',0.,'INVALID_GAP'),('negative_gap',-.01,'INVALID_GAP')]:
    add(name,e,steps=1,expected=status,group='rejection')
# Explicit unsupported-source rejection through the real executable's configuration path.
for source in ['GRADIENT_PROXY','PALABOS_VALIDATED']:
    name='reject_'+source.lower();add(name,.02,steps=1,expected='UNSUPPORTED_SOURCE',group='rejection')
    p=ROOT/'cases'/name/'case.cfg';p.write_text(p.read_text().replace('shear_source ANALYTIC','shear_source '+source))
(ROOT/'contracts/CASE_INDEX.json').write_text(json.dumps(cases,indent=2)+'\n')
print(json.dumps({'cases':len(cases),'formal_gaps':formal,'endpoint_policy':'At most one input center ULP inward before writing case; actual geometry always independently checked. Engine uses strict rejection.'}))
