"""Thin deterministic adapter. No sampling, resizing or distribution fitting.
Diameter input is the frozen sampler CSV in micrometres, converted exactly once.
Sphere density here is a technical engine test value, not SonoVue material data.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np

DUMP_FIELDS='id type diameter radius mass x y z xu yu zu vx vy vz omegax omegay omegaz fx fy fz proc ix iy iz'

def load_population(path):
    p=Path(path)
    a=np.atleast_1d(np.genfromtxt(p,delimiter=',',names=True,dtype=None,encoding='utf-8'))
    if not {'bubble_id','diameter_um'} <= set(a.dtype.names or ()):
        raise ValueError('Missing required population columns')
    ids=np.asarray(a['bubble_id'],dtype=np.int64)
    if not np.array_equal(ids,a['bubble_id']) or not np.array_equal(ids,np.arange(len(ids))):
        raise ValueError('Frozen sampler population must have ordered zero-based integer IDs')
    d=np.asarray(a['diameter_um'],dtype=np.float64)
    if not np.all(np.isfinite(d)&(d>0)): raise ValueError('Invalid diameter')
    return ids+1, d*1e-6

def initial_state(d, spec):
    n=len(d); box=np.asarray(spec['box_m'],dtype=np.float64)
    if spec['placement']=='grid':
        shape=np.asarray(spec['grid_shape'],dtype=int)
        if np.prod(shape)!=n: raise ValueError('Grid must contain exactly N particles')
        ijk=np.array(np.unravel_index(np.arange(n),tuple(shape))).T
        xyz=(ijk+0.5)*spec['spacing_m']
        if spec['spacing_m']<=np.max(d)+spec['grid_clear_margin_m']:
            raise ValueError('Grid spacing must exceed max diameter plus margin')
    elif spec['placement']=='single':
        xyz=np.array([box*0.5])
    elif spec['placement']=='two':
        if n!=2: raise ValueError('Expected two spheres')
        separation=0.5*(d[0]+d[1])+spec['surface_gap_m']
        xyz=np.array([box*0.5,box*0.5]);xyz[:,0]+=np.array([-0.5,0.5])*separation
    else: raise ValueError('Unsupported deterministic placement')
    velocity=np.broadcast_to(np.array(spec['velocity_m_s'],dtype=float),(n,3)).copy()
    omega=np.broadcast_to(np.array(spec['omega_rad_s'],dtype=float),(n,3)).copy()
    if not np.all(np.isfinite(xyz)) or not np.all((xyz>0)&(xyz<box)): raise ValueError('Invalid initial coordinates')
    return xyz,velocity,omega

def write_case(population, spec, destination):
    out=Path(destination); out.mkdir(parents=True,exist_ok=True)
    ids,d=load_population(population)
    if len(ids)!=spec['N']: raise ValueError('N mismatch')
    xyz,v,om=initial_state(d,spec)
    density=spec['technical_density_kg_m3']
    with (out/'particles.data').open('x') as f:
        f.write('Frozen SonoVue sampler input; technical rigid-sphere engine case\n\n')
        f.write(f'{len(ids)} atoms\n1 atom types\n\n')
        for size,axis in zip(spec['box_m'],'xyz'): f.write(f'0 {size:.17g} {axis}lo {axis}hi\n')
        f.write('\nAtoms # sphere\n\n')
        atoms=np.column_stack((ids,np.ones(len(ids),dtype=int),d,np.full(len(ids),density),xyz))
        np.savetxt(f,atoms,fmt=['%d','%d']+['%.17g']*5)
        f.write('\nVelocities\n\n')
        np.savetxt(f,np.column_stack((ids,v,om)),fmt=['%d']+['%.17g']*6)
    lines=['units si','dimension 3','boundary p p p','atom_style sphere','newton off',
           'processors '+str(spec['mpi_ranks'])+' 1 1','read_data particles.data','comm_modify vel yes',
           'neighbor 1e-6 bin','neigh_modify every 1 delay 0 check no']
    if spec['contact']:
        lines += ['pair_style gran/hooke/history 0.001 0 0 0 0 0','pair_coeff * *']
    else: lines += ['pair_style zero 6e-6','pair_coeff * *']
    lines += ['fix integrate all nve/sphere',f'timestep {spec["dt_s"]:.17g}',
              'thermo '+str(max(1,spec['steps']//10)),
              'thermo_style custom step atoms time ke','thermo_modify lost error norm no format float %.17g','run 0',
              'write_dump all custom initial.dump '+DUMP_FIELDS+' modify sort id format float %.17g']
    if spec['steps']:
        lines += [f'run {spec["steps"]}', 'write_dump all custom final.dump '+DUMP_FIELDS+' modify sort id format float %.17g']
    else: lines += ['write_dump all custom final.dump '+DUMP_FIELDS+' modify sort id format float %.17g']
    lines += ['print "ENGINE_CASE_COMPLETED"']
    (out/'in.lammps').write_text('\n'.join(lines)+'\n')
    (out/'CASE_SPEC.json').write_text(json.dumps(spec,indent=2)+'\n')
    return {'N':len(ids),'diameter_m_min':float(d.min()),'diameter_m_max':float(d.max()),
            'population_sha256':hashlib.sha256(Path(population).read_bytes()).hexdigest(),
            'id_mapping':'LAMMPS id = sampler bubble_id + 1','conversion':'diameter_um * 1e-6',
            'diameter_changed':'NO','resampling':'NO'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--population',required=True);ap.add_argument('--spec',required=True);ap.add_argument('--output',required=True)
    a=ap.parse_args();print(json.dumps(write_case(a.population,json.loads(Path(a.spec).read_text()),a.output),indent=2))
