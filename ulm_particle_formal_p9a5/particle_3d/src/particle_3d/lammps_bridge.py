"""Serial storage/neighbor bridge only. No LAMMPS particle physics or time stepping."""
from pathlib import Path
import ctypes,shlex,sys
import numpy as np
from .lammps_state import BridgeParticle,PROPERTIES,PROPERTY_COMMAND
from .lammps_neighbors import canonical_pairs


def new_lammps():
    # Official wheel needs MPICH ABI. Load the project-venv library only, with
    # no system links, environment changes, or writes to the Frozen FEM venv.
    library=Path(sys.prefix)/'lib/libmpi.so.12'
    if library.exists():ctypes.CDLL(str(library),mode=ctypes.RTLD_GLOBAL)
    from lammps import lammps
    obj=lammps(cmdargs=['-log','none','-screen','none','-nocite'])
    if obj.extract_setting('world_size')!=1:
        obj.close();raise ValueError('P6_V0_SINGLE_MPI_RANK_REQUIRED')
    return obj


def validate_command(command):
    if any(c in command for c in ['\n','\r',';','$','&']):raise ValueError('Only one literal bridge command permitted')
    p=shlex.split(command)
    allowed={'units','atom_style','atom_modify','boundary','region','create_box','mass','fix','pair_style','pair_coeff','neighbor','neigh_modify','run','read_restart','write_restart','write_dump','thermo'}
    if not p or p[0] not in allowed:raise ValueError('FORBIDDEN_LAMMPS_DYNAMICS_COMMAND')
    if p[0]=='fix' and command!=PROPERTY_COMMAND:raise ValueError('Only frozen property/atom fix is allowed')
    if p[0]=='run' and p not in [['run','0'],['run','0','post','no']]:raise ValueError('LAMMPS_NONZERO_RUN_FORBIDDEN')
    if p[0]=='pair_style' and (len(p)<3 or p[1]!='zero'):raise ValueError('Only pair_style zero permitted')
    if p[0]=='mass' and p!=['mass','*','1']:raise ValueError('Mass is a placeholder only')
    if p[0]=='atom_style' and p!=['atom_style','atomic']:raise ValueError('Physical geometry is not a LAMMPS atom style')


class LammpsParticleBridge:
    def __init__(self,particles,policy,*,box_half_width_m=.001):
        particles=list(particles)
        if not particles or len({p.particle_id for p in particles})!=len(particles):raise ValueError('Nonempty unique stable particle IDs required')
        self._start(policy,box_half_width_m)
        try:
            for c in ['units si','atom_style atomic','atom_modify map array sort 0 0.0','boundary f f f',
                      f'region bridge_box block {-box_half_width_m:.17g} {box_half_width_m:.17g} {-box_half_width_m:.17g} {box_half_width_m:.17g} {-box_half_width_m:.17g} {box_half_width_m:.17g} units box',
                      'create_box 2 bridge_box','mass * 1',PROPERTY_COMMAND]:self.command(c)
            if any(np.any(abs(p.position)>=box_half_width_m) for p in particles):raise ValueError('Particle outside fixed storage box')
            created=self._lmp.create_atoms(len(particles),[p.particle_id for p in particles],[p.type_code for p in particles],np.concatenate([p.position for p in particles]).tolist())
            if created!=len(particles):raise ValueError('LAMMPS atom insertion count mismatch')
            self._lmp.numpy.extract_atom('f')[:len(particles)]=0. # initialize unused storage, never clear after a force audit
            self.write(particles)
            self._query_setup();self.rebuild()
        except BaseException:
            self.close();raise

    def _start(self,policy,box_half_width_m):
        self.policy=policy;self.box_half_width_m=box_half_width_m
        self.commands=[];self.force_audits=[];self.closed=False;self._lmp=new_lammps()

    def command(self,text):
        validate_command(text);self.commands.append(text);self._lmp.command(text)

    def _query_setup(self):
        for c in [f'pair_style zero {self.policy.center_cutoff_m:.17g}','pair_coeff * *',f'neighbor {self.policy.skin_m:.17g} nsq','neigh_modify delay 0 every 1 check no','thermo 0']:
            self.command(c)

    @property
    def version(self):return self._lmp.version()

    def close(self):
        if not self.closed:self._lmp.close();self.closed=True

    def __enter__(self):return self
    def __exit__(self,*args):self.close()

    def _indices(self):
        count=self._lmp.extract_setting('nlocal');ids=self._lmp.numpy.extract_atom('id')[:count]
        if count!=int(self._lmp.get_natoms()):raise ValueError('Serial atom count mismatch')
        return {int(tag):k for k,tag in enumerate(ids)}

    def write(self,particles):
        records={p.particle_id:p for p in particles};indices=self._indices()
        if set(records)!=set(indices):raise ValueError('STATE_WRITE_MUST_PRESERVE_STABLE_IDS')
        for p in records.values():
            if np.any(abs(p.position)>=self.box_half_width_m):raise ValueError('Particle outside fixed storage box')
        x=self._lmp.numpy.extract_atom('x');types=self._lmp.numpy.extract_atom('type')
        arrays={name:self._lmp.numpy.extract_atom(name) for name,_,_ in PROPERTIES}
        for tag,index in indices.items():
            p=records[tag];x[index]=p.position;types[index]=p.type_code
            for name,key,_ in PROPERTIES:arrays[name][index]=getattr(p,key)

    def read(self):
        indices=self._indices();arrays={name:self._lmp.numpy.extract_atom(name) for name,_,_ in PROPERTIES}
        x=self._lmp.numpy.extract_atom('x');types=self._lmp.numpy.extract_atom('type');records=[]
        for tag,index in sorted(indices.items()):
            fields={key:(arrays[name][index].copy() if width>1 else arrays[name][index].item()) for name,key,width in PROPERTIES}
            if fields['type_code']!=int(types[index]):raise ValueError('LAMMPS type and custom particle type disagree')
            records.append(BridgeParticle(tag,position=x[index].copy(),**fields))
        return records

    def audit_unused_force_storage(self,phase):
        count=self._lmp.extract_setting('nlocal')
        force=self._lmp.numpy.extract_atom('f')[:count]
        torque=self._lmp.numpy.extract_atom('torque')
        max_force=float(np.max(np.abs(force),initial=0.));max_torque=0. if torque is None else float(np.max(np.abs(torque[:count]),initial=0.))
        record=dict(phase=phase,max_lammps_force=max_force,max_lammps_torque=max_torque,
            torque_storage='UNDEFINED_FOR_ATOM_STYLE_ATOMIC' if torque is None else 'DEFINED_BUT_UNUSED',
            used_as_physics_input=False,velocity_source='PARTICLE_3D_RESISTANCE_SOLVER_OR_UPSTREAM_P4_VALIDATION',lammps_step=int(self._lmp.extract_global('ntimestep')))
        if max_force!=0 or max_torque!=0:raise ValueError('NONZERO_UNUSED_LAMMPS_FORCE_OR_TORQUE')
        self.force_audits.append(record);return record

    def rebuild(self):
        self.audit_unused_force_storage('BEFORE_NEIGHBOR_BUILD')
        self.command('run 0 post no')
        self.audit_unused_force_storage('AFTER_NEIGHBOR_BUILD')
        index=self._lmp.find_pair_neighlist('zero')
        if index<0:raise ValueError('LAMMPS_NEIGHBOR_LIST_UNAVAILABLE')
        tags=self._lmp.numpy.extract_atom('id');nl=self._lmp.numpy.get_neighlist(index)
        pairs=[]
        for local,neighbors in nl:
            for neighbor in neighbors:
                other=int(neighbor)&0x3fffffff # strip LAMMPS special-neighbor flag bits
                pairs.append((int(tags[local]),int(tags[other])))
        self.raw_pairs=canonical_pairs(pairs)
        return self.raw_pairs

    def dump(self,path):
        columns=['id','type','x','y','z']
        for name,_,width in PROPERTIES:columns.extend([name] if width==1 else [f'{name}[{i}]' for i in range(1,width+1)])
        self.command('write_dump all custom '+shlex.quote(str(Path(path).resolve()))+' '+' '.join(columns)+' modify sort id')

    @classmethod
    def from_binary_restart(cls,path,policy,box_half_width_m):
        obj=cls.__new__(cls);obj._start(policy,box_half_width_m)
        try:
            obj.command('atom_modify map array sort 0 0.0')
            obj.command('read_restart '+shlex.quote(str(Path(path).resolve())))
            # Exact same ID, property order, names and types, AFTER read_restart.
            obj.command(PROPERTY_COMMAND)
            obj._query_setup()
            # Force arrays are transient, not physics/restart state.
            obj._lmp.numpy.extract_atom('f')[:obj._lmp.extract_setting('nlocal')]=0.
            obj.rebuild();obj.read()
            return obj
        except BaseException:obj.close();raise
