"""P7 stable-ID dynamic insertion/deletion extension; P6 bridge remains unchanged."""
import numpy as np
from .lammps_bridge import LammpsParticleBridge
from .lammps_state import PROPERTY_COMMAND


class DynamicParticleBridge(LammpsParticleBridge):
    def __init__(self,particles,policy,*,box_half_width_m=.001):
        self.used_ids=set()
        particles=list(particles)
        if particles:
            super().__init__(particles,policy,box_half_width_m=box_half_width_m)
            self.used_ids.update(p.particle_id for p in particles)
            return
        self._start(policy,box_half_width_m)
        try:
            w=box_half_width_m
            for c in ['units si','atom_style atomic','atom_modify map array sort 0 0.0','boundary f f f',
                f'region bridge_box block {-w:.17g} {w:.17g} {-w:.17g} {w:.17g} {-w:.17g} {w:.17g} units box',
                'create_box 2 bridge_box','mass * 1',PROPERTY_COMMAND]: self.command(c)
            self._query_setup(); self.rebuild()
        except BaseException: self.close(); raise
    def _indices(self):
        if self._lmp.extract_setting('nlocal')==0: return {}
        return super()._indices()
    def read(self):
        return super().read() if self._indices() else []
    def audit_unused_force_storage(self,phase):
        if self._lmp.extract_setting('nlocal')==0:
            record=dict(phase=phase,max_lammps_force=0.,max_lammps_torque=0.,used_as_physics_input=False,lammps_step=int(self._lmp.extract_global('ntimestep')))
            self.force_audits.append(record); return record
        return super().audit_unused_force_storage(phase)
    def insert(self,particle):
        if particle.particle_id in self.used_ids: raise ValueError('STABLE_ID_NEVER_REUSED')
        if np.any(abs(particle.position)>=self.box_half_width_m): raise ValueError('Particle outside fixed storage box')
        old=self.read(); self.audit_unused_force_storage('BEFORE_INSERT')
        n=self._lmp.create_atoms(1,[particle.particle_id],[particle.type_code],particle.position.tolist())
        if n!=1: raise ValueError('Insertion failed')
        i=self._indices()[particle.particle_id]; self._lmp.numpy.extract_atom('f')[i]=0.
        self.write(old+[particle]); self.used_ids.add(particle.particle_id); self.rebuild()
    def remove(self,particle_id):
        if particle_id not in self._indices(): raise ValueError('Unknown active stable ID')
        self.audit_unused_force_storage('BEFORE_DELETE')
        # Only these literal storage commands bypass P6's frozen command whitelist.
        for c in [f'group p7_delete id {int(particle_id)}','delete_atoms group p7_delete compress no','group p7_delete delete']:
            self.commands.append(c); self._lmp.command(c)
        self.rebuild()
    @classmethod
    def from_binary_restart(cls,path,policy,box_half_width_m):
        obj=super().from_binary_restart(path,policy,box_half_width_m)
        obj.used_ids={p.particle_id for p in obj.read()}
        return obj
