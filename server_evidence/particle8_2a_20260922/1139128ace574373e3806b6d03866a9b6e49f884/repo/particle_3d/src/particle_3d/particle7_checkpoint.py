"""Actual LAMMPS binary plus P6 schema/provenance and full P7 sidecar state."""
from pathlib import Path
from collections import deque,Counter
import json,shlex,hashlib,subprocess
from .audit import sha256,check_hash
from .particle3_cases import write_json
from .lammps_state import schema_contract
from .lammps_neighbors import ValidationNeighborPolicy
from .particle7_bridge import DynamicParticleBridge
from .injection_population import InjectionScheduler
from .particle7_lifecycle import PopulationLifecycle

SCHEMA='PARTICLE7_INJECTION_CHECKPOINT_V0'
REPO=Path(__file__).resolve().parents[3]


def scientific_identity():
    from .particle6_checkpoint import frozen_provenance
    paths=list((REPO/'particle_3d/src/particle_3d').glob('*.py'))+list((REPO/'particle_3d/contracts').glob('*.json'))
    return dict(frozen=frozen_provenance(REPO/'formal_3D_flow_solver/FEM_SimVascular'),
        model_source_sha256={str(p.relative_to(REPO)):sha256(p) for p in sorted(paths)})


def admission_identity(admission):
    def digest(a): return hashlib.sha256(a.tobytes()).hexdigest()
    return dict(class_name=type(admission).__name__,guard=admission.guard,width_m=getattr(admission,'width',None),
        velocity_m_s=admission.velocity.tolist(),inlet_triangles_sha256=digest(admission.sampler.triangles),
        inlet_normal_velocity_sha256=digest(admission.sampler.q),
        wall_triangles_sha256=None if admission.wall is None else digest(admission.wall.triangles))


def write_checkpoint(folder,engine,provenance):
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=False)
    b=engine.bridge
    if b is None: raise ValueError('Actual LAMMPS bridge required')
    engine.sync(); b.command('write_restart '+shlex.quote(str((folder/'state.restart').resolve())))
    write_json(folder/'sidecar.json',dict(schema=SCHEMA,p6_schema=schema_contract(),provenance=provenance,
        state=engine.state(),scientific_identity=scientific_identity(),admission_identity=admission_identity(engine.admission),
        sonovue_contract=engine.scheduler.source.mb_contract,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),used_ids=sorted(b.used_ids),policy=b.policy.to_dict(),box_half_width_m=b.box_half_width_m,lammps_version=b.version))
    write_json(folder/'manifest.json',dict(schema=SCHEMA,files={name:sha256(folder/name) for name in ['state.restart','sidecar.json']}))


def read_checkpoint(folder,admission_factory,classifier,provenance,*,mover=None):
    folder=Path(folder); m=json.loads((folder/'manifest.json').read_text())
    if m['schema']!=SCHEMA or set(m['files'])!={'state.restart','sidecar.json'}: raise ValueError('Checkpoint manifest schema mismatch')
    for name,digest in m['files'].items(): check_hash(folder/name,digest)
    d=json.loads((folder/'sidecar.json').read_text())
    if d['schema']!=SCHEMA or d['provenance']!=provenance or d['p6_schema']!=schema_contract(): raise ValueError('Checkpoint provenance/schema mismatch')
    if d['scientific_identity']!=scientific_identity(): raise ValueError('Checkpoint scientific source/Frozen integrity mismatch')
    s=d['state']; scheduler=InjectionScheduler.restore(s['scheduler']); admission=admission_factory(scheduler.source)
    if scheduler.source.mb_contract!=d['sonovue_contract'] or admission_identity(admission)!=d['admission_identity']: raise ValueError('Checkpoint inlet/admission/SonoVue model mismatch')
    raw=d['policy']; policy=ValidationNeighborPolicy(**{k:raw[k] for k in ['center_cutoff_m','skin_m','reason','role']})
    bridge=DynamicParticleBridge.from_binary_restart(folder/'state.restart',policy,d['box_half_width_m'])
    try:
        if bridge.version!=d['lammps_version']: raise ValueError('LAMMPS version mismatch')
        bridge.used_ids=set(d['used_ids']); obj=PopulationLifecycle(scheduler,admission,classifier,bridge=bridge,mover=mover)
        obj.time_s=s['time_s']; obj.active={p.particle_id:p for p in bridge.read()}
        from .particle3_cases import json_safe
        # Compare exact round-tripped binary/custom properties with saved Python records.
        actual=json.loads(json.dumps([p.to_dict() for p in sorted(obj.active.values(),key=lambda p:p.particle_id)],default=json_safe))
        if actual!=sorted(s['active'],key=lambda p:p['particle_id']): raise ValueError('Binary/sidecar active state mismatch')
        obj.pending={k:deque(v) for k,v in s['pending'].items()}; obj.births={int(k):v for k,v in s['births'].items()}
        for k in ['events','exits','timeline','neighbor_comparisons','neighbor_mismatches']: setattr(obj,k,s[k])
        admission.counts=Counter(s['admission_counts']); admission.candidates=s['admission_candidates']
        return obj
    except BaseException: bridge.close(); raise
