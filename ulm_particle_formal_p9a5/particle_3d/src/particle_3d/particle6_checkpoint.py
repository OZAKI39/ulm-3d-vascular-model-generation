"""Binary LAMMPS state plus verified global sidecar. Never regenerate particles."""
from pathlib import Path
import json,shlex,subprocess
from .audit import sha256,check_hash
from .particle3_cases import write_json
from .lammps_state import SCHEMA_VERSION,schema_contract
from .lammps_bridge import LammpsParticleBridge
from .lammps_neighbors import ValidationNeighborPolicy
from .particle6_stepper import Particle6Stepper

DEPENDENCIES=['7cfe5141382600e28582f05bff712a6f09c38a39','6e59c605cb8251b9fa2e80d2dbed0cfa6daaf03c',
'c704e08d39f6134da300fc91cc93c8561314a17b','95e54fe474349c35aaad2e2754aff0bad7108c42',
'a1e09cb3e2db7c91d7adc75adf40a0f939a880bf','adbde27d31824e1a3ec5c47c9b297e1c0a248d6d']
MODEL_FLAGS=dict(lammps_role='STATE_NEIGHBOR_CHECKPOINT_ONLY',physics_authority='PARTICLE_3D',
 lammps_force_integration=False,lammps_time_integration=False,production_neighbor_cutoff_frozen=False,
 production_neighbor_skin_frozen=False,production_lubrication_cutoff_frozen=False,production_particle_timestep_frozen=False,
 nonspherical_lubrication='NOT_FROZEN',real_rbc_lammps_dynamics='DEFERRED_DUE_TO_UPSTREAM_PHYSICS',
 continuum_validity_at_sub_nanometer_gap='NOT_ESTABLISHED',no_cfd_executed=True,particle7_started=False)


def frozen_provenance(fem):
    fem=Path(fem);hashes={}
    for line in (fem/'frozen_reference/SHA256SUMS.txt').read_text().splitlines():
        digest,name=line.split('  ',1);check_hash(fem/name,digest);hashes[name]=digest
    return dict(frozen_fem_commit='c83dccea0f1fe5cd9fd3f5939e2eba883aeba9c2',frozen_sha256=hashes,
        dependency_commits={f'particle{i}':c for i,c in enumerate(DEPENDENCIES)})


def source_identity(repo):
    repo=Path(repo);root=repo/'particle_3d';paths=[]
    for folder in ['src','scripts','tests','contracts']:
        paths.extend(p for p in (root/folder).rglob('*') if p.is_file() and p.suffix in ['.py','.json'])
    paths.append(root/'PARTICLE6_README.md')
    return dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        source_sha256={str(p.relative_to(repo)):sha256(p) for p in sorted(paths)},
        source_commit_status='WORKTREE_SHA_SNAPSHOT; BIND_EXACT_COMMIT_AFTER_FULL_TESTS')


def write_checkpoint(folder,stepper,provenance,repo):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    if stepper.bridge is None:raise ValueError('Checkpoint requires actual LAMMPS storage')
    bridge=stepper.bridge;binary=folder/'state.restart';sidecar=folder/'particle_sidecar.json'
    bridge.rebuild();bridge.command('write_restart '+shlex.quote(str(binary.resolve())))
    global_state=dict(stepper.global_metadata(),schema_version=SCHEMA_VERSION,custom_property_schema=schema_contract(),
        provenance=provenance,model_flags=MODEL_FLAGS,validation_neighbor_policy=bridge.policy.to_dict(),
        box_half_width_m=bridge.box_half_width_m,dynamic_viscosity_pa_s=stepper.mu,particle_count=len(bridge.read()),
        stable_ids=[p.particle_id for p in bridge.read()],lammps_version=bridge.version)
    write_json(sidecar,global_state)
    manifest=dict(schema_version=SCHEMA_VERSION,lammps_version=bridge.version,**source_identity(repo),
        files={binary.name:sha256(binary),sidecar.name:sha256(sidecar)})
    write_json(folder/'checkpoint_manifest.json',manifest)
    return manifest


def read_checkpoint(folder,provenance,velocity_provider,mu,*,wall=None,boundary_classifier=None):
    folder=Path(folder);manifest=json.loads((folder/'checkpoint_manifest.json').read_text())
    if manifest['schema_version']!=SCHEMA_VERSION or set(manifest['files'])!={'state.restart','particle_sidecar.json'}:
        raise ValueError('CHECKPOINT_SCHEMA_MISMATCH')
    for name,digest in manifest['files'].items():check_hash(folder/name,digest)
    sidecar=json.loads((folder/'particle_sidecar.json').read_text())
    if sidecar['schema_version']!=SCHEMA_VERSION or sidecar['custom_property_schema']!=schema_contract():raise ValueError('CHECKPOINT_PROPERTY_SCHEMA_MISMATCH')
    if sidecar['provenance']!=provenance or sidecar['model_flags']!=MODEL_FLAGS:raise ValueError('CHECKPOINT_PHYSICS_PROVENANCE_MISMATCH')
    if sidecar['dynamic_viscosity_pa_s']!=mu:raise ValueError('CHECKPOINT_VISCOSITY_MISMATCH')
    raw=sidecar['validation_neighbor_policy'];policy=ValidationNeighborPolicy(**{k:raw[k] for k in ['center_cutoff_m','skin_m','reason','role']})
    bridge=LammpsParticleBridge.from_binary_restart(folder/'state.restart',policy,sidecar['box_half_width_m'])
    try:
        if bridge.version!=manifest['lammps_version'] or bridge.version!=sidecar['lammps_version']:raise ValueError('CHECKPOINT_LAMMPS_VERSION_MISMATCH')
        particles=bridge.read()
        if len(particles)!=sidecar['particle_count'] or [p.particle_id for p in particles]!=sidecar['stable_ids']:raise ValueError('CHECKPOINT_STABLE_ID_MISMATCH')
        return Particle6Stepper(particles,policy,velocity_provider,mu,bridge=bridge,model=sidecar['model'],wall=wall,
            boundary_classifier=boundary_classifier,physical_time_s=sidecar['physical_time_s'],step_index=sidecar['particle_step_index'],
            boundary_event=sidecar['boundary_event'],boundary_events={int(k):v for k,v in sidecar['boundary_events'].items()})
    except BaseException:bridge.close();raise
