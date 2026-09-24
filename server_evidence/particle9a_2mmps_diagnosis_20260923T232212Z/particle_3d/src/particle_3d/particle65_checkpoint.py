"""V1-specific restart metadata over the unchanged P6 binary storage bridge."""
from pathlib import Path
import json,shlex
from .audit import sha256,check_hash
from .particle3_cases import write_json
from .particle6_checkpoint import source_identity
from .lammps_state import SCHEMA_VERSION,schema_contract
from .lammps_neighbors import ValidationNeighborPolicy
from .lammps_bridge import LammpsParticleBridge
from .nearfield_regularization import NearFieldRegularizationV1,contract
from .particle65_motion import Particle65Stepper


def write_checkpoint(folder,stepper,provenance,repo):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
    if stepper.bridge is None:raise ValueError('Actual LAMMPS bridge required')
    bridge=stepper.bridge;bridge.rebuild();binary=folder/'state.restart'
    bridge.command('write_restart '+shlex.quote(str(binary.resolve())))
    side=dict(stepper.global_metadata(),contract=contract(),schema_version=SCHEMA_VERSION,custom_property_schema=schema_contract(),
        provenance=provenance,validation_neighbor_policy=bridge.policy.to_dict(),box_half_width_m=bridge.box_half_width_m,
        dynamic_viscosity_pa_s=stepper.mu,stable_ids=[p.particle_id for p in bridge.read()],lammps_version=bridge.version)
    write_json(folder/'particle_sidecar.json',side)
    manifest=dict(**source_identity(repo),files={name:sha256(folder/name) for name in ['state.restart','particle_sidecar.json']})
    write_json(folder/'checkpoint_manifest.json',manifest)
    return manifest


def read_checkpoint(folder,provenance,provider,mu,*,wall=None,boundary_classifier=None):
    folder=Path(folder);m=json.loads((folder/'checkpoint_manifest.json').read_text())
    if set(m['files'])!={'state.restart','particle_sidecar.json'}:raise ValueError('CHECKPOINT_FILE_SCHEMA_MISMATCH')
    for name,digest in m['files'].items():check_hash(folder/name,digest)
    d=json.loads((folder/'particle_sidecar.json').read_text())
    if d['contract']!=contract() or d['provenance']!=provenance or d['dynamic_viscosity_pa_s']!=mu:raise ValueError('CHECKPOINT_V1_PROVENANCE_MISMATCH')
    if d['schema_version']!=SCHEMA_VERSION or d['custom_property_schema']!=schema_contract():raise ValueError('CHECKPOINT_STATE_SCHEMA_MISMATCH')
    raw=d['validation_neighbor_policy'];query=ValidationNeighborPolicy(**{k:raw[k] for k in ['center_cutoff_m','skin_m','reason','role']})
    policy=NearFieldRegularizationV1(d['h_molecular_floor_m'],d['parameter_role'])
    bridge=LammpsParticleBridge.from_binary_restart(folder/'state.restart',query,d['box_half_width_m'])
    try:
        if bridge.version!=d['lammps_version'] or [p.particle_id for p in bridge.read()]!=d['stable_ids']:raise ValueError('CHECKPOINT_IDS_OR_VERSION_MISMATCH')
        result=Particle65Stepper(bridge.read(),query,provider,mu,wall=wall,bridge=bridge,policy=policy,boundary_classifier=boundary_classifier,
            physical_time_s=d['physical_time_s'],step_index=d['particle_step_index'])
        result.boundary_event=d['boundary_event'];result.boundary_events={int(k):v for k,v in d['boundary_events'].items()}
        return result
    except BaseException:bridge.close();raise
