"""Fail-closed checkpoint contracts and absolute-time chains (CPU only on import).

Native files are archived only after normal process exit. An archive being intact
does not make missing SDPD RNG state recoverable. Diagnostic-only restoration is
an explicit experiment, never permission to join a production trajectory.
"""
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT, read_json, write_json, sha256_file, fingerprint


NATIVE_FILES = [
    'tests/restart/interactions.py', 'tests/restart/particle_vector.py',
    'src/mirheo/core/interactions/pairwise/sdpd.h', 'src/mirheo/core/interactions/pairwise/sdpd.cu',
    'src/mirheo/core/interactions/pairwise/base_pairwise.h',
    'src/mirheo/core/interactions/interface.h', 'src/mirheo/core/mirheo_object.cpp',
    'src/mirheo/core/interactions/pairwise/dpd.cu',
    'src/mirheo/core/interactions/pairwise/kernels/sdpd.h',
    'src/mirheo/core/interactions/pairwise/kernels/stress_wrapper.h',
    'src/mirheo/core/interactions/pairwise/stress.h',
    'src/mirheo/core/interactions/utils/step_random_gen.cpp',
    'src/mirheo/core/simulation.cpp', 'src/mirheo/core/task_scheduler.cpp','src/mirheo/core/task_scheduler.h',
    'src/mirheo/core/mirheo.cpp',
    'src/mirheo/core/pvs/particle_vector.cpp', 'src/mirheo/core/pvs/checkpoint/helpers.cpp',
    'src/mirheo/core/integrators/vv.cu', 'src/mirheo/plugins/particle_channel_saver.cpp',
    'src/mirheo/plugins/stats.cu','src/mirheo/plugins/virial_pressure.cu',
    'src/mirheo/bindings/mirheo.cpp']


def native_contract(root=None):
    root = Path(root or PROJECT_ROOT/'vendor/Mirheo')
    sources = {name:(root/name).read_text() for name in NATIVE_FILES}
    # This is an audit of a pinned source/build, not a general C++ parser. No
    # implementation is considered supported solely because text was added.
    no_sdpd_hooks = all(not re.search(r'\b(checkpoint|restart)\s*\(', sources[name]) for name in
                       ['src/mirheo/core/interactions/pairwise/sdpd.h',
                        'src/mirheo/core/interactions/pairwise/sdpd.cu',
                        'src/mirheo/core/interactions/pairwise/base_pairwise.h',
                        'src/mirheo/core/interactions/interface.h'])
    library=PROJECT_ROOT/'.venv/lib/python3.12/site-packages/libmirheo.cpython-312-x86_64-linux-gnu.so'
    symbols=subprocess.run(['nm','-C',str(library)],capture_output=True,text=True,check=True,timeout=15).stdout
    relevant=[line for line in symbols.splitlines() if re.search(r'(Pairwise(SDPD|DPD)Interaction|MirObject)::(checkpoint|restart)\(',line)]
    return {
        'status':'RESTART_NOT_VALIDATED', 'real_GPU_validation':'NOT_RUN',
        'native_rng_persistence_proven':False, 'formal_chain_allowed':False,
        'static_finding':'SDPD_INTERACTION_HAS_NO_CHECKPOINT_RESTART_OVERRIDE' if no_sdpd_hooks else 'SOURCE_REQUIRES_NEW_MANUAL_AUDIT',
        'reason':'PairwiseSDPD has writeState/readState, but PairwiseSDPDInteraction inherits no-op MirObject checkpoint/restart. The stress wrapper owns another copied kernel RNG; its state and StressManager scheduling are also unproven.',
        'native_seed':42424242, 'velocity_seed_is_a_different_generator':True,
        'official_examples':'interactions.py exercises DPD at dt=0; particle_vector.py exercises PV-only restart at dt=0. Neither validates nonzero-dt SDPD.',
        'state_contract':{
            'particles':'PV XMF/HDF5: global positions, velocities, int64 IDs and persistent channels; mass must be reconstructed from the exact frozen specification.',
            'coordinator':'u.restart restores simulation.state.txt: currentTime, currentStep, checkpointId. PV InitialConditions.Restart alone does not.',
            'dt':'Not in coordinator state file; exact frozen dt must be checked and passed to run.',
            'RNG':'Required: both active SDPD kernel generators, lastTime, lastSample, mt19937 state and stress scheduling. NOT PERSISTED by this SDPD interaction.',
            'integrator':'Current vv uses x and v plus freshly computed force; no hidden half-velocity serialized. Native checkpoint is after cell lists, before force clearing/integration of current_step.',
            'persistent_channels':['saved_density','saved_positions','saved_velocities','saved_forces','saved_stresses'],
            'ephemeral_channels':['densities','__forces','stresses'],
            'channel_phase':'Saved fields belong to the preceding force step; recompute ephemeral channels before advancing, do not report stale fields as a new measurement.',
            'plugins':'Same names/order. PostprocessStats restart copies historical native_stats.csv then appends; that file is secondary and never merged as new samples. VirialPressureDumper opens a fresh file on handshake. The Python observer exclusively creates raw_statistics.csv per segment and records restored boundary metadata separately.',
            'file_links':'MirObject createCheckpointSymlink actually uses ln -f (hard links). Materialize every file independently; record aliases and XMF references.',
            'save_phase':'No public manual checkpoint binding. Each Mirheo.run calls Simulation.init, recreating RunData and TaskScheduler(nExecutions=0). With checkpoint enabled, each 200-step chunk therefore saves at its start. Native scratch uses PingPong; boundary N is triggered by advancing once from N to N+1. Archive uses state.txt, never the directory label. Extra advances are charged and excluded from committed samples.'},
        'required_components':['time_state','particle_state','interaction_rng_state','frozen_parameters','integrator_phase','plugin_contract'],
        'sources':[{'path':str(root/name),'sha256':sha256_file(root/name),
                    'anchor_lines':[i+1 for i,line in enumerate(sources[name].splitlines()) if re.search(r'checkpoint|restart|writeState|readState|nExecutions_|make_unique<RunData>',line)]} for name in NATIVE_FILES],
        'compiled_library_symbol_audit':{'library':str(library),'sha256':sha256_file(library),
            'relevant_symbols':relevant,'SDPD_override_symbols_present':any('PairwiseSDPDInteraction::' in line for line in relevant),
            'scope':'CPU nm symbol inspection only; not GPU behavior validation.'},
        'selection':None}


def inside(root, relative):
    p = Path(relative)
    if p.is_absolute() or '..' in p.parts:
        raise ValueError('UNSAFE_CHECKPOINT_REFERENCE')
    resolved = (root/p).resolve(strict=True)
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError('CHECKPOINT_REFERENCE_ESCAPES_ARCHIVE')
    return resolved


def inventory(directory):
    root = Path(directory).resolve()
    hashes, references, aliases, inodes = {}, {}, {}, {}
    for p in sorted(root.rglob('*')):
        if p.is_dir():
            if p.is_symlink():raise ValueError('CHECKPOINT_DIRECTORY_SYMLINK_FORBIDDEN')
            continue
        rel = str(p.relative_to(root))
        resolved = inside(root,rel)
        if not resolved.is_file():raise ValueError('CHECKPOINT_NONREGULAR_FILE')
        hashes[rel] = sha256_file(resolved)
        stat = resolved.stat(); key = (stat.st_dev,stat.st_ino)
        if key in inodes:aliases[rel] = inodes[key]
        else:inodes[key] = rel
        if p.is_symlink():aliases[rel] = str(resolved.relative_to(root))
        if p.suffix.lower() in ['.xmf','.xdmf']:
            for node in ET.parse(p).iter('DataItem'):
                if node.get('Format','').upper() == 'HDF':
                    text = (node.text or '').strip()
                    if ':' not in text:raise ValueError('INVALID_XMF_HDF_REFERENCE')
                    ref = text.split(':',1)[0]
                    target = inside(root,str(p.parent.relative_to(root)/ref))
                    references.setdefault(rel,[]).append({'file':str(target.relative_to(root)), 'dataset':text.split(':',1)[1]})
    return hashes, references, aliases


def native_time(directory, dt):
    files = list(Path(directory).glob('*.state.txt')) + list(Path(directory).glob('state.txt'))
    if len(files) != 1:raise ValueError('MISSING_OR_AMBIGUOUS_COORDINATOR_TIME')
    values = files[0].read_text().split()
    if len(values) != 3:raise ValueError('INVALID_COORDINATOR_TIME')
    t, step, ident = float(values[0]), int(values[1]), int(values[2])
    if step < 0 or not math.isfinite(t) or not math.isclose(t,step*dt,abs_tol=2e-10,rel_tol=0):
        raise ValueError('CHECKPOINT_TIME_STEP_DT_MISMATCH')
    return {'absolute_step':step,'absolute_time_star':t,'checkpoint_id':ident,'time_file':files[0].name,'dt_star':dt}


def archive_checkpoint(source, destination, *, spec, contract, completion, parent=None):
    source, destination = Path(source), Path(destination)
    if completion.get('status') != 'COMPLETED_PLANNED_STEPS':
        raise ValueError('CHECKPOINT_WRITER_DID_NOT_EXIT_NORMALLY')
    if destination.exists():raise FileExistsError('IMMUTABLE_CHECKPOINT_ALREADY_EXISTS')
    before, references, aliases = inventory(source)
    if not any(p.endswith('.xmf') for p in before):raise ValueError('SNAPSHOT_IS_NOT_NATIVE_CHECKPOINT')
    state = native_time(source,spec['dt_star'])
    if state['absolute_step'] > completion['steps']:raise ValueError('CHECKPOINT_AFTER_ACTUAL_COMPLETION')
    destination.mkdir(parents=True)
    for name in before:
        p = destination/name;p.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(inside(source.resolve(),name),p)  # break hard links too
    if inventory(source)[0] != before or inventory(destination)[0] != before:
        raise ValueError('CHECKPOINT_CHANGED_DURING_COPY')
    manifest = {**state, 'kind':'NATIVE_CHECKPOINT_ARCHIVE', 'copy_complete':True,
                'files':before, 'references':references, 'aliases_materialized':aliases,
                'frozen_parameters_sha256':fingerprint(spec), 'mass_star':spec['m_star'],
                'parent_checkpoint':parent, 'source':str(source.resolve()),
                'source_contract_sha256':fingerprint(contract),
                'restart_validity':'RESTART_NOT_VALIDATED',
                'interaction_rng_state':None,
                'warning':'Integrity is not restart validity; current SDPD native RNG persistence is absent.'}
    write_json(destination/'checkpoint_manifest.json',manifest)
    # A separate final marker means partial copies cannot be used.
    write_json(destination/'COMPLETE.json',{'manifest_sha256':sha256_file(destination/'checkpoint_manifest.json')})
    return destination


def verify_checkpoint(directory, *, require_rng=True, expected_manifest=None):
    root = Path(directory).resolve()
    if not (root/'COMPLETE.json').is_file():raise ValueError('CHECKPOINT_INCOMPLETE_OR_SNAPSHOT_ONLY')
    digest = sha256_file(root/'checkpoint_manifest.json')
    if read_json(root/'COMPLETE.json')['manifest_sha256'] != digest or (expected_manifest and digest != expected_manifest):
        raise ValueError('CHECKPOINT_MANIFEST_HASH_MISMATCH')
    m = read_json(root/'checkpoint_manifest.json')
    if not m.get('copy_complete') or m.get('kind') != 'NATIVE_CHECKPOINT_ARCHIVE':raise ValueError('NOT_NATIVE_CHECKPOINT')
    for name,h in m['files'].items():
        try:p=inside(root,name)
        except (FileNotFoundError,RuntimeError) as exc:raise ValueError('CHECKPOINT_FILE_MISSING') from exc
        if sha256_file(p) != h:raise ValueError('CHECKPOINT_FILE_HASH_MISMATCH')
    actual, refs, _ = inventory(root)
    if set(actual)-{'COMPLETE.json','checkpoint_manifest.json'} != set(m['files']) or refs != m['references']:
        raise ValueError('CHECKPOINT_FILE_SET_OR_REFERENCE_CHANGED')
    if native_time(root,m['dt_star'])['absolute_step'] != m['absolute_step']:
        raise ValueError('CHECKPOINT_TIME_STATE_CHANGED')
    if require_rng and (not m.get('interaction_rng_state') or m.get('restart_validity') != 'PASS'):
        raise ValueError('RESTART_NOT_VALIDATED_MISSING_RANDOM_STATE')
    return m


def restore_coordinator(u, directory, *, diagnostic_only=False, expected_manifest=None, compute=True):
    m = verify_checkpoint(directory,require_rng=not diagnostic_only,expected_manifest=expected_manifest)
    require_readable_checkpoint_channels(directory)
    u.restart(str(directory))  # never draw/set velocities, remove COM or use PV-only restart here
    state = u.getState() if compute else None
    if state is not None and (state.current_step != m['absolute_step'] or
            not math.isclose(state.current_time,m['absolute_time_star'],abs_tol=2e-10,rel_tol=0)):
        raise ValueError('RESTORED_COORDINATOR_TIME_MISMATCH')
    return m


def checkpoint_channel_forms(directory):
    """CPU audit of persistent XMF forms accepted by this pinned native reader.

    This narrower check is neither an HDF integrity check nor RNG validation.
    In particular Force/Other is serialized with one component by this build;
    relabeling it as a vector cannot recover the missing components.
    """
    # Only particle channel forms used in this experiment are approved here.
    # Rigid-object/polyline layouts require their own native layout audit.
    forms = {'Scalar':1, 'Vector':3, 'Tensor6':6, 'Vector4':4}
    rows = []
    for path in sorted(Path(directory).glob('*.xmf')):
        for attribute in ET.parse(path).iter('Attribute'):
            info, data = attribute.find('Information'), attribute.find('DataItem')
            form = info.get('Value') if info is not None else None
            dims = [int(x) for x in data.get('Dimensions','').split()] if data is not None else []
            valid = form in forms and bool(dims) and dims[-1] == forms[form]
            rows.append({'xmf':path.name, 'name':attribute.get('Name'), 'form':form,
                         'datatype':info.get('Datatype') if info is not None else None,
                         'dimensions':dims, 'HDF_reference':(data.text or '').strip() if data is not None else None,
                         'status':'PASS' if valid else 'UNSUPPORTED_CHECKPOINT_CHANNEL'})
    return {'status':'FAIL' if any(r['status']!='PASS' for r in rows) else 'PASS',
            'scope':'Persistent XMF channel forms only; not full checkpoint or RNG validity.', 'channels':rows}


def require_readable_checkpoint_channels(directory):
    audit = checkpoint_channel_forms(directory)
    invalid = [row for row in audit['channels'] if row['status'] != 'PASS']
    if invalid:
        raise ValueError('UNSUPPORTED_CHECKPOINT_CHANNEL '+', '.join(
            f"{r['name']}={r['form']}/{r['datatype']} dims={r['dimensions']}" for r in invalid))
    return audit


def compare_by_id(a, b, domain, *, position_atol=2e-5, velocity_atol=2e-5):
    orders = []
    for x in [a,b]:
        ids=np.asarray(x['ids'])
        if ids.dtype.kind not in 'iu' or len(np.unique(ids)) != len(ids):raise ValueError('INVALID_PARTICLE_IDS')
        orders.append(np.argsort(ids))
    if not np.array_equal(np.asarray(a['ids'])[orders[0]],np.asarray(b['ids'])[orders[1]]):
        raise ValueError('PARTICLE_IDS_DIFFER')
    dx=np.asarray(a['positions'])[orders[0]]-np.asarray(b['positions'])[orders[1]]
    dx-=np.asarray(domain)*np.rint(dx/np.asarray(domain))
    dv=np.asarray(a['velocities'])[orders[0]]-np.asarray(b['velocities'])[orders[1]]
    if not np.isfinite(dx).all() or not np.isfinite(dv).all():raise ValueError('NONFINITE_RESTORED_STATE')
    return {'status':'PASS' if np.max(abs(dx))<=position_atol and np.max(abs(dv))<=velocity_atol else 'FAIL',
            'max_position_error_star':float(np.max(abs(dx))), 'max_velocity_error_star':float(np.max(abs(dv))),
            'RMS_velocity_error_star':float(np.sqrt(np.mean(dv**2))),
            'position_atol':position_atol,'velocity_atol':velocity_atol,'N':len(dx),'matched_by':'int64 particle ID'}


def require_formal_restart(contract, validation):
    if not contract.get('native_rng_persistence_proven') or validation.get('status') != 'PASS' or validation.get('evidence') != 'REAL_GPU_SEPARATE_PROCESSES':
        raise ValueError('RESTART_NOT_VALIDATED_NO_FORMAL_CHAIN')


def merge_segments(segments, dt, every, final_step):
    """Only completed, hash-verified, absolute-time segments; boundary row is metadata."""
    merged=[];last=0;parent=None
    for index, segment in enumerate(segments):
        if segment['status']!='COMPLETED' or segment['segment_index']!=index:
            raise ValueError('INCOMPLETE_OR_UNORDERED_SEGMENT')
        if segment['parent_checkpoint']!=parent or segment['start_step']!=last:
            raise ValueError('BROKEN_CHECKPOINT_OR_TIME_CHAIN')
        if segment['end_step']<=last or segment['end_step']-last!=segment['segment_steps']:
            raise ValueError('SEGMENT_ABSOLUTE_STEPS_MISMATCH')
        rows=segment['rows']
        expected=list(range(last if index==0 else last+every,segment['end_step']+1,every))
        # A restored t=0-like row must be marked boundary, never silently dropped.
        kept=[r for r in rows if not r.get('restore_boundary',False)]
        boundary=[r for r in rows if r.get('restore_boundary',False)]
        if boundary and (index==0 or len(boundary)!=1 or boundary[0]['step']!=last):raise ValueError('INVALID_BOUNDARY_ROW')
        if [r['step'] for r in kept]!=expected:raise ValueError('MISSING_DUPLICATED_OR_NONABSOLUTE_SAMPLES')
        if any(not math.isclose(r['time_star'],r['step']*dt,abs_tol=2e-10,rel_tol=0) for r in rows):
            raise ValueError('ROW_ABSOLUTE_TIME_MISMATCH')
        merged.extend(kept);last=segment['end_step'];parent=segment['checkpoint_sha256']
    if last!=final_step:raise ValueError('FORMAL_WINDOW_NOT_REACHED')
    return merged
