"""Recover serialization from preserved records, without running a scientific solver."""
import ast
import copy
import hashlib
import json
from pathlib import Path

A = Path(__file__).resolve().parent
ROOT = A.parents[2]
OLD = A.parent / 'repository_cleanup_v1'
ROWS = json.loads((OLD / 'RECOVERY_UNAVAILABLE.json').read_text())
INVENTORY = {r['path']: r for r in json.loads((OLD / 'INVENTORY_RAW.json').read_text())}
TEMPLATE_ROOT = Path('/tmp/clot_repository_cleanup_v1_fnxy7lzr/final_pytest')
RECORDS = []

def sha(data):
    return hashlib.sha256(data).hexdigest()

def verified_read(path, evidence):
    data = path.read_bytes()
    rel = str(path.relative_to(ROOT))
    assert sha(data) == INVENTORY[rel]['hash'], rel
    evidence.append(dict(path=rel, sha256=sha(data)))
    return data

def read_json(path, evidence):
    return json.loads(verified_read(path, evidence))

def save_stage(row, value):
    path = A / 'semantic_recovery_staging' / row['path']
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    return str(path.relative_to(A))

def recover_identity(row):
    target = ROOT / row['path']
    evidence = []
    config = read_json(target.parent / 'CONFIG.json', evidence)
    template = json.loads((TEMPLATE_ROOT / 'test_static_resolved_field_use0/resolved_run/IDENTITY.json').read_text())
    template['config_sha256'] = sha(json.dumps(config, sort_keys=True).encode())
    template['source_sha256'] = {k: sha(verified_read(target.parent / 'source_snapshot' / k, evidence))
                                 for k in template['source_sha256']}
    data = (json.dumps(template, indent=2, allow_nan=False) + '\n').encode()
    # Expected historical SHA is the acceptance criterion, including all metadata.
    assert sha(data) == row['hash'] and len(data) == row['size_bytes'], row['path']
    staged = A / 'recovery_staging' / row['hash']
    staged.parent.mkdir(exist_ok=True)
    staged.write_bytes(data)
    assert sha(staged.read_bytes()) == row['hash']
    with target.open('xb') as out:
        out.write(staged.read_bytes())
    assert sha(target.read_bytes()) == row['hash']
    record = dict(path=row['path'], sha256=row['hash'], status='EXACT_RECONSTRUCTED_FROM_HISTORICAL_CONFIG_AND_SOURCE',
                  staged=str(staged.relative_to(A)), input_evidence=evidence,
                  template_role='serialization order and invariant fields only; full original SHA verified')
    with (A / 'EXACT_RECOVERY_LOG.jsonl').open('a') as out:
        out.write(json.dumps(record) + '\n')
    return record

def recover_adr_summary(row):
    folder = (ROOT / row['path']).parent
    evidence = []
    config = read_json(folder / 'CONFIG.json', evidence)
    history = read_json(folder / 'history.json', evidence)
    jumps = read_json(folder / 'damage_jumps.json', evidence)
    equilibrations = []
    for i, path in enumerate(sorted(folder.glob('equilibrium_*/SUMMARY.json'))):
        d = read_json(path, evidence)
        equilibrations.append(dict(index=i, exposure_s=d['erosion_exposure_time_s'], phase=d['phase'],
                                   converged=d['converged'], iterations=d['iterations'],
                                   relative_force_residual=d['final']['relative_force_residual']))
    source = folder / 'source_snapshot/pd_clot/regularization/erosion_equilibrium_runner.py'
    tree = ast.parse(verified_read(source, evidence))
    expression = next(n.value for n in ast.walk(tree) if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == 'summary' for t in n.targets))
    expression = copy.deepcopy(expression)
    omitted = ['wall_time_s', 'cpu_time_s', 'peak_process_rss_bytes']
    expression.keywords = [k for k in expression.keywords if k.arg not in omitted]
    failure = ('INCOMPLETE_MECHANICAL_EQUILIBRIUM' if not equilibrations[-1]['converged'] else
               'INCOMPLETE_TOPOLOGY_OUTSIDE_EQUILIBRIUM_SCOPE' if any(j['broken_bond_count'] for j in jumps) else None)
    scope = dict(history=history, jumps=jumps, equilibrations=equilibrations, failure=failure,
                 count=config['streaming_erosion']['number_of_erosion_steps'],
                 exposure=jumps[-1]['erosion_exposure_end_s'] if jumps else 0.,
                 energy_negative=any(j['delta_U_damage_J'] < 0 for j in jumps))
    # Only a pure dict expression from the preserved historical writer is evaluated.
    # No function call into mechanics, damage, or a solver is made.
    value = eval(compile(ast.Expression(expression), str(source), 'eval'),
                 {'__builtins__': {}, 'dict': dict, 'len': len, 'sum': sum}, scope)
    fixture = Path(*folder.parts[folder.parts.index('tests') + 2:])
    fresh = json.loads((TEMPLATE_ROOT / fixture / 'SUMMARY.json').read_text())
    fresh_projection = {k: v for k, v in fresh.items() if k not in omitted}
    assert set(value) == set(fresh_projection)
    # Agreement with an already-existing validation output is supplementary, not the
    # source of the recovered values. Values above come only from hashed old records.
    agreement = value == fresh_projection
    assert agreement, row['path']
    staged = save_stage(row, value)
    return dict(path=row['path'], expected_sha256=row['hash'], classification='SCIENTIFIC_EQUIVALENT',
                original_path_restored=False, staged=staged, scientific_keys_verified=list(value),
                scientific_key_count=len(value), omitted_non_scientific_keys=omitted,
                full_container_semantic_exact=False, input_evidence=evidence,
                comparison=dict(scientific_projection_equal_to_existing_validation_output=agreement,
                                historical_inputs_sha_verified=True, new_simulations=0),
                old_hash_test_satisfied=False,
                limitation='Original timing/RSS unknown; stage does not replace old hash dependency.')

def partial_streaming_summary(row):
    folder = (ROOT / row['path']).parent
    evidence = []
    history = read_json(folder / 'history.json', evidence)
    components = read_json(folder / 'components.json', evidence)
    settling = read_json(folder / 'settling_windows.json', evidence)
    writer = folder / 'source_snapshot/pd_clot/regularization/erosion_runner.py'
    verified_read(writer, evidence)
    value = dict(final=history[-1], maximum_driver_relative_drift=max(h['driver_relative_drift'] for h in history[1:]),
                 damage_energy_status=history[-1]['damage_energy_status'])
    stage = save_stage(row, value)
    return dict(path=row['path'], expected_sha256=row['hash'], classification='NOT_EQUIVALENT',
                original_path_restored=False, staged=stage, verified_partial_keys=list(value),
                input_evidence=evidence, new_simulations=0, old_hash_test_satisfied=False,
                unrecoverable_scientific_fields=['minimum_J_or_intrinsic_stretch', 'maximum_acceleration_m_s2'],
                unrecoverable_telemetry=['elapsed_wall_s'],
                reason='Global extrema updated at every mechanical step; saved states and settling samples do not retain every step. Reusing a fresh run would not verify historical extrema.')

if __name__ == '__main__':
    for row in ROWS:
        if (ROOT / row['path']).exists():
            continue
        if row['path'].endswith('/IDENTITY.json'):
            RECORDS.append(recover_identity(row))
        elif row['path'].endswith('/SUMMARY.json'):
            if ((ROOT / row['path']).parent / 'damage_jumps.json').exists():
                RECORDS.append(recover_adr_summary(row))
            else:
                RECORDS.append(partial_streaming_summary(row))
    (A / 'SEMANTIC_RECOVERY_COMPARISON.json').write_text(json.dumps(RECORDS, indent=2) + '\n')
    from collections import Counter
    print(json.dumps(dict(Counter(r.get('classification', r.get('status')) for r in RECORDS))), flush=True)
