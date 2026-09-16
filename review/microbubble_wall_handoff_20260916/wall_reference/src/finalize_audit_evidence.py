"""Read-only numerical/source finalizer; never rerun or repair a package result.

Verifies the independent SI conversion, serialization, first-output identity,
source checkouts/build inputs, and original baseline manifests. Generates only
new audit evidence. Run before report generation and final SHA256SUMS.
"""
from pathlib import Path
import csv, hashlib, json, subprocess
import numpy as np
import h5py

R = Path(__file__).resolve().parents[1]
W = Path(json.loads((R/'provenance/TASK_PATHS.json').read_text())['work'])

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def save(name, obj):
    (R/name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def csvout(name, rows):
    with (R/name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def rel(a, b):
    return float(np.max(abs(a-b)/np.maximum(np.maximum(abs(a), abs(b)), 1e-30)))

def main():
    before_raw = {str(p.relative_to(R)):sha(p) for p in sorted((R/'raw').glob('*'))}
    freeze_path = R/'provenance/FIRST_OUTPUT_IDENTITY.json'
    if freeze_path.exists():
        assert json.loads(freeze_path.read_text())['sha256'] == before_raw
    else:
        save('provenance/FIRST_OUTPUT_IDENTITY.json', {'scope':'first formal sweep and first historical diagnostic outputs; no repeat', 'sha256':before_raw})
    initial = json.loads((R/'provenance/INITIAL_SOURCE_IDENTITY.json').read_text())
    source = {}
    for project, identity in initial.items():
        root = W/'upstream'/project
        git = lambda *args: subprocess.check_output(['git','-C',str(root),*args], text=True).strip()
        checked = 0; bad = []; build_bad = []
        for line in (R/'provenance'/(project+'_SOURCE_SHA256SUMS')).read_text().splitlines():
            digest, name = line.split(None,1); name = name.lstrip('* '); checked += 1
            if not (root/name).is_file() or sha(root/name) != digest: bad.append(name)
            if project == 'pystokes':
                p = W/'build/pystokes'/name
                if not p.is_file() or sha(p) != digest: build_bad.append(name)
        source[project] = {'commit':git('rev-parse','HEAD'), 'expected_commit':identity['commit'],
            'git_status_porcelain':git('status','--porcelain'), 'tracked_files_checked':checked,
            'tracked_hash_mismatches':bad, 'build_copy_tracked_mismatches':build_bad,
            'submodules':git('submodule','status'), 'status':'PASS'}
        assert source[project]['commit'] == identity['commit']
        assert not bad and not build_bad and not source[project]['git_status_porcelain']
    assert sha(W/'build/rmbw/Lubrication_Class.so') == sha(R/'provenance/Lubrication_Class.so')
    pprov = json.loads((R/'PYSTOKES_PROVENANCE.json').read_text())
    for name, digest in pprov['installed_native_extensions'].items():
        assert sha(W/'envs/wallref_pystokes/lib/python3.12/site-packages/pystokes'/name) == digest
    save('provenance/SOURCE_INTEGRITY_FINAL.json', {'status':'PASS','sources':source,'compiled_bindings_identity':'PASS'})

    baselines = [json.loads((R/'provenance'/('BASELINE_IDENTITY_'+suffix+'.json')).read_text()) for suffix in ['BEFORE','AFTER']]
    assert baselines[0] == baselines[1] and baselines[0]['status'] == 'PASS'
    labels = ['RMBW_SPHERE','RMBW_LUBRICATION','PYSTOKES']
    original = [json.loads(line) for label in labels for line in (R/'raw'/(label+'.jsonl')).read_text().splitlines()]
    convention = json.loads((R/'WALL_REFERENCE_CONVENTION.json').read_text())
    expected = {(label,size,e,mu) for label in labels for size in convention['radii_m'] for e in convention['gaps_h_over_a'] for mu in convention['viscosities_Pa_s']}
    keys = [(x['implementation'],x['size'],x['epsilon'],x['viscosity_Pa_s']) for x in original]
    assert set(keys) == expected and len(keys) == len(expected) == 378
    assert all(x['first_formal_output'] and x['finite'] for x in original)
    warning_rows = [{k:x[k] for k in ['implementation','size','epsilon','viscosity_Pa_s','warnings']}
                    for x in original if x['warnings']]
    tables = list(csv.DictReader((R/'WALL_REFERENCE_RAW_RESULTS.csv').open()))
    errors = {'rmbw_native_to_SI_resistance':0., 'pystokes_independent_self_formula':0., 'scaled_resistance_identity':0.}
    with h5py.File(R/'WALL_REFERENCE_MOBILITY_MATRICES.h5') as hm, h5py.File(R/'WALL_REFERENCE_RESISTANCE_MATRICES.h5') as hr:
        assert hm.attrs['convention_sha256'] == sha(R/'WALL_REFERENCE_CONVENTION.json')
        assert hr.attrs['convention_sha256'] == hm.attrs['convention_sha256']
        M = hm['matrix_SI'][:]; B = hm['work_conjugate_normalized_matrix'][:]
        Rs = hr['matrix_SI'][:]; C = hr['work_conjugate_normalized_matrix'][:]
        assert np.array_equal(M,np.array([x['matrix_SI'] for x in original]))
        for k, row in enumerate(original):
            a, mu = row['radius_m'], row['viscosity_Pa_s']
            assert row['radius_m'] == convention['radii_m'][row['size']]
            assert row['gap_m'] == a*row['epsilon']
            scale = np.sqrt([6*np.pi*mu*a]*3+[8*np.pi*mu*a**3]*3)
            assert np.allclose(B[k], M[k]*np.outer(scale,scale),rtol=2e-15,atol=0)
            residual = float(np.max(abs(C[k]@B[k]-np.eye(6))))
            errors['scaled_resistance_identity'] = max(errors['scaled_resistance_identity'],residual)
            assert residual < 1e-12
            assert bool(hm['physical_matrix_valid'][k]) == (row['implementation']=='RMBW_LUBRICATION')
            assert tables[k]['implementation'] == row['implementation']
            if row['implementation'] == 'RMBW_LUBRICATION':
                meta = row['metadata']; L = meta['native_length_unit_m']; eta0 = meta['native_viscosity_unit_Pa_s']
                native = np.array(meta['excess_resistance_native'])+np.array(meta['bulk_resistance_native'])
                # Resistance units derived directly from F/V, F/Omega, T/V, T/Omega.
                unit = eta0*L*np.outer([1,1,1,L,L,L],[1,1,1,L,L,L])
                err = rel(Rs[k],native*unit)
                errors['rmbw_native_to_SI_resistance'] = max(errors['rmbw_native_to_SI_resistance'],err)
                assert err < 1e-12 and not meta['debye_cut_active']
                assert not meta['sign_correction_applied'] and not meta['matrix_symmetrized']
            elif row['implementation'] == 'PYSTOKES':
                # Independent expressions read from frozen wallBounded.pyx self terms.
                s = a/(a+row['gap_m']); mt = 1/(6*np.pi*mu*a); mr = 1/(8*np.pi*mu*a**3)
                predicted = np.diag([mt*(1-9*s/16+s**3/8-s**5/16)]*2+
                    [mt*(1-9*s/8+s**3/2-s**5/8)]+[mr*(1-5*s**3/16)]*2+[mr*(1-s**3/8)])
                c = s**4/(64*np.pi*mu*a*a)
                predicted[0,4]=-c; predicted[1,3]=c; predicted[4,0]=c; predicted[3,1]=-c
                err = rel(predicted,M[k]); errors['pystokes_independent_self_formula'] = max(errors['pystokes_independent_self_formula'],err)
                assert err < 1e-12
        for filename, values, normalized in [('WALL_MOBILITY_MATRIX.csv',M,B),('WALL_RESISTANCE_MATRIX.csv',Rs,C)]:
            rows = list(csv.DictReader((R/filename).open())); assert len(rows) == len(original)*36
            for idx,r in enumerate(rows):
                k,rem = divmod(idx,36); i,j = divmod(rem,6)
                assert int(r['row'])==i and int(r['column'])==j and r['implementation']==original[k]['implementation']
                assert float(r['value_SI']) == values[k,i,j] and float(r['normalized']) == normalized[k,i,j]
    vis = json.loads((R/'VISUALIZATION_PROVENANCE.json').read_text())
    assert sha(R/vis['generator']) == vis['generator_sha256']
    for name,digest in vis['input_files'].items(): assert sha(R/name)==digest
    for item in vis['artifacts']: assert sha(R/item['file'])==item['sha256']
    csvout('validation/EXTREME_GAP_STABILITY.csv',[r for r in tables if float(r['epsilon'])<=.005])
    csvout('validation/FAR_FIELD_TREND.csv',[r for r in tables if float(r['epsilon'])>=5 and r['size']=='d50' and float(r['viscosity_Pa_s'])==.001])
    after_raw = {str(p.relative_to(R)):sha(p) for p in sorted((R/'raw').glob('*'))}
    assert before_raw == after_raw
    save('validation/INDEPENDENT_FINALIZER.json', {'status':'PASS', 'package_runs_started':0,
        'first_raw_outputs_unchanged':True, 'unique_formal_cases':len(original), 'matrix_entries_per_output':len(original)*36,
        'all_formal_rows_finite':True, 'warning_rows':warning_rows, 'CSV_HDF5_raw_identity':'PASS',
        'convention_hash_identity':'PASS', 'SI_conversion_independent_checks':errors,
        'source_and_binding_integrity':'PASS','baseline_before_after_identity':'PASS',
        'baseline_manifest_files_checked':sum(x['checked_files'] for x in baselines[0]['baselines']),
        'visualization_inputs_identity':'PASS', 'scientific_status':'See primary selection; finalizer PASS does not certify all software matrices'})
    print(json.dumps({'FINALIZER':'PASS','errors':errors,'formal_cases':len(original)},indent=2))

if __name__=='__main__': main()
