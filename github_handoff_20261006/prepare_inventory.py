"""Read-only scientific inventory; write handoff metadata only. Never stage/delete."""
from pathlib import Path
import csv, hashlib, json, subprocess

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / 'github_handoff_20261006'
STAGES = [
    'pfpd_baqp_finite_domain_boundary_consistency_v1',
    'pfpd_baqp_fully_closed_force_stencil_audit_v1',
    'pfpd_baqp_affine_coefficient_reference_audit_v1',
    'nonlocal_boundary_formulation_v1',
    'nonlocal_boundary_nonaffine_localization_v1',
    'nonlocal_boundary_nonaffine_shear_localization_v1',
    'baqp_near_boundary_kinematic_energy_audit_v1',
    'higher_order_baqp_rk_oracle_v1',
    'higher_order_baqp_energy_force_candidate_v1',
    'higher_order_baqp_boundary_reintegration_v1',
]
FIGURES = {
    'higher_order_baqp_rk_oracle_v1': ['c1_vs_pd2_boundary_bias.png', 'pd2_polynomial_reproduction.png'],
    'higher_order_baqp_energy_force_candidate_v1': ['c1_variational_force_regression.png', 'pd2_energy_force_directional_derivative.png'],
    'higher_order_baqp_boundary_reintegration_v1': ['c1_vs_pd2_boundary_field_quality.png', 'pd2_manufactured_recovery.png'],
}
TEXT = {'.py', '.jl', '.md', '.json', '.csv', '.txt', '.log', '.yaml', '.yml', '.toml', '.patch', '.sh'}
SKIP_PARTS = {'__pycache__', '.pytest_cache', '.mypy_cache', '.ruff_cache', 'mpl_cache', 'mpl', '.git', 'external', 'node_modules', '.venv', 'venv', 'env'}
SUPERSEDED_LOGS = {
    'TARGETED_TESTS_INITIAL_API_TYPO.log', 'SOLVES_ATTEMPT_1_INTERRUPTED.log',
    'TARGETED_TESTS_INITIAL.log', 'TARGETED_TESTS_PRE_SOLVE.log',
    'FOCUSED_TESTS_PRE_SOLVE.log', 'PLOT_REFERENCES.log',
    'PLOT_REFERENCES_LAYOUT_REVIEW.log', 'PLOT_RECOVERY.log',
    'DECOMPOSITION_REPRODUCIBLE_GRID.log', 'N27_VERIFIED_PREREQUISITES.log',
    'neighbor_redraw.log',
}

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def tracked():
    return set(subprocess.check_output(['git', 'ls-files', '-z'], cwd=REPO).decode().split('\0')) - {''}

def select(p, inherited):
    rel = p.relative_to(REPO).as_posix(); size = p.stat().st_size
    stage = next((s for s in STAGES if s in p.parts), None)
    replacement = f'clot_clearance/verification/{stage}/' if stage else 'clot_clearance/analysis/'
    if inherited:
        return 'INHERITED', 'Already tracked; retained unchanged by selective handoff', ''
    if p.suffix in {'.py', '.jl', '.sh'}:
        return 'INCLUDE', 'Current scientific source, test, or reproducibility script', ''
    if p.name in FIGURES.get(stage, []):
        assert size <= 3*2**20
        return 'INCLUDE', 'One of two representative PNGs for this current PD2 stage', ''
    if p.name in SUPERSEDED_LOGS or (stage=='nonlocal_boundary_nonaffine_localization_v1' and p.name in {'ORACLE_RUN.log','REPORT_BUILD.log'}):
        return 'EXCLUDE', 'Resolved retry/debug log; final test/run/report logs retained', replacement + 'final REVIEW and final logs'
    if p.suffix == '.npz':
        if stage == 'higher_order_baqp_energy_force_candidate_v1' and size <= 2*2**20:
            return 'INCLUDE', 'Compact accepted state / FD direction / N7 test fixture needed by unchanged focused tests', ''
        if stage == 'higher_order_baqp_boundary_reintegration_v1' and p.name.endswith('_RECOVERED_STATE.npz'):
            return 'INCLUDE', 'One compact accepted solution for each of the five final cases; omit repeated field arrays and solve inputs', ''
        return 'EXCLUDE', 'Raw array or reference cache; compact accepted metrics retained; local original untouched', replacement + 'DATA_DICTIONARY.md and final JSON/CSV tables'
    if p.suffix in {'.png', '.pdf', '.svg', '.html', '.pptx'}:
        return 'EXCLUDE', 'Redundant visual format, historical visualization, or downloaded paper; only six current PNGs selected', replacement + ('references/SOURCE_PROVENANCE.json' if 'references' in p.parts else 'README/REVIEW and compact tables')
    if 'references' in p.parts and p.suffix == '.txt':
        return 'EXCLUDE', 'Extracted full paper text; provenance and scoped formula attribution retained instead', replacement + 'references/SOURCE_PROVENANCE.json'
    if 'assembly_implementation_check' in p.parts and p.suffix != '.md':
        return 'EXCLUDE', 'Superseded implementation-check evidence; final implementation results retained', replacement + 'TARGETED_TESTS.log and final summary tables'
    if p.suffix in {'.csv', '.json', '.log', '.txt'} and size > 2**20:
        return 'EXCLUDE', 'Large per-node table or inventory; compact scientific aggregates retained', replacement + 'DATA_DICTIONARY.md / summary JSON and github_handoff_20261006/LARGE_TABLE_AGGREGATES.json'
    if p.suffix in TEXT:
        return 'INCLUDE', 'Scientific report, configuration, compact data, final log, or provenance', ''
    return 'EXCLUDE', 'Nonessential auxiliary artifact; scientific source and final evidence retained', replacement

def main():
    known = tracked(); candidates = set()
    for name in ['pd_clot', 'scripts', 'tests', 'configs', 'analysis', 'vendor']:
        base = REPO/'clot_clearance'/name
        for p in base.rglob('*'):
            if p.is_file() and not any(x in SKIP_PARTS for x in p.relative_to(base).parts):
                if p.suffix in TEXT | {'.npz', '.png', '.pdf', '.svg', '.html', '.pptx'}:
                    candidates.add(p)
    for kind in ['verification', 'results']:
        for stage in STAGES:
            base=REPO/'clot_clearance'/kind/stage
            for p in base.rglob('*'):
                if p.is_file() and not any(x in SKIP_PARTS for x in p.relative_to(base).parts): candidates.add(p)
    rows=[]
    for p in sorted(candidates):
        rel=p.relative_to(REPO).as_posix(); decision,reason,replacement=select(p,rel in known)
        rows.append(dict(path=rel,size=p.stat().st_size,type=p.suffix,sha256=sha(p),decision=decision,reason=reason,compact_replacement=replacement))
    fields=list(rows[0])
    with (OUT/'SCIENTIFIC_FILE_INVENTORY.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
    # Frozen input inventory is created once and never overwritten by later handoff work.
    freeze=OUT/'PRE_HANDOFF_SCIENTIFIC_HASHES.json'
    hashes={r['path']:r['sha256'] for r in rows}
    if freeze.exists():
        for path,expected in json.loads(freeze.read_text()).items():
            assert hashes.get(path)==expected, 'Scientific input changed during handoff: '+path
    else:freeze.write_text(json.dumps(hashes,indent=2)+'\n')
    additions=[r for r in rows if r['decision']=='INCLUDE']
    (OUT/'SCIENTIFIC_SELECTED_PATHS.json').write_text(json.dumps([r['path'] for r in additions],indent=2)+'\n')
    (OUT/'EXCLUDED_RELEVANT_FILES.json').write_text(json.dumps([r for r in rows if r['decision']=='EXCLUDE'],indent=2)+'\n')
    print(json.dumps({'inventory_files':len(rows),'new_selected_files':len(additions),'new_selected_bytes':sum(r['size'] for r in additions),'excluded_relevant_files':sum(r['decision']=='EXCLUDE' for r in rows),'stages':STAGES,'figures':FIGURES},indent=2))
    print('Largest proposed new files:')
    for r in sorted(additions,key=lambda r:r['size'],reverse=True)[:15]:print(r['size'],r['path'])

if __name__=='__main__':main()
