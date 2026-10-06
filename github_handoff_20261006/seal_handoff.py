"""Prepare reviewable manifests and a size-checked explicit staging allowlist.

Does not stage, commit, push, modify scientific files, or run simulations.
"""
from pathlib import Path
import collections, csv, hashlib, json, re, subprocess
from urllib.parse import unquote

R=Path(__file__).resolve().parents[1];O=R/'github_handoff_20261006'
ROOT_DOCS=['PROJECT_ACTIVE_CONTEXT.md','HISTORICAL_INDEX.md','HANDOFF_CONTEXT_20261006.md','HANDOFF_ENVIRONMENT.md','HANDOFF_SYNC_MANIFEST_20261006.md','LOCAL_ARTIFACTS_NOT_SYNCED.md']
SUPPORT=['SOURCE_GIT_STATE.json','ENVIRONMENT.json','EXTERNAL_DEPENDENCIES.md','SCIENTIFIC_FILE_INVENTORY.csv','SCIENTIFIC_SELECTED_PATHS.json','EXCLUDED_RELEVANT_FILES.json','PRE_HANDOFF_SCIENTIFIC_HASHES.json','LARGE_TABLE_AGGREGATES.json','PD2_POLYNOMIAL_FIGURE_DATA.csv','SERVER_TARGETED_TESTS.log','SERVER_TARGETED_TESTS.json','PORTABLE_TARGETED_TESTS.log','PORTABLE_TARGETED_TESTS.json','PRE_HANDOFF_ACTIVE_CONTEXT.md','TASK_REQUEST.txt','prepare_inventory.py','check_portable_snapshot.py','seal_handoff.py','STAGING_INVENTORY.csv','PROPOSED_STAGING_PATHS.txt','HANDOFF_QA.json']
paths=sorted(set(json.loads((O/'SCIENTIFIC_SELECTED_PATHS.json').read_text())+ROOT_DOCS+['github_handoff_20261006/'+x for x in SUPPORT]))
known=set(subprocess.check_output(['git','ls-files','-z'],cwd=R).decode().split('\0'))
excluded=json.loads((O/'EXCLUDED_RELEVANT_FILES.json').read_text())
state=json.loads((O/'SOURCE_GIT_STATE.json').read_text())
source=state['commands']['git branch --show-current'].strip();head=state['commands']['git rev-parse HEAD'].strip()
branch='clot-clearance/pd2-boundary-handoff-20261006'
assert subprocess.check_output(['git','branch','--show-current'],cwd=R,text=True).strip()==branch
assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()==head, 'Run sealing before commit only'
science_rows=list(csv.DictReader((O/'SCIENTIFIC_FILE_INVENTORY.csv').open()))
science={r['path']:r for r in science_rows}
changed_science=[p for p,h in json.loads((O/'PRE_HANDOFF_SCIENTIFIC_HASHES.json').read_text()).items() if hashlib.sha256((R/p).read_bytes()).hexdigest()!=h]
assert not changed_science,changed_science
assert not set(paths)&{r['path'] for r in excluded}
assert json.loads((O/'PORTABLE_TARGETED_TESTS.json').read_text())['PASS']
assert json.loads((O/'SERVER_TARGETED_TESTS.json').read_text())['exit_code']==0
figs=[p for p in paths if p.endswith('.png')];assert len(figs)==6
stages=sorted({p.split('/')[2] for p in paths if p.startswith('clot_clearance/verification/')})
for rel in paths:
    if not (R/rel).exists():
        assert rel in ['HANDOFF_SYNC_MANIFEST_20261006.md','github_handoff_20261006/STAGING_INVENTORY.csv','github_handoff_20261006/PROPOSED_STAGING_PATHS.txt','github_handoff_20261006/HANDOFF_QA.json']
        (R/rel).write_text('')
(O/'PROPOSED_STAGING_PATHS.txt').write_text('\n'.join(paths)+'\n')
allowed=set(paths)|known
broken=[]
for doc in ROOT_DOCS[:4]+['github_handoff_20261006/EXTERNAL_DEPENDENCIES.md']:
    # Historical index keeps historical links, whose absence is documented separately.
    if doc=='HISTORICAL_INDEX.md':continue
    for ref in re.findall(r'\]\(([^)]+)\)',(R/doc).read_text()):
        if re.match(r'\w+://',ref) or ref.startswith('#'):continue
        target=((R/doc).parent/unquote(ref.split('#')[0])).resolve()
        if not target.exists():broken.append({'file':doc,'target':ref})
        elif target.is_file() and target.relative_to(R).as_posix() not in allowed:broken.append({'file':doc,'not_selected':ref})
assert not broken,broken
secret_patterns=[rb'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----',rb'gh[pousr]_[A-Za-z0-9]{30,}',rb'github_pat_[A-Za-z0-9_]{40,}',rb'AKIA[0-9A-Z]{16}']
secret_hits=[]
for rel in paths:
    p=R/rel
    if p.suffix in {'.npz','.png'}:continue
    for pattern in secret_patterns:
        if re.search(pattern,p.read_bytes()):secret_hits.append(rel)
assert not secret_hits, 'Credential-pattern match in proposed staging files (contents not printed)'
tracked_updates=[p for p in paths if p in known];assert set(tracked_updates)=={'PROJECT_ACTIVE_CONTEXT.md','HISTORICAL_INDEX.md'}

def inventory():
    rows=[]
    for p in paths:
        reason=science[p]['reason'] if p in science else 'Handoff navigation, provenance, reproduction check, or explicit scope/size audit'
        rows.append(dict(path=p,size=(R/p).stat().st_size,type=Path(p).suffix,git_change='M' if p in known else 'A',reason_for_inclusion=reason))
    return rows

def manifest(rows):
    total=sum(x['size'] for x in rows);largest=max(rows,key=lambda x:x['size'])
    lines=['# PD2 boundary handoff sync manifest — 2026-10-06','',
      f'- Source branch: `{source}`',f'- Source HEAD: `{head}`',f'- New branch: `{branch}`',
      f"- Initial state recorded at: `{state['timestamp_utc']}`",
      '- Scientific status: **HIGHER_ORDER_BOUNDARY_MECHANICS=VERIFIED_FOR_INTACT_BENCHMARKS**.',
      '- Boundary reintegration: **PASSED**; four manufactured recoveries plus exactly one classical-shear solve were already complete before handoff.',
      f'- Files added: **{len(rows)-len(tracked_updates)}**; updated: **{len(tracked_updates)}**; total added/updated: **{len(rows)}**; deleted: **0**.',
      f'- Total new staged content: **{total} bytes ({total/2**20:.6f} MiB)**. This is the sum of complete added/updated blobs, including handoff metadata; Git pack compression and pre-existing ancestry are not counted.',
      f"- Largest synchronized file: `{largest['path']}`, **{largest['size']} bytes**.",
      '- Recommended next scientific task: **CURVED_CLOT_WALL_INTACT_MECHANICS_V1**, not started.',
      '', '## Scope and inventories','',
      'The source working tree included untracked completed research. This branch was created directly from its local HEAD and selectively adds that state. There was no reset to an old remote, cleanup, scientific edit, historical deletion, history rewrite, PR, or new scientific campaign.',
      '', '[Complete proposed staging paths](github_handoff_20261006/PROPOSED_STAGING_PATHS.txt) · [per-file staging size/type/reason](github_handoff_20261006/STAGING_INVENTORY.csv) · [scientific include/exclude inventory with SHA256](github_handoff_20261006/SCIENTIFIC_FILE_INVENTORY.csv). The staged Git index is audited independently before commit.',
      '', 'Source trees: `clot_clearance/pd_clot/pfpd_candidate/` including `higher_order_rk/` and `nonlocal_boundary/`; current related `scripts/`, `tests/`, and analysis source/data. Existing material, spectral split, geometry, regularization solver, configurations and the two-file `vendor/` helper remain inherited. No clone is copied.',
      '', 'Included verification stages:', '']
    lines += ['- `clot_clearance/verification/'+s+'/` — final reports/configuration/compact evidence; raw-file exclusions below.' for s in stages]
    lines += ['', '## Representative figures','']+['- `'+p+'`' for p in figs]
    lines += ['', 'Only these six PNGs are newly synchronized: two per current stage, all under 3 MiB; no new PDFs or redundant raster variants. Their small JSON/CSV data and source plotting scripts are included. The polynomial heatmap has an additional exact 36-row table in `github_handoff_20261006/PD2_POLYNOMIAL_FIGURE_DATA.csv`; extracting this table did not rerun the oracle.',
      '', '## Local-only artifacts and size policy','',
      f'**{len(excluded)} relevant files, {sum(x["size"] for x in excluded)} bytes**, remain local. Each path, exact byte size, SHA256, exclusion reason and compact replacement is listed in [LOCAL_ARTIFACTS_NOT_SYNCED.md](LOCAL_ARTIFACTS_NOT_SYNCED.md), with machine-readable records in [EXCLUDED_RELEVANT_FILES.json](github_handoff_20261006/EXCLUDED_RELEVANT_FILES.json). This linked per-file manifest is part of this sync manifest.',
      '', 'Excluded categories: large N19 reference caches and geometry; large/repeated node/bond field arrays; duplicate solve inputs and historical state copies with retained metrics; detailed large per-node tables with compact summaries; non-selected historical/duplicate/high-resolution figures; paper PDFs and extracted full texts; superseded implementation-check output. Unrelated raw data, environments, cloned repositories, compiler products, old recovery indexes and scratch caches were outside selection and remain untouched.',
      '', 'No staged file exceeds 5 MiB; no staged raw array exceeds 2 MiB; total content is below 50 MiB; no large-file exception is needed. Twenty-one compact NPZ files are selected only as evidence/tests: sixteen energy-force fixtures (including a 613108-byte N7 reference-operator fixture needed by unchanged tests) and five final boundary accepted states. The small N7 file is intentional fixed test data, not a disposable runtime cache; no N19 cache is included. Exact quantities are checked from the staging inventory.',
      '', 'The repository-wide NPZ ignore rule was inspected and left unchanged. Only the explicitly listed NPZ paths are force-added individually through the allowlist; this does not widen ignore rules or stage a directory blindly.',
      '', 'Original scientific manifests describe full server outputs and are preserved byte-for-byte. Excluded files may still be named in those manifests and linked in original reports. Their absence on GitHub is declared here, and does not change the original science or assertions. Existing historical tracked binaries/papers in ancestry were not removed.',
      '', '## Tests, integrity and limitations','',
      '- No active project writer was found before creating the branch; all five boundary solve outputs and final report were complete. Boundary FINAL_STATUS and FINAL_ARTIFACT_QA report success.',
      '- Existing focused suites on server: **36 passed in 4.42 s** (RK 12, energy-force 12, boundary 12). [Log](github_handoff_20261006/SERVER_TARGETED_TESTS.log). No full pytest or simulation campaign.',
      '- Isolated copy containing only selected/inherited required files: **20 passed** (11 energy-force, 9 boundary; 4 explicitly deselected). [Log](github_handoff_20261006/PORTABLE_TARGETED_TESTS.log). First copy omitted an already tracked vendor helper; adding that existing dependency to the copy resolved the import failure without editing tests or scientific code. No new branch dependency was missing.',
      '- Full 36-test replay and all original plots require listed local-only artifacts. The portable subset checks live small mechanics and saved-result gates; it does not rerun N19 equilibria or certify every excluded array.',
      '- All 1204 files frozen before handoff retain SHA256 identity, including omitted scientific artifacts; the expanded inventory additionally includes two inherited vendor files. No physics, production source, old scientific result, old test assertion, old manifest, or .gitignore is modified.',
      '- New handoff documentation local links resolve within the selected/ancestral branch files. Historical reports retain original local-only links. Source import inspection and the isolated focused tests check required local module coverage.',
      '- Credential-pattern scan found no matches; credentials and tokens were not printed. No environments, paper PDFs, clones or build artifacts are newly staged.',
      '', 'See [HANDOFF_CONTEXT_20261006.md](HANDOFF_CONTEXT_20261006.md) for scientific caveats and next-task scope; [HANDOFF_ENVIRONMENT.md](HANDOFF_ENVIRONMENT.md) for versions and commands. Current context was updated from accepted reviews; pre-handoff active context is archived for provenance and the historical index text remains intact.',
      '', '## Commit and remote verification procedure','',
      'Create one handoff commit, push only the new branch to origin without force, then compare local HEAD with `git ls-remote origin refs/heads/clot-clearance/pd2-boundary-handoff-20261006`. Recheck branch nonexistence immediately before push; stop if it has appeared. Compare the three protected remote refs with the recorded initial values. Post-push command evidence is kept locally under `github_handoff_20261006/` to avoid a second self-referential commit; the final response provides the verified SHA. This pre-commit manifest does not pretend the push had already occurred.',
      '', 'Deliverable URL: [handoff branch](https://github.com/OZAKI39/ulm-3d-vascular-model-generation/tree/clot-clearance/pd2-boundary-handoff-20261006). Stop after remote verification; no PR or further science.', '']
    return '\n'.join(lines)

# Explicit self-size fixed point: the manifest/inventory/QA sizes include themselves,
# without claiming a recursively self-containing SHA256.
for _ in range(20):
    before=[(p,(R/p).stat().st_size) for p in paths]
    rows=inventory()
    with (O/'STAGING_INVENTORY.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0].keys());w.writeheader();w.writerows(rows)
    (R/'HANDOFF_SYNC_MANIFEST_20261006.md').write_text(manifest(rows))
    qa={'PASS':True,'phase':'PRE_COMMIT','source_HEAD':head,'new_branch':branch,'added_updated_files':len(rows),'staged_content_bytes':sum(x['size'] for x in rows),'original_scientific_hashes_checked':1204,'original_scientific_files_changed':changed_science,'all_inventory_scientific_files':len(science_rows),'new_figure_count':len(figs),'new_NPZ_count':sum(p.endswith('.npz') for p in paths),'new_handoff_broken_links':broken,'credential_pattern_hits':secret_hits,'server_targeted_tests':'36 passed','portable_targeted_tests':'20 passed, 4 deselected','old_manifests_modified':False,'production_modified':False,'new_simulations':0,'tracked_deletions':0,'gitignore_modified':False,'remote_push_status':'Not yet performed; verify separately after commit'}
    (O/'HANDOFF_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
    if before==[(p,(R/p).stat().st_size) for p in paths]:break
else:raise RuntimeError('Self-size records did not stabilize')
rows=inventory();assert sum(x['size'] for x in rows)<=50*2**20
assert all(x['size']<=5*2**20 for x in rows)
assert all(x['size']<=2*2**20 for x in rows if x['type']=='.npz')
assert sum(p.endswith('.npz') for p in paths)==21
assert qa['staged_content_bytes']==sum(x['size'] for x in rows)
print(json.dumps(qa,indent=2))
print('TOP 30 PROPOSED FILES')
for row in sorted(rows,key=lambda x:x['size'],reverse=True)[:30]:print(row['size'],row['path'])

if '--print-paths' in __import__('sys').argv:
    print('COMPLETE PROPOSED STAGING SET (directory: explicit basenames)')
    groups=collections.defaultdict(list)
    for path in paths:groups[str(Path(path).parent)].append(Path(path).name)
    for directory,names in sorted(groups.items()):print(directory+': '+', '.join(names))
