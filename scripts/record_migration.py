"""One-time read-only comparison against original files; not a runtime dependency.

Original paths are taken only from the pre-copy provenance manifest. This audit
does not alter, move, remove, or import anything in those original workspaces.
"""
import ast
from collections import Counter
import csv
import difflib
import json
from pathlib import Path
from human_vessel_flow.paths import ROOT, read_json, write_json, sha256

def symbols(path):
    if path.suffix!='.py':return []
    tree=ast.parse(path.read_text())
    return [n.name for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef))]

def main():
    provenance=ROOT/'docs/provenance'
    manifest=read_json(provenance/'SOURCE_MANIFEST_BEFORE_COPY.json')+read_json(provenance/'DIAGNOSTIC_SOURCES_BEFORE_COPY.json')
    moves={r['original_destination']:r['destination'] for r in read_json(provenance/'LAYOUT_MOVES.json')}
    rows=[]; immutable=[]; source_checks=[];patches=[]
    by_dest={}
    for row in manifest:
        source=Path(row['source']);actual=sha256(source)
        source_checks.append(dict(source=str(source),sha256=actual,unchanged=actual==row['sha256']))
        relative=moves.get(row['destination'],row['destination']);dest=ROOT/relative
        copied=dict(row,destination=relative)
        if dest.is_file():
            copied.update(new_sha256=sha256(dest),byte_identical=sha256(dest)==row['sha256'])
            if row['category'] not in {'source','test','adapted_source'}:immutable.append(copied)
            else:by_dest[relative]=copied
        else:copied.update(status='PRUNED_UNUSED_CONTROLLER_IN_NEW_PROJECT_ONLY')
        rows.append(copied)
    # Refactored entry points originate in scripts rather than the copied modules.
    for name in ['NATIVE_RUNNER_REFACTOR','MESH_AUDIT_REFACTOR','POSTPROCESS_REFACTOR']:
        r=read_json(provenance/(name+'.json'))
        by_dest[r['destination']]=dict(source=r['source'],destination=r['destination'],category='source',sha256=r['source_sha256'])
    b=ROOT.parent/'computation_examples/brava_flow_roi_18mlmin/scripts'
    for old,new in [('prepare_geometry.py','src/human_vessel_flow/geometry.py'),
                    ('design_boundary.py','src/human_vessel_flow/boundary.py'),
                    ('prepare_case.py','src/human_vessel_flow/fem.py'),
                    ('generate_mesh.py','native/tetgen_mesh.py')]:
        by_dest[new]=dict(source=str(b/old),destination=new,category='refactored_entry',sha256=sha256(b/old))
    module_rows=[];function_rows=[]; adaptations=[]
    for relative,r in sorted(by_dest.items()):
        source=Path(r['source']);dest=ROOT/relative
        r=dict(r,new_sha256=sha256(dest));r['byte_identical']=r['sha256']==r['new_sha256']
        adaptations.append(r)
        old_names=symbols(source);new_names=symbols(dest)
        module_rows.append(dict(old_file=str(source),new_file=relative,
            status='BYTE_IDENTICAL' if r['byte_identical'] else 'REFACTORED_WITH_DIFF',
            original_sha256=r['sha256'],new_sha256=r['new_sha256'],
            old_symbols=';'.join(old_names),new_symbols=';'.join(new_names)))
        for name in new_names:
            function_rows.append(dict(old_file=str(source),old_symbol=name if name in old_names else '<top-level recipe / explicit wrapper>',
                                      new_file=relative,new_symbol=name,
                                      mapping='same symbol; inspect diff for imports/paths' if name in old_names else 'extracted or added explicit entry; inspect module diff'))
        if not r['byte_identical']:
            patches.extend(difflib.unified_diff(source.read_text().splitlines(True),dest.read_text().splitlines(True),fromfile=str(source),tofile=relative))
    # New orchestration modules have no pre-existing one-to-one file.
    mapped=set(by_dest)
    for p in sorted((ROOT/'src/human_vessel_flow').rglob('*.py')):
        relative=str(p.relative_to(ROOT))
        if relative in mapped:continue
        names=symbols(p)
        module_rows.append(dict(old_file='<new explicit orchestration>',new_file=relative,status='NEW_WRAPPER',
            original_sha256='',new_sha256=sha256(p),old_symbols='',new_symbols=';'.join(names)))
        function_rows.extend(dict(old_file='<new explicit orchestration>',old_symbol='',new_file=relative,new_symbol=n,mapping='new wrapper; numerical operations call retained modules') for n in names)
    for name,data in [('MODULE_MAPPING.csv',module_rows),('FUNCTION_CLASS_MAPPING.csv',function_rows)]:
        with (ROOT/'docs'/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    write_json(provenance/'FINAL_SOURCE_ADAPTATIONS.json',adaptations)
    (provenance/'FINAL_SOURCE_ADAPTATIONS.patch').write_text(''.join(patches))
    write_json(provenance/'FINAL_FILE_MAPPING.json',rows)
    write_json(provenance/'ORIGINAL_FILES_UNCHANGED.json',dict(status='PASS' if all(r['unchanged'] for r in source_checks) else 'FAIL',files=source_checks))
    write_json(ROOT/'validation/FINAL_IMMUTABLE_DATA.json',dict(status='PASS' if all(r['byte_identical'] for r in immutable) else 'FAIL',
                  files=len(immutable),bytes=sum(r['size'] for r in immutable),records=immutable))
    summary=dict(originals_unchanged=all(r['unchanged'] for r in source_checks),original_files_checked=len(source_checks),
                 immutable_files=len(immutable),immutable_bytes=sum(r['size'] for r in immutable),
                 migrated_module_rows=len(module_rows),mapped_symbols=len(function_rows),
                 module_status=dict(Counter(r['status'] for r in module_rows)),
                 copied_then_pruned_unused_controllers=sum(not (ROOT/moves.get(r['destination'],r['destination'])).exists() for r in manifest),
                 deleted_originals=0,new_formal_CFD_runs=0,new_formal_population_runs=0)
    write_json(provenance/'FINAL_MIGRATION_SUMMARY.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
