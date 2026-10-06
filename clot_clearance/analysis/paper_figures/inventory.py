"""Read-only metadata inventory of project evidence; never traverses dependencies."""
import csv,json,os
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
FORMATS={'.csv','.tsv','.json','.jsonl','.yaml','.yml','.pkl','.pickle','.npz','.npy','.h5','.hdf5','.parquet','.md','.log','.txt','.png','.pdf','.svg','.py','.ipynb'}
PRUNE={'.git','__pycache__','.venv','node_modules','external','julia_depot','julia_env','site-packages','dist-packages'}
CURRENT=('pfpd_baqp_finite_domain_boundary_consistency_v1','pfpd_at1_activation_root_cause_and_fix_v1','pfpd_at1_c5_stability_v1','pfpd_at1_c7_early_branch_v1','pfpd_at1_c5_activation_stationarity_audit_v1','pfpd_at1_petsc_vi_solver_v1')
def classify(path):
 s=str(path)
 if any(w in s for w in ('cleanup','recovery','MANIFEST','PRESERVATION','DELETE_')):return 'administrative','No scientific observable; provenance/administration only','no'
 if any(w in path.parts for w in ('tests','test_artifacts','fixtures','source_snapshot')):return 'legacy_or_fixture','Test/fixture or source snapshot; not accepted-state evidence','no'
 if any(w in s for w in CURRENT):
  quantity=('Force layers, continuum traction and decomposition' if 'boundary_consistency' in s else 'Common-strain localization and loading controls' if 'root_cause_and_fix' in s else 'Accepted states, activation, stationarity and stability')
  return 'latest_for_respective_scope',quantity,'candidate_requires_field_validation'
 if 'pfpd_' in s:return 'legacy_or_model_ambiguous','Earlier candidate/discretization; not interchangeable with current model','review_only'
 return 'legacy_or_unreviewed','Earlier clot/CFD/erosion output or supporting material','review_only'
def main():
 rows=[];skipped=[]
 for top in ('verification','results','analysis','output','figures','reports','diagnostics','audit','runs','scripts','notebooks','visualization'):
  start=ROOT/top
  if not start.exists():continue
  for folder,dirs,files in os.walk(start,followlinks=False):
   for d in list(dirs):
    p=Path(folder)/d
    if d in PRUNE or p==HERE or p.is_symlink():dirs.remove(d);skipped.append(str(p.relative_to(ROOT)))
   for name in files:
    p=Path(folder)/name
    if p.suffix.lower() not in FORMATS or not p.is_file():continue
    rel=p.relative_to(ROOT);status,quantity,usable=classify(rel)
    if p.suffix=='.py' and 'plot' not in name and 'figure' not in name:continue
    rows.append(dict(path=rel.as_posix(),data_type=p.suffix.lstrip('.'),bytes=p.stat().st_size,quantity=quantity,directly_usable=usable,scientific_status=status,inspection='metadata only; selected inputs validated by plotting loader'))
 out=HERE/'output/manifests';out.mkdir(parents=True,exist_ok=True)
 with (out/'data_inventory.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(sorted(rows,key=lambda r:r['path']))
 counts=Counter(r['scientific_status'] for r in rows)
 (out/'inventory_scope.json').write_text(json.dumps({'scope':'Available first-party scientific result/report/log/figure files, including legacy/cache/test outputs; metadata only. Dependency trees and symlink recursion excluded. No scientific calculation or repository-wide hashing.','files':len(rows),'by_status':dict(counts),'excluded_directories':skipped},indent=2)+'\n')
 print('Inventory:',len(rows),'files;',dict(counts))
if __name__=='__main__':main()
