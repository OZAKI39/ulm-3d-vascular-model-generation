from pathlib import Path
import os, json, shutil, hashlib, subprocess, csv
from collections import Counter
from datetime import datetime, timezone
SRC=Path('/home/lzy/projects/ulm_3D_vascular')
DST=Path('/home/lzy/projects/.github-print-sync-20260925/ulm_3D_vascular')
META=DST/'reports/github_print_sync'
META.mkdir(parents=True,exist_ok=True)
LIMIT=50*1024**2
selected={}; excluded=[]; links=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''): h.update(b)
 return h.hexdigest()
def reject(p,reason):
 try: size=p.stat().st_size if p.is_file() else 0
 except OSError:size=0
 excluded.append({'path':str(p.relative_to(SRC)),'bytes':size,'reason':reason})
def add(path,dest=None):
 p=SRC/path
 if not p.exists() and not p.is_symlink():return
 rel=p.relative_to(SRC)
 if p.is_symlink(): links.append((p,dest or rel));return
 if any(x in rel.parts for x in ['.git','__pycache__','.pytest_cache','.venv','build-wsl','build','.ipynb_checkpoints']):return
 if p.name.endswith((':Zone.Identifier','.pyc','.pyo')) or p.name in ['.DS_Store']:return
 if p.is_dir():
  for child in sorted(p.iterdir()):add(child.relative_to(SRC))
 elif p.stat().st_size>LIMIT:reject(p,'single_file_over_50_MiB')
 else:selected[str(dest or rel)]=p
for name in ['vascular_processing','config','configs','tools','utils','test_data','docs/research_history','docs/ULTRALISER_PIPELINE.md',
             'requirements.txt','requirements-vascularmd.txt','requirements-fixture.txt','cfd_lumen_config.yaml','main.py',
             's1-1_swc_roi_generate_mouse.py','s1-2_swc_roi_generate_human.py','s1-3_swc_roi_generate_MeVO.py',
             's1-4_sacrificial_box_and_ports.py','s2_swc_stl_model_generate.py','TOPBRAIN_BRAVA_LABEL_TRANSFER.md']:
 add(name)
# Small shared utility modules retained to preserve import closure; no unrelated CFD results.
for p in (SRC/'tests').glob('*.py'):
 if p.name.startswith(('test_brava','test_compact','test_manufacturing','test_mevo','test_nn_production','test_roi_display',
                       'test_s1_4','test_swc','test_topbrain','test_vascularmd','test_boundary_review',
                       'test_ultraliser','test_rodent','test_sampling','test_schmid','test_mesh_cleanup',
                       'test_pipeline','test_hierarchical','test_voxelization','test_nne2')) or p.name in ['conftest.py','schmid_test_data.py','inspect_swc_rois.py']:
  add(p.relative_to(SRC))
for n in ['tests/data','tests/README_swc_roi_inspection.md','tests/swc_roi_inspection_output']:add(n)
for p in (SRC/'references').glob('*.md'):
 if any(t in p.name for t in ['BraVa','TopBrain','VascularMD']):add(p.relative_to(SRC))
for n in ['topbrain_brava_transfer','topbrain_mevo','topbrain_mevo_cache','topbrain_validation','mevo_validation',
          'vascularmd_validation','vascularmd_mouse_probe','mouse_brava_vmd_validation','mouse_vmd_branch_aic',
          'human_brava','human_brava_vmd','model_generate','reference_validation']:
 add('outputs/'+n)
add('vessel_model/T - Brava')
for p in (SRC/'vessel_model').glob('*vasculature_vmd*'):add(p.relative_to(SRC))
# Preserve the single raw mouse sample corresponding to the AIC processed sample, including paired images.
raw=SRC/'vessel_model/T - A high-resolution dataset of mouse brain vasculature'
for p in raw.rglob('fMOST_0_5_6_0_0_6_0001_02_01.*'):
 if 'BVLab-Annotation' not in p.parts:add(p.relative_to(SRC))
# Retain third-party source and license without nested git metadata, compiled binaries or example meshes/videos.
for p in (SRC/'third_party').iterdir():
 if p.is_file():add(p.relative_to(SRC))
for p in (SRC/'third_party/vascularmd').iterdir():
 if p.name not in ['.git','__pycache__','Output','Documentation']:add(p.relative_to(SRC))
for p in (SRC/'Ultraliser').iterdir():
 if p.name not in ['.git','build-wsl','build']:add(p.relative_to(SRC))
# Dataset provenance and official labels are small and sufficient to audit the cached exact masks.
source_data=SRC/'data/external/TopBrain'
add('data/external/TopBrain/dataset_provenance.json','data/topbrain_provenance/dataset_provenance.json')
for p in source_data.rglob('License.txt'):add(p.relative_to(SRC),'data/topbrain_provenance/License.txt')
for p in source_data.rglob('labelmap_topbrain_*.txt'):add(p.relative_to(SRC),Path('data/topbrain_provenance/itksnap_labelmap_txt')/p.name)
for p in source_data.rglob('labelsTr_topbrain_mr/*.nii.gz'):add(p.relative_to(SRC),Path('data/topbrain_provenance/labelsTr_topbrain_mr')/p.name)
# Latest device discovery is a diagnostic only, not a completed union or slice.
add('.codex_tmp/union_printer_probe.json','reports/github_print_sync/prior_union_printer_probe.json')
add('README.md','docs/SOURCE_WORKSPACE_README.md')
# Record omitted bulk files, including exact sizes. Never traverse symlinked datasets twice.
for rel,reason in [('data/external/TopBrain','bulk_original_TopBrain_archive_and_images'),
                   ('vessel_model/Schmid','bulk_dataset_not_required_by_BG001_print_fixture'),
                   ('third_party/vascularmd/Output','third_party_example_outputs'),
                   ('third_party/vascularmd/Documentation','third_party_demo_media')]:
 for p in (SRC/rel).rglob('*'):
  if p.is_file() and p not in selected.values():reject(p,reason)
for rel,reason in [('outputs/cfd_flow','unrelated_CFD_solver_runs'),('outputs/cfd_surface_prepare','unrelated_CFD_solver_runs'),
                   ('network_1d0d','unrelated_network_solver'),('Ultraliser/build-wsl','compiled_build_cache'),
                   ('.git','git_metadata_not_payload'),('.codex_tmp','scratch_not_public_evidence_except_printer_probe')]:
 excluded.append({'path':rel+'/', 'bytes':None,'reason':reason})
# Remove inherited files retired from this selected snapshot; only touches the isolated worktree.
base_paths=subprocess.check_output(['git','ls-files','-z'],cwd=DST).decode().split('\0')
removed=[]
for rel in base_paths:
 if rel and rel not in selected and rel not in ['.gitignore','.gitattributes','README.md']:
  p=DST/rel
  if p.is_file() or p.is_symlink():p.unlink();removed.append(rel)
entries=[]
for rel,p in sorted(selected.items()):
 target=DST/rel;target.parent.mkdir(parents=True,exist_ok=True)
 if target.is_symlink():target.unlink()
 shutil.copy2(p,target)
 h=sha(p)
 assert h==sha(target)
 entries.append({'path':rel,'source_path':str(p.relative_to(SRC)),'bytes':p.stat().st_size,'sha256':h,'kind':'file'})
for p,rel in links:
 try: target_rel=p.resolve(strict=True).relative_to(SRC)
 except (ValueError,OSError):reject(p,'external_or_missing_symlink_target');continue
 target=DST/rel; actual=DST/target_rel
 if not actual.exists():reject(p,'symlink_target_not_in_snapshot');continue
 target.parent.mkdir(parents=True,exist_ok=True)
 if target.is_symlink() or target.is_file():target.unlink()
 target.symlink_to(os.path.relpath(actual,target.parent),target_is_directory=actual.is_dir())
 entries.append({'path':str(rel),'source_path':str(p.relative_to(SRC)),'kind':'symlink','target':os.readlink(target),'source_target':os.readlink(p)})
# Ordinary Git blobs: all individual files are <=50 MiB, no external LFS data required.
(DST/'.gitattributes').write_text('* -filter\n*.py text eol=lf\n*.md text eol=lf\n*.yaml text eol=lf\n*.stl binary\n*.vtp binary\n*.vtk binary\n*.png binary\n*.3mf binary\n*.zip binary\n*.nii.gz binary\n')
ignore=(SRC/'.gitignore').read_text()
ignore+='\n# This branch carries the curated printing snapshot as ordinary Git files.\n.codex_tmp/\n.env\n.env.*\n*.pem\n*.key\n.venv/\nUltraliser/build*/\nthird_party/vascularmd/Output/\nthird_party/vascularmd/Documentation/\n'
(DST/'.gitignore').write_text(ignore)
third=[]
for rel in ['Ultraliser','third_party/vascularmd']:
 d={'path':rel,'upstream_commit':subprocess.check_output(['git','-C',str(SRC/rel),'rev-parse','HEAD'],text=True).strip(),
    'upstream_origin':subprocess.check_output(['git','-C',str(SRC/rel),'remote','get-url','origin'],text=True).strip(),
    'source_vendored':True,'nested_git_omitted':True}
 third.append(d)
(META/'source_file_manifest.json').write_text(json.dumps({'schema_version':1,'created_utc':datetime.now(timezone.utc).isoformat(),
 'source_root':str(SRC),'source_branch':subprocess.check_output(['git','branch','--show-current'],cwd=SRC,text=True).strip(),
 'source_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=SRC,text=True).strip(),
 'base_commit':'d8df343e129eb958e1860f5ece6f6b1da6e9f48d','branch':'sync/vascular-print-models-20260925',
 'maximum_file_bytes':LIMIT,'third_party':third,'entries':entries,'removed_retired_base_files':removed,
 'state':'SNAPSHOT_ONLY_O3_SMOOTHNESS_UNRESOLVED_UNION_NOT_IMPLEMENTED'},ensure_ascii=False,indent=2)+'\n')
with (META/'excluded_files.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','reason']);w.writeheader();w.writerows(excluded)
print(json.dumps({'copied':len(entries),'bytes':sum(x.get('bytes',0) for x in entries),'MiB':sum(x.get('bytes',0) for x in entries)/1048576,
 'excluded_records':len(excluded),'max_MiB':max(x.get('bytes',0) for x in entries)/1048576,'removed_retired_base_files':removed},ensure_ascii=False,indent=2))
