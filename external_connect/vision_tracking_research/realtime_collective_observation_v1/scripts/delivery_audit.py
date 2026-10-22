import sys,json,hashlib,subprocess,re,os
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT,dump,sha

def main():
 base=Path('/home/lzy/projects/computation_examples/_archives/mevo_vision_f10x_20261017');history=json.loads((ROOT/'evidence/history_before.json').read_text());bad=[]
 for rel,digest in history['unchanged_files'].items():
  if sha(base/rel)!=digest:bad.append(rel)
 camera=json.loads((ROOT/'evidence/camera_source_identity.json').read_text());cbad=[r['path'] for r in camera if sha(r['path'])!=r['sha256']]
 before_status=subprocess.check_output(['git','status','--porcelain'],cwd=base,text=True);assert not bad and not cbad and not before_status
 dump(ROOT/'evidence/history_after.json',dict(f10x_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=base,text=True).strip(),f10x_status=before_status,source_mismatches=bad,camera_mismatches=cbad,protected_history_unchanged=True))
 from pypdf import PdfReader
 from PIL import Image
 figures=sorted((ROOT/'figures').glob('*.png'));pdfs=sorted((ROOT/'figures').glob('*.pdf'))
 assert len(figures)>=6 and len(figures)==len(pdfs)
 for p in figures:Image.open(p).verify()
 for p in pdfs:assert len(PdfReader(p).pages)==1
 v=json.loads((ROOT/'evidence/video_validation.json').read_text());assert len(v)==3
 for a in v:
  assert sha(ROOT/a['path'])==a['sha256'];assert a['full_decode'] and a['pts_monotonic'] and a.get('max_pts_error_s',0)<1e-5
 gate=json.loads((ROOT/'FINAL_GATE.json').read_text());assert gate['LIVE_CAMERA_PASS'] is False and gate['TRAINING_READY'] is False and gate['SCIENTIFIC_CONTROL_ADMITTED'] is False
 assert gate['STABILITY_10_MIN_PASS'] and gate['VIDEOS_VALIDATED'] and gate['FIGURES_VALIDATED']
 # Delivery metadata excludes itself and an out-of-tree self-referential commit ID.
 excluded={'.venv','vendor','runtime','__pycache__','.pytest_cache'};rows=[]
 for p in ROOT.rglob('*'):
  if not p.is_file() or excluded.intersection(p.relative_to(ROOT).parts) or p.name in ['INPUT_OUTPUT_SHA256.json','DELIVERY_IDENTITY.json'] or p.suffix=='.tar':continue
  rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
 dump(ROOT/'INPUT_OUTPUT_SHA256.json',dict(scope='Local F10Y inputs and outputs; environments/runtime caches excluded; remote-only long states indexed separately',source_video={'path':'/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4','sha256':'3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c'},files=rows,remote_recovery_manifest='evidence/remote_recovery_manifest.json'))
 print(json.dumps({'files_hashed':len(rows),'figure_pairs':len(figures),'videos':len(v),'history_unchanged':True,'scientific_pass_claimed':False}))
if __name__=='__main__':main()
