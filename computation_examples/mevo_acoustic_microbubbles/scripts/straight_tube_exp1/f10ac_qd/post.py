"""Ordered CPU-only analysis and presentation; no automatic new physical jobs."""
import sys,time,subprocess
from .common import *
def main():
 stage=sys.argv[1];start=time.perf_counter()
 if stage=='delivery':
  from .report import main as report
  from .finalize import main as finalize
  report();finalize();return
 if stage=='components_report':
  from .group_diagnostics import main as groups
  from .report import main as report
  groups();report();return
 if stage=='review_media':
  from .figures import main as figures
  from .videos import main as videos
  from .resource_monitor import main as resources
  write(BASE/'evidence/FIRST_MEDIA_GENERATION.json',dict(figures=read(BASE/'figures/MANIFEST.json'),videos=read(BASE/'videos/MANIFEST.json')))
  resources();figures();videos('Video_F10AC_1_shape_diversity_discovery')
  subprocess.run(['.venv/bin/python','-m','scripts.straight_tube_exp1.f10ab.validate_media',str(BASE)],check=True)
  check=subprocess.run(['.venv/bin/python','-m','pytest','-q','tests/f10aa','tests/f10ab','tests/f10ac'],text=True,capture_output=True)
  write(BASE/'evidence/FINAL_TESTS.json',dict(returncode=check.returncode,stdout=check.stdout,stderr=check.stderr,identity=identity()))
  if check.returncode:raise RuntimeError(check.stdout+'\n'+check.stderr)
  write(BASE/'evidence/MEDIA_REVIEW_TIMING.json',dict(wall_s=time.perf_counter()-start,identity=identity()));return
 if stage=='analyze_media':
  from .audit import main as audit
  from .interpret import main as interpret
  from .figures import main as figures
  from .videos import main as videos
  audit();interpret();figures();videos()
  subprocess.run(['.venv/bin/python','-m','scripts.straight_tube_exp1.f10ab.validate_media',str(BASE)],check=True)
  write(BASE/'evidence/ANALYZE_MEDIA_TIMING.json',dict(wall_s=time.perf_counter()-start,identity=identity()));return
 if stage=='audit':
  from .audit import main as fn
 elif stage=='interpret':
  from .interpret import main as fn
 elif stage=='report':
  from .report import main as fn
 elif stage=='finalize':
  from .finalize import main as fn
 elif stage=='media':
  from .figures import main as figures
  from .videos import main as videos
  figures();videos()
  subprocess.run(['.venv/bin/python','-m','scripts.straight_tube_exp1.f10ab.validate_media',str(BASE)],check=True)
  write(BASE/'evidence/MEDIA_TIMING.json',dict(wall_s=time.perf_counter()-start,identity=identity()));return
 else:raise ValueError(stage)
 fn();write(BASE/'evidence'/(stage.upper()+'_TIMING.json'),dict(wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
