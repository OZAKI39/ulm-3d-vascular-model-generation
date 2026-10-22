import subprocess,sys,json
from particle_3d.audit import sha256

def test_all_figures_reproducible_from_data(repo,report,tmp_path):
 subprocess.run([sys.executable,str(repo/'particle_3d/scripts/generate_particle65_report.py'),'--figures-only','--output-dir',str(tmp_path)],check=True)
 d=json.loads((report/'FIGURE_MANIFEST.json').read_text())
 assert len(d['figures'])==12
 for r in d['figures']:
  assert sha256(tmp_path/r['filename'])==r['sha256']
  assert all(sha256(report/p)==h for p,h in r['data_sha256'].items())
