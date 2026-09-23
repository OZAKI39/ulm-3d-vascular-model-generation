import json,sys,subprocess
from PIL import Image
from particle_3d.audit import sha256

def test_reproducible_eleven_figures_from_saved_data(repo,tmp_path):
    subprocess.run([sys.executable,'-B',str(repo/'particle_3d/scripts/generate_particle6_report.py'),'--output',str(tmp_path),'--no-notes'],cwd=repo,check=True)
    manifest=json.loads((repo/'particle_3d/reports/particle6/FIGURE_MANIFEST.json').read_text());assert len(manifest)==11
    for row in manifest:
        assert sha256(tmp_path/row['figure'])==row['sha256']
        for path,digest in row['sources'].items():assert sha256(repo/path)==digest
        with Image.open(tmp_path/row['figure']) as im:assert im.width>=1500 and im.height>=700
