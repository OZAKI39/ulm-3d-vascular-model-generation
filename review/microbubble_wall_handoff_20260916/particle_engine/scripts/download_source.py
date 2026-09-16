from pathlib import Path
from datetime import datetime,timezone
import json,urllib.request,hashlib,tarfile,subprocess,shutil
W=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
R=Path('/workspace/microbubble_lammps/results/lammps_particle_engine_20260915_214759')
assert shutil.disk_usage('/workspace').free>=20*2**30
metadata=json.loads((W/'build_provenance/OFFICIAL_RELEASE_METADATA.json').read_text())
archive=W/'downloads'/metadata['source_asset_name']
assert not archive.exists() and not (W/'upstream/lammps').exists()
archive.parent.mkdir(parents=True,exist_ok=True)
started=datetime.now(timezone.utc).isoformat()
with urllib.request.urlopen(metadata['download_url'],timeout=60) as response,archive.open('wb') as f:
 shutil.copyfileobj(response,f,2**20)
h=hashlib.sha256()
with archive.open('rb') as f:
 for block in iter(lambda:f.read(2**20),b''):h.update(block)
assert h.hexdigest()==metadata['source_expected_sha256']
metadata.update(source_sha256=h.hexdigest(),download_started_utc=started,download_finished_utc=datetime.now(timezone.utc).isoformat(),download_bytes=archive.stat().st_size,status='PASS')
(W/'build_provenance/LAMMPS_SOURCE_PROVENANCE.json').write_text(json.dumps(metadata,indent=2)+'\n')
src=W/'upstream/lammps';src.mkdir(parents=True)
subprocess.run(['tar','-xzf',str(archive),'-C',str(src),'--strip-components=1'],check=True)
assert (src/'cmake/CMakeLists.txt').is_file()
for name in ['src/version.h','src/gitversion.h']:
 p=src/name
 if p.is_file():print(name,p.read_text())
print('SOURCE_DOWNLOAD_AND_SHA256=PASS',flush=True)
print('SOURCE_ROOT',src,flush=True)
