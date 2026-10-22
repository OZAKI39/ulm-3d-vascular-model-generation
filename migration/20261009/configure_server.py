"""Configure project-owned runtime and explicit historical path aliases."""
from pathlib import Path
import json,os,subprocess
root=Path('/workspace/vascular_printing');migration=root/'migration/20261009'
aliases={
 '/home/lzy/projects/vascular_printing':root,
 '/home/lzy/projects/ulm_3D_vascular':root,
 '/home/lzy/projects/ulm_particle_3d_particle0/ulm_3D_vascular':root,
 '/mnt/c/Users/lzy/AppData/Roaming/BambuStudio':migration/'external/bambu_profiles',
 '/mnt/c/Windows/Fonts':migration/'external/fonts',
 '/root/.local/share/fonts/vascular_printing':migration/'external/fonts',
}
for path,target in aliases.items():
 p=Path(path)
 if p.is_symlink() and p.resolve()==target.resolve():continue
 if p.exists() or p.is_symlink():raise RuntimeError('PATH_CONFLICT: '+path)
 p.parent.mkdir(parents=True,exist_ok=True);p.symlink_to(target,target_is_directory=True)
tools=root/'tools/bin';tools.mkdir(parents=True,exist_ok=True)
wrapper=tools/'bambu-studio'
wrapper.write_text('''#!/usr/bin/env bash
set -euo pipefail
PRINTING_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
export LC_ALL=C.UTF-8
exec xvfb-run -a "$PRINTING_ROOT/tools/bambu-studio-2.7.1.57/squashfs-root/AppRun" "$@"
''');wrapper.chmod(0o755)
(root/'config').mkdir(exist_ok=True)
(root/'config/server_runtime.env').write_text('''# Generated for the isolated vascular_printing deployment, 2026-10-09.
export VASCULAR_BAMBU_BIN=/workspace/vascular_printing/tools/bin/bambu-studio
export VASCULAR_BAMBU_NATIVE_BIN=/workspace/vascular_printing/tools/bambu-studio-2.7.1.57/squashfs-root/bin/bambu-studio
export VASCULAR_BAMBU_PROFILE_ROOT=/workspace/vascular_printing/migration/20261009/external/bambu_profiles
export PATH="/workspace/vascular_printing/tools/bin:$PATH"
export MPLBACKEND=Agg
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
''')
(migration/'path_aliases.json').write_text(json.dumps({k:str(v) for k,v in aliases.items()},indent=2))
subprocess.run(['fc-cache','-f',str(migration/'external/fonts')],check=True)
print('Project:',root,'Independent Python:',root/'.venv/bin/python')
print('Compatibility aliases:',len(aliases),'No copied credentials or Windows executables')
