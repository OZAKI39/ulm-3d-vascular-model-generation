"""Final immutable execution identity and independent source/input verification."""
import sys,shlex,importlib.metadata as md
from .common import *
def main():
 verify();(BASE/'RESOURCE_STOP').touch();import mevo_mb.fast_sim.accelerated_water as water;import mevo_mb.fast_sim.accelerated_sequential as sqp
 source={str(p.relative_to(ROOT)):sha(p) for folder in ['scripts/straight_tube_exp1/f10ag','configs/straight_tube_exp1/f10ag','tests/f10ag'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in str(p)}
 server=BASE.resolve().parents[2]
 (server/'CURRENT_SOURCE_PATH.txt').write_text(str(ROOT)+'\n')
 (server/'activate.sh').write_text('source '+shlex.quote(str(server/'.venv/bin/activate'))+'\nexport PYTHONPATH='+shlex.quote(str(ROOT/'src')+':'+str(ROOT))+'\nexport MPLBACKEND=Agg\n')
 (server/'README.md').write_text('# F10AG execution handoff\n\nCurrent immutable source: `'+str(ROOT)+'`.\n\nActivate with `source '+str(server/'activate.sh')+'`. Individual scientific jobs retain their original source commits under `execution/`; do not edit those snapshots. Start the next study in a new Git checkout and result namespace.\n\nScientific entry point: `results/straight_tube_exp1/f10ag/F10AG_DEVELOPMENT_REPORT.md`. Large raw HF files and input hashes: `REMOTE_RAW_INDEX.json`.\n')
 write(BASE/'FINAL_EXECUTION_SEAL.json',dict(identity=identity(),execution_root=str(ROOT),environment=sys.prefix,python=sys.executable,source_sha256=source,imports=dict(water=water.__file__,sequential=sqp.__file__),versions={n:md.version(n) for n in ['numpy','scipy','warp-lang','h5py']},protected_inputs_pass=True,scope='Final source/environment verification; per-run scientific source commits remain authoritative. No new trajectory.'))
if __name__=='__main__':main()
