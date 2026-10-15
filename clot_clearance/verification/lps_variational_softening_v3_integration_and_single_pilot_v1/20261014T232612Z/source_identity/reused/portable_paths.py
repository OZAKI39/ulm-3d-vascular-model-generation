"""Resolve archived asset paths inside the isolated deployment directory."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PREFIX=ROOT.parent
_MAP={
    '/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance':ROOT,
    '/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles':PREFIX/'deps/y_workflow',
    '/workspace/mevo_acoustic_microbubbles':PREFIX/'deps/y_workflow',
    '/mnt/c/Windows/Fonts':PREFIX/'deps/fonts',
}
def relocate(value):
    p=Path(value)
    for old,new in _MAP.items():
        try:return new/p.relative_to(old)
        except ValueError:pass
    return p
