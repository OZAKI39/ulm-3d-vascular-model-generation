"""Resolve archived asset paths inside the isolated deployment directory."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[1]
PREFIX=ROOT.parent
DEPENDENCIES=Path(os.environ.get("Y_CLOT_DEPENDENCY_ROOT", ROOT/"github_dependencies" if (ROOT/"github_dependencies").is_dir() else PREFIX/"deps"))
FONTS=Path(os.environ.get("Y_CLOT_FONT_DIR", DEPENDENCIES/"fonts"))
_MAP={
    '/home/lzy/projects/temp_storage/clot_clearance_github_sync_20260930/clot_clearance':ROOT,
    '/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles':DEPENDENCIES/'y_workflow',
    '/workspace/mevo_acoustic_microbubbles':DEPENDENCIES/'y_workflow',
    '/mnt/c/Windows/Fonts':FONTS,
}
def relocate(value):
    p=Path(value)
    for old,new in _MAP.items():
        try:
            target=new/p.relative_to(old)
            if old=="/mnt/c/Windows/Fonts" and not target.exists():
                substitute={"arial.ttf":"LiberationSans-Regular.ttf","arialbd.ttf":"LiberationSans-Bold.ttf","ariali.ttf":"LiberationSans-Italic.ttf","arialbi.ttf":"LiberationSans-BoldItalic.ttf"}.get(p.name)
                if substitute:
                    for folder in [Path("/usr/share/fonts/truetype/liberation2"),Path("/usr/share/fonts/truetype/liberation")]:
                        if (folder/substitute).exists():return folder/substitute
            return target
        except ValueError:pass
    return p
