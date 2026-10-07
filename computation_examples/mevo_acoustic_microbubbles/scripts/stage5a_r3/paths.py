"""Read-only path translation for local postprocessing of server artifacts."""

from pathlib import Path
from scripts.stage5a_r3.common import ROOT, R2ROOT


def local(path):
    text = str(path)
    for marker, root in [
        ("/workspace/mevo_acoustic_microbubbles_stage5a_r3/", ROOT),
        ("/workspace/mevo_acoustic_microbubbles_stage5a_r2/", R2ROOT),
    ]:
        if text.startswith(marker):
            relative = text[len(marker) :]
            if relative.startswith(".execution/"):
                relative = relative.split("/", 2)[2]
            return root / relative
    return Path(path)
