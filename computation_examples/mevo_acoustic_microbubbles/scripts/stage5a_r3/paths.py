"""Read-only path translation for local postprocessing of server artifacts."""

from pathlib import Path
from scripts.stage5a_r3.common import ROOT, R2ROOT, BASE


def local(path):
    text = str(path)
    archive = "/workspace/mevo_acoustic_microbubbles_stage4l/results/stage4/lubrication_upgrade/dense_probe/"
    if text.startswith(archive):
        return BASE / "inputs/frozen_stage4l" / text[len(archive) :]
    for marker, root in [
        ("/workspace/mevo_acoustic_microbubbles_stage5a_r3/", ROOT),
        ("/workspace/mevo_acoustic_microbubbles_stage5a_r2/", R2ROOT),
        (
            "/home/lzy/projects/computation_examples/_archives/mevo_stage5a_r3_20261007/computation_examples/mevo_acoustic_microbubbles/",
            ROOT,
        ),
        (
            "/home/lzy/projects/computation_examples/_archives/mevo_stage5a_r2_20261006/computation_examples/mevo_acoustic_microbubbles/",
            R2ROOT,
        ),
    ]:
        if text.startswith(marker):
            relative = text[len(marker) :]
            if relative.startswith(".execution/"):
                relative = relative.split("/", 2)[2]
            return root / relative
    return Path(path)
