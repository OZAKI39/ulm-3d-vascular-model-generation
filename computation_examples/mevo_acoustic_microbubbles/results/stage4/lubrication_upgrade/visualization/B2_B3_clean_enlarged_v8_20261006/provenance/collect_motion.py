"""Collect independent case renders into one local viewing delivery."""

import argparse
import csv
import json
from pathlib import Path
import shutil

from scripts.stage4l.render_motion import ROOT, sha, viewer, write


def collect(out):
    assert not (out / "RENDER_MANIFEST.json").exists(), "Preserve completed deliveries."
    parts = sorted((out / "parts").glob("*/RENDER_MANIFEST.json"))
    assert len(parts) == 4, "All four independent renders must finish."
    videos, checks, states, protected, sources, provenance = {}, {}, {}, {}, {}, {}
    rows = []
    (out / "keyframes").mkdir()
    (out / "provenance").mkdir()
    first = None
    for p in parts:
        data = json.loads(p.read_text())
        assert data["status"] == "RENDERED_PENDING_DELIVERY_CHECK"
        assert not (set(videos) & set(data["videos"]))
        first = first or data
        for target, incoming in [(protected, data["protected_inputs"]), (sources, data["renderer_sources"])]:
            for key, digest in incoming.items():
                assert key not in target or target[key] == digest
                target[key] = digest
        for key, video in data["videos"].items():
            source = p.parent / video["path"]
            assert sha(source) == video["sha256"]
            shutil.copy2(source, out / video["path"])
            videos[key] = video
        for png in (p.parent / "keyframes").glob("*.png"):
            shutil.copy2(png, out / "keyframes" / png.name)
        for script in (p.parent / "provenance").glob("*"):
            target = out / "provenance" / script.name
            if target.exists():
                assert sha(target) == sha(script)
            else:
                shutil.copy2(script, target)
        checks.update(data["screen_checks"])
        states.update(json.loads((p.parent / "STATE_GROUPS.json").read_text()))
        with (p.parent / "frame_map.csv").open() as f:
            rows.extend(csv.DictReader(f))
        provenance[str(p.relative_to(out))] = sha(p)
    assert set(videos) == {f"{case}_{kind}" for case in ["B2", "B3"] for kind in ["dense_probe", "primary"]}
    # Show aggregation probes first in the viewer, then the complete primary transport.
    order = [f"{case}_{kind}" for kind in ["dense_probe", "primary"] for case in ["B2", "B3"]]
    videos = {key: videos[key] for key in order}
    rows.sort(key=lambda row: (order.index(row["video"]), int(row["frame"])))
    for p, digest in (protected | sources).items():
        assert sha(p) == digest, p
    sources[str(Path(__file__).resolve())] = sha(__file__)
    shutil.copy2(__file__, out / "provenance" / "collect_motion.py")
    write(out / "STATE_GROUPS.json", states)
    write(out / "SCREEN_CHECKS.json", checks)
    with (out / "frame_map.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    first.update(
        videos=videos,
        screen_checks=checks,
        protected_inputs=protected,
        renderer_sources=sources,
        case_identifiers_in_video=False,
        part_manifest_sha256=provenance,
        elapsed_s=None,
        render_execution="Four independent local render processes; see part manifests for timings.",
    )
    write(out / "RENDER_MANIFEST.json", first)
    viewer(out, videos)
    (out / "REPRODUCE.txt").write_text(
        f"cd {ROOT}\n"
        "PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python "
        "-m scripts.stage4l.render_motion --output <new-directory>\n"
        "PYTHONPATH=src:. /home/lzy/projects/computation_examples/_shared/python_environment/bin/python "
        "-m scripts.stage4l.check_motion <new-directory>\n"
        "Independent --kind/--case renders may instead be placed in <new-directory>/parts/<case>_<kind>, "
        "then combined with python -m scripts.stage4l.collect_motion <new-directory> before verification.\n"
    )
    print(json.dumps(dict(status="COLLECTED", videos=list(videos), output=str(out))), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    collect(parser.parse_args().output.resolve())
