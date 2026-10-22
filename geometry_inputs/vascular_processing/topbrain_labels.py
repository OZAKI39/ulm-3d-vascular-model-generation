"""Name-based ITK-SNAP label resolver with an explicit, sourced 2025 fallback."""
from dataclasses import dataclass
from pathlib import Path
import re
import shlex

from .topbrain_qc import TopBrainError, sha256

# Official common arterial labels, verified against the challenge data page.
# https://topbrain2025.grand-challenge.org/data/ (accessed 2026-09-24)
# This is NOT a complete MR/CT vein label map and is NOT a guessed TA36 map.
OFFICIAL_2025_ARTERIES = dict(enumerate([
    "BA", "R-P1P2", "L-P1P2", "R-ICA", "R-M1", "L-ICA", "L-M1",
    "R-Pcom", "L-Pcom", "Acom", "R-A1A2", "L-A1A2", "R-A3", "L-A3",
    "3rd-A2", "3rd-A3", "R-M2", "R-M3", "L-M2", "L-M3", "R-P3P4",
    "L-P3P4", "R-VA", "L-VA", "R-SCA", "L-SCA", "R-AICA", "L-AICA",
    "R-PICA", "L-PICA", "R-AChA", "L-AChA", "R-OA", "L-OA"], start=1))


def normalized(name):
    return re.sub(r"[\s_\-–—]+", "", name).casefold()


@dataclass
class LabelMap:
    values: dict[int, str]
    path: Path | None
    warnings: list[str]

    def resolve(self, name):
        matches = [value for value, label in self.values.items() if normalized(label) == normalized(name)]
        # TA36 subdivides the ICA; only its supraclinoid context is an alias here.
        if not matches and name in {"R-ICA", "L-ICA"}:
            matches = [v for v, label in self.values.items() if normalized(label) == normalized(name+"-C6-C7")]
        if len(matches) != 1:
            raise TopBrainError("LABEL_NOT_FOUND", f"Label {name!r} has {len(matches)} matches in {self.path}")
        return matches[0]

    def has(self, name):
        try:
            self.resolve(name)
            return True
        except TopBrainError:
            return False

    def report(self):
        return {"path": str(self.path) if self.path else None,
                "sha256": sha256(self.path) if self.path else None,
                "values": self.values, "warnings": self.warnings,
                "fallback_source": None if self.path else "https://topbrain2025.grand-challenge.org/data/"}


def parse_itksnap(path):
    path = Path(path)
    values, names = {}, set()
    for number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        tokens = shlex.split(line, comments=True)
        if not tokens:
            continue
        if len(tokens) < 8:
            raise TopBrainError("LABELMAP_INVALID", f"Expected ITK-SNAP ID RGB A VIS MESH NAME: {path}:{number}")
        try:
            value = int(tokens[0])
        except ValueError as exc:
            raise TopBrainError("LABELMAP_INVALID", f"Invalid ID at {path}:{number}") from exc
        name = " ".join(tokens[7:])
        if value < 0 or value in values or normalized(name) in names:
            raise TopBrainError("LABELMAP_INVALID", f"Duplicate/invalid label at {path}:{number}")
        values[value] = name
        names.add(normalized(name))
    if not values:
        raise TopBrainError("LABELMAP_INVALID", f"Empty label map: {path}")
    return LabelMap(values, path.resolve(), [])


def load_labelmap(root, modality, version, explicit=None):
    if explicit:
        return parse_itksnap(explicit)
    paths = sorted(p for directory in Path(root).rglob("itksnap_labelmap_txt") for p in directory.rglob("*.txt"))
    if paths:
        version_paths = [p for p in paths if bool(re.search(r"ta36|topaneu|v2", p.name, re.I)) == (version == "ta36")]
        if version_paths:
            paths = version_paths
        elif version == "ta36":
            paths = [p for p in paths if any("C6-C7" in name for name in parse_itksnap(p).values.values())]
            if not paths:
                raise TopBrainError("LABELMAP_NOT_FOUND", "No identifiable TA36 map; supply its actual --labelmap explicitly")
        modality_paths = [p for p in paths if re.search(rf"(?:^|[_\-.]){modality}(?:a|_|\.|-|$)", p.name, re.I)]
        if modality_paths:
            paths = modality_paths
        else:
            other = "ct" if modality == "mr" else "mr"
            paths = [p for p in paths if not re.search(rf"(?:^|[_\-.]){other}(?:a|_|\.|-|$)", p.name, re.I)]
        if len(paths) != 1:
            raise TopBrainError("LABELMAP_AMBIGUOUS", f"Use --labelmap to choose the correct map: {paths}")
        return parse_itksnap(paths[0])
    if version != "2025":
        raise TopBrainError("LABELMAP_NOT_FOUND", "TA36 requires its supplied label map; 2025 integers are not assumed")
    return LabelMap(OFFICIAL_2025_ARTERIES.copy(), None,
                    ["LABELMAP_FALLBACK: using verified 2025 common artery mapping; vein values are not interpreted"])
