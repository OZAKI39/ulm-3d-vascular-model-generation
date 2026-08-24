"""Dataset discovery, pairing, grouping, and leakage auditing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


@dataclass(frozen=True, slots=True)
class RodentSampleRecord:
    sample_id: str
    cohort: str
    source_stem: str
    parent_group_id: str
    tile_index_zyx: tuple[int, int, int] | None
    split: str | None
    image_path: Path | None
    mask_path: Path | None
    swc_path: Path | None
    eligible: bool
    skip_reason: str | None

    def to_dict(self) -> dict[str, Any]:
        optional_auxiliary_notes = []
        if self.image_path is None:
            optional_auxiliary_notes.append("raw image TIFF unavailable; image display disabled")
        if self.mask_path is None:
            optional_auxiliary_notes.append("Mask TIFF unavailable; registration QC disabled")
        return {
            "sample_id": self.sample_id,
            "cohort": self.cohort,
            "source_stem": self.source_stem,
            "parent_group_id": self.parent_group_id,
            "tile_index_zyx": self.tile_index_zyx,
            "split": self.split,
            "image_path": str(self.image_path) if self.image_path else None,
            "mask_path": str(self.mask_path) if self.mask_path else None,
            "swc_path": str(self.swc_path) if self.swc_path else None,
            "image_size_bytes": self.image_path.stat().st_size if self.image_path else None,
            "mask_size_bytes": self.mask_path.stat().st_size if self.mask_path else None,
            "swc_size_bytes": self.swc_path.stat().st_size if self.swc_path else None,
            "eligible": self.eligible,
            "skip_reason": self.skip_reason,
            "optional_auxiliary_notes": optional_auxiliary_notes,
        }


@dataclass(frozen=True, slots=True)
class RodentCatalog:
    input_dir: Path
    records: tuple[RodentSampleRecord, ...]
    split_parent_overlaps: dict[str, list[str]]

    @property
    def eligible_records(self) -> tuple[RodentSampleRecord, ...]:
        return tuple(record for record in self.records if record.eligible)

    def report(self) -> dict[str, Any]:
        cohorts = sorted({record.cohort for record in self.records})
        return {
            "input_dir": str(self.input_dir),
            "record_count": len(self.records),
            "eligible_record_count": len(self.eligible_records),
            "skipped_record_count": len(self.records) - len(self.eligible_records),
            "cohorts": {
                cohort: {
                    "record_count": sum(record.cohort == cohort for record in self.records),
                    "eligible_count": sum(
                        record.cohort == cohort and record.eligible for record in self.records
                    ),
                    "parent_group_count": len(
                        {record.parent_group_id for record in self.records if record.cohort == cohort}
                    ),
                }
                for cohort in cohorts
            },
            "split_parent_overlaps": self.split_parent_overlaps,
            "important_interpretation": (
                "Eligibility is SWC-centric: a non-empty SWC is required, while raw image "
                "and Mask TIFFs are optional QC/display inputs. Author-provided train/"
                "validation/test labels are tile-level; parent-group overlap is reported "
                "and must be considered before downstream learning."
            ),
        }


def _raw_parent_and_tile(stem: str) -> tuple[str, tuple[int, int, int] | None]:
    parts = stem.split("_")
    if len(parts) >= 4:
        try:
            tile = tuple(int(value) for value in parts[-3:])
        except ValueError:
            tile = None
        else:
            return "_".join(parts[:-3]), tile
    return stem, None


def _train_parent_and_tile(stem: str) -> tuple[str, tuple[int, int, int] | None]:
    if stem.startswith("cut_"):
        return _raw_parent_and_tile(stem)
    return stem.split("_", 1)[0], None


def _split_index(train_root: Path) -> dict[str, str]:
    index: dict[str, str] = {}
    for split in ("train", "val", "test"):
        path = train_root / f"{split}.txt"
        if not path.is_file():
            continue
        for line in path.read_text(encoding="utf-8-sig").splitlines():
            name = line.strip()
            if name:
                index[Path(name).stem] = split
    return index


def _records_for_triplet(
    cohort: str,
    image_dir: Path,
    mask_dir: Path,
    swc_dir: Path,
    split_index: dict[str, str] | None = None,
) -> list[RodentSampleRecord]:
    images = {path.stem: path for path in image_dir.glob("*.tif")} if image_dir.is_dir() else {}
    masks = {path.stem: path for path in mask_dir.glob("*.tif")} if mask_dir.is_dir() else {}
    swcs = {path.stem: path for path in swc_dir.glob("*.swc")} if swc_dir.is_dir() else {}
    stems = sorted(set(images) | set(masks) | set(swcs))
    output: list[RodentSampleRecord] = []
    for stem in stems:
        image = images.get(stem)
        mask = masks.get(stem)
        swc = swcs.get(stem)
        reasons: list[str] = []
        if swc is None:
            reasons.append("missing SWC")
        elif swc.stat().st_size == 0:
            reasons.append("empty SWC")
        if cohort == "train":
            parent, tile = _train_parent_and_tile(stem)
        else:
            parent, tile = _raw_parent_and_tile(stem)
        output.append(
            RodentSampleRecord(
                sample_id=f"{cohort}__{stem}",
                cohort=cohort,
                source_stem=stem,
                parent_group_id=parent,
                tile_index_zyx=tile,
                split=(split_index or {}).get(stem),
                image_path=image,
                mask_path=mask,
                swc_path=swc,
                eligible=not reasons,
                skip_reason="; ".join(reasons) if reasons else None,
            )
        )
    return output


def _split_parent_overlaps(records: Iterable[RodentSampleRecord]) -> dict[str, list[str]]:
    groups = {
        split: {record.parent_group_id for record in records if record.split == split}
        for split in ("train", "val", "test")
    }
    return {
        "train_vs_val": sorted(groups["train"] & groups["val"]),
        "train_vs_test": sorted(groups["train"] & groups["test"]),
        "val_vs_test": sorted(groups["val"] & groups["test"]),
    }


def build_catalog(input_dir: Path, cohort: str = "all") -> RodentCatalog:
    input_dir = input_dir.resolve()
    definitions = {
        "raw-total": input_dir / "raw_data" / "total_vascular_data" / "total_vascular_data",
        "raw-analysis": input_dir / "raw_data" / "analysis_data" / "analysis_data",
        "train": input_dir / "train_data",
    }
    selected = definitions if cohort == "all" else {cohort: definitions[cohort]}
    split_index = _split_index(definitions["train"])
    records: list[RodentSampleRecord] = []
    for name, root in selected.items():
        records.extend(
            _records_for_triplet(
                name,
                root / "images",
                root / "mask",
                root / "swc",
                split_index if name == "train" else None,
            )
        )
    records.sort(key=lambda item: item.sample_id)
    return RodentCatalog(input_dir, tuple(records), _split_parent_overlaps(records))


def select_records(
    catalog: RodentCatalog,
    *,
    sample_id: str | None,
    parent_group_id: str | None,
    split: str | None,
    max_samples: int | None,
) -> list[RodentSampleRecord]:
    records = list(catalog.eligible_records)
    if sample_id:
        records = [
            record
            for record in records
            if record.sample_id == sample_id or record.source_stem == sample_id
        ]
    if parent_group_id:
        records = [record for record in records if record.parent_group_id == parent_group_id]
    if split:
        records = [record for record in records if record.split == split]
    if max_samples is not None:
        records = records[:max_samples]
    return records
