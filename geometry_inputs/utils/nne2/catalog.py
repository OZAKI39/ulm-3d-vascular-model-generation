"""Read vdb.mat and keep only complete, locally resolvable NNE2 records."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.io import loadmat


@dataclass(frozen=True, slots=True)
class NNE2Record:
    record_id: int
    subject_id: str | None
    tree_id: int | None
    branching_order: int | None
    depth_um: float | None
    ref_pointer: str
    stack_pointer: str
    stack_index: int | None
    ref_um_per_px: float | None
    stack_um_per_px: float | None
    stack_step_um: float | None
    map_pointer: str
    reference_file: Path | None
    stack_dir: Path | None
    map_file: Path | None
    stack_name: str | None
    complete: bool
    skip_reasons: tuple[str, ...] = field(default_factory=tuple)

    @property
    def tree_key(self) -> str | None:
        if self.subject_id is None or self.tree_id is None:
            return None
        return f"{self.subject_id}_tree_{self.tree_id}"

    def row(self) -> dict[str, Any]:
        payload = asdict(self)
        for key in ("reference_file", "stack_dir", "map_file"):
            value = payload[key]
            payload[key] = str(value) if value else ""
        payload["skip_reasons"] = ";".join(self.skip_reasons)
        payload["tree_key"] = self.tree_key or ""
        return payload


@dataclass(frozen=True, slots=True)
class NNE2Catalog:
    records: tuple[NNE2Record, ...]

    @property
    def complete_records(self) -> tuple[NNE2Record, ...]:
        return tuple(item for item in self.records if item.complete)

    @property
    def skipped_records(self) -> tuple[NNE2Record, ...]:
        return tuple(item for item in self.records if not item.complete)

    def select(
        self,
        *,
        subject_id: str | None = None,
        tree_id: int | None = None,
        stack_name: str | None = None,
    ) -> tuple[NNE2Record, ...]:
        output = []
        for item in self.complete_records:
            if subject_id is not None and item.subject_id != subject_id:
                continue
            if tree_id is not None and item.tree_id != tree_id:
                continue
            if stack_name is not None and item.stack_name != stack_name:
                continue
            output.append(item)
        return tuple(output)

    def report(self) -> dict[str, Any]:
        reasons: dict[str, int] = {}
        for item in self.skipped_records:
            for reason in item.skip_reasons:
                reasons[reason] = reasons.get(reason, 0) + 1
        return {
            "total_record_count": len(self.records),
            "complete_record_count": len(self.complete_records),
            "skipped_record_count": len(self.skipped_records),
            "complete_tree_count": len(
                {item.tree_key for item in self.complete_records if item.tree_key}
            ),
            "complete_stack_count": len(
                {item.stack_name for item in self.complete_records if item.stack_name}
            ),
            "skip_reason_counts": reasons,
            "missing_data_policy": "skip_record_without_processing",
        }


def _scalar(value: Any) -> Any:
    while isinstance(value, np.ndarray) and value.size == 1:
        value = value.reshape(-1)[0]
    if isinstance(value, np.ndarray) and value.dtype.kind in "US":
        return "".join(str(part) for part in value.reshape(-1))
    return value


def _text(value: Any) -> str:
    value = _scalar(value)
    if value is None:
        return ""
    if isinstance(value, np.ndarray) and value.size == 0:
        return ""
    return str(value).strip()


def _number(value: Any, kind: type[int] | type[float]) -> int | float | None:
    value = _scalar(value)
    if isinstance(value, np.ndarray) and value.size == 0:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(result):
        return None
    return int(result) if kind is int else float(result)


def _payload_files(input_dir: Path, folder_name: str, suffixes: Iterable[str]) -> list[Path]:
    root = input_dir / folder_name
    suffixes_lower = {value.lower() for value in suffixes}
    if not root.is_dir():
        return []
    return sorted(
        item
        for item in root.rglob("*")
        if item.is_file() and item.suffix.lower() in suffixes_lower
    )


def _path_leaf(pointer: str) -> str:
    return pointer.rstrip("/\\").replace("\\", "/").rsplit("/", 1)[-1]


def _resolve_reference(pointer: str, references: list[Path]) -> Path | None:
    leaf = _path_leaf(pointer)
    if not leaf or leaf == "0":
        return None
    matches = [item for item in references if item.name.startswith(leaf)]
    return matches[0] if len(matches) == 1 else None


def _resolve_stack(pointer: str, stack_index: dict[str, Path]) -> Path | None:
    leaf = _path_leaf(pointer)
    if not leaf or leaf == "0":
        return None
    return stack_index.get(leaf)


def load_nne2_catalog(input_dir: Path) -> NNE2Catalog:
    input_dir = Path(input_dir).expanduser().resolve()
    matrix = loadmat(input_dir / "vdb.mat", squeeze_me=False, struct_as_record=False).get(
        "vdb_gnu"
    )
    if matrix is None or matrix.ndim != 2 or matrix.shape[1] < 24:
        raise ValueError("vdb.mat does not contain the expected vdb_gnu 24-column table")

    references = _payload_files(input_dir, "hana_refs", {".tif", ".tiff"})
    stack_root = input_dir / "hana_stk"
    stack_dirs = [
        item
        for item in stack_root.rglob("*")
        if item.is_dir() and any(item.glob("*Ch2*.tif"))
    ] if stack_root.is_dir() else []
    stacks_by_name = {item.name: item for item in stack_dirs}
    maps = _payload_files(input_dir, "maps", {".jpg", ".jpeg", ".png", ".bmp"})
    maps_by_subject = {item.stem: item for item in maps}

    records: list[NNE2Record] = []
    for row_index in range(1, matrix.shape[0]):
        row = matrix[row_index]
        subject = _text(row[0]) or None
        tree_id = _number(row[1], int)
        branching_order = _number(row[5], int)
        depth_um = _number(row[8], float)
        ref_pointer = _text(row[15])
        stack_pointer = _text(row[16])
        stack_frame = _number(row[18], int)
        ref_scale = _number(row[20], float)
        stack_scale = _number(row[21], float)
        map_pointer = _text(row[22])
        stack_step = _number(row[23], float)
        reference_file = _resolve_reference(ref_pointer, references)
        stack_dir = _resolve_stack(stack_pointer, stacks_by_name)
        map_file = maps_by_subject.get(subject or "")
        reasons: list[str] = []
        required = {
            "missing_subject_id": subject,
            "missing_tree_id": tree_id,
            "missing_branching_order": branching_order,
            "missing_depth": depth_um,
            "missing_reference_image": reference_file,
            "missing_stack": stack_dir,
            "missing_stack_index": stack_frame,
            "missing_reference_scale": ref_scale,
            "missing_stack_scale": stack_scale,
            "missing_stack_step": stack_step,
            "missing_subject_map": map_file,
        }
        for reason, value in required.items():
            if value is None or value == "":
                reasons.append(reason)
        if stack_frame is not None and stack_frame < 1:
            reasons.append("invalid_stack_index")
        if ref_scale is not None and ref_scale <= 0:
            reasons.append("invalid_reference_scale")
        if stack_scale is not None and stack_scale <= 0:
            reasons.append("invalid_stack_scale")
        if stack_step is not None and stack_step <= 0:
            reasons.append("invalid_stack_step")
        records.append(
            NNE2Record(
                record_id=row_index,
                subject_id=subject,
                tree_id=int(tree_id) if tree_id is not None else None,
                branching_order=int(branching_order) if branching_order is not None else None,
                depth_um=depth_um,
                ref_pointer=ref_pointer,
                stack_pointer=stack_pointer,
                stack_index=int(stack_frame) if stack_frame is not None else None,
                ref_um_per_px=ref_scale,
                stack_um_per_px=stack_scale,
                stack_step_um=stack_step,
                map_pointer=map_pointer,
                reference_file=reference_file,
                stack_dir=stack_dir,
                map_file=map_file,
                stack_name=stack_dir.name if stack_dir else None,
                complete=not reasons,
                skip_reasons=tuple(sorted(set(reasons))),
            )
        )
    return NNE2Catalog(tuple(records))
