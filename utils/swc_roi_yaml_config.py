"""Strict YAML configuration loader for the SWC-to-ROI workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml

from .rodent_vasculature.config import RodentVasculatureConfig
from .sampling.sampling_config import SamplingConfig


@dataclass(frozen=True, slots=True)
class SWCROIRunConfig:
    """Validated settings needed by the command-line orchestration layer."""

    source_path: Path
    rodent: RodentVasculatureConfig
    sampling: SamplingConfig
    sampling_enabled: bool
    show_gui: bool
    verbose: bool


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a YAML mapping")
    return dict(value)


def _section(
    parent: Mapping[str, Any],
    key: str,
    *,
    allowed: set[str],
    required: set[str] | None = None,
) -> dict[str, Any]:
    if key not in parent:
        raise ValueError(f"Missing YAML section: {key}")
    values = _mapping(parent[key], key)
    unknown = sorted(set(values) - allowed)
    if unknown:
        raise ValueError(f"Unknown keys in {key}: {', '.join(unknown)}")
    missing = sorted((required or allowed) - set(values))
    if missing:
        raise ValueError(f"Missing keys in {key}: {', '.join(missing)}")
    return values


def _boolean(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be true or false")
    return value


def _integer(value: Any, label: str, *, nullable: bool = False) -> int | None:
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        suffix = " or null" if nullable else ""
        raise ValueError(f"{label} must be an integer{suffix}")
    return int(value)


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    return float(value)


def _string(value: Any, label: str, *, nullable: bool = False) -> str | None:
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not value.strip():
        suffix = " or null" if nullable else ""
        raise ValueError(f"{label} must be a non-empty string{suffix}")
    return value.strip()


def _sequence(value: Any, label: str, *, length: int | None = None) -> list[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError(f"{label} must be a YAML sequence")
    result = list(value)
    if length is not None and len(result) != length:
        raise ValueError(f"{label} must contain exactly {length} values")
    return result


def _float_tuple(value: Any, label: str, *, length: int | None = None) -> tuple[float, ...]:
    return tuple(
        _number(item, f"{label}[{index}]")
        for index, item in enumerate(_sequence(value, label, length=length))
    )


def _int_tuple(value: Any, label: str, *, length: int | None = None) -> tuple[int, ...]:
    values: list[int] = []
    for index, item in enumerate(_sequence(value, label, length=length)):
        converted = _integer(item, f"{label}[{index}]")
        assert converted is not None
        values.append(converted)
    return tuple(values)


def _resolve_path(value: Any, label: str, project_root: Path, *, nullable: bool = False) -> Path | None:
    text = _string(value, label, nullable=nullable)
    if text is None:
        return None
    path = Path(text).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def load_swc_roi_yaml_config(path: Path, *, project_root: Path) -> SWCROIRunConfig:
    """Load, type-check, and translate the human-facing YAML configuration."""

    source_path = Path(path).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"SWC/ROI YAML configuration does not exist: {source_path}")
    payload = yaml.safe_load(source_path.read_text(encoding="utf-8"))
    root = _mapping(payload, "YAML root")
    expected_root = {
        "schema_version",
        "paths",
        "pipeline",
        "coordinates",
        "preprocessing",
        "exports",
        "visualization",
        "sampling",
        "runtime",
    }
    unknown_root = sorted(set(root) - expected_root)
    missing_root = sorted(expected_root - set(root))
    if unknown_root:
        raise ValueError(f"Unknown top-level YAML keys: {', '.join(unknown_root)}")
    if missing_root:
        raise ValueError(f"Missing top-level YAML keys: {', '.join(missing_root)}")
    schema_version = _integer(root["schema_version"], "schema_version")
    if schema_version != 1:
        raise ValueError(f"Unsupported schema_version: {schema_version}; expected 1")

    paths = _section(
        root,
        "paths",
        allowed={"input_dir", "output_dir"},
    )
    pipeline = _section(
        root,
        "pipeline",
        allowed={
            "stage",
            "source_run",
            "cohort",
            "sample_id",
            "parent_group_id",
            "split",
            "max_samples",
        },
    )
    coordinates = _section(
        root,
        "coordinates",
        allowed={"spacing_xyz_um", "expected_shape_zyx", "enforce_expected_shape"},
    )
    preprocessing = _section(
        root,
        "preprocessing",
        allowed={
            "smooth_centerlines",
            "smoothing_window_points",
            "resample_step_um",
            "strict_nonpositive_radius",
            "analysis_component_id",
        },
    )
    exports = _section(
        root,
        "exports",
        allowed={"save_graphml", "save_vtp", "save_npz"},
    )
    visualization = _section(
        root,
        "visualization",
        allowed={"enabled", "max_samples", "max_direction_arrows", "figure2a"},
    )
    figure2a = _section(
        visualization,
        "figure2a",
        allowed={"enabled", "volume_opacity", "window_size", "show_gui"},
    )
    sampling = _section(
        root,
        "sampling",
        allowed={
            "enabled",
            "seed",
            "anchor",
            "roi",
            "features",
            "clustering",
            "selection",
            "comparison",
            "visualization",
        },
    )
    sampling_anchor = _section(
        sampling,
        "anchor",
        allowed={"mode", "min_distance_um", "max_candidates"},
    )
    sampling_roi = _section(
        sampling,
        "roi",
        allowed={"size_um", "min_branch_count", "max_cut_ports"},
    )
    sampling_features = _section(
        sampling,
        "features",
        allowed={
            "mode",
            "radius_quantiles",
            "scaler",
            "radius_weight",
            "structure_weight",
        },
    )
    sampling_clustering = _section(
        sampling,
        "clustering",
        allowed={"method", "n_clusters", "exploratory_k", "kmeans_max_iter"},
    )
    sampling_selection = _section(
        sampling,
        "selection",
        allowed={
            "mode",
            "target_count",
            "representatives_per_cluster",
            "max_overlap",
            "min_representative_distance_um",
        },
    )
    sampling_comparison = _section(
        sampling,
        "comparison",
        allowed={"compare_feature_modes"},
    )
    sampling_visualization = _section(
        sampling,
        "visualization",
        allowed={"max_roi_previews"},
    )
    runtime = _section(root, "runtime", allowed={"verbose"})

    project_root = Path(project_root).resolve()
    input_dir = _resolve_path(paths["input_dir"], "paths.input_dir", project_root)
    output_dir = _resolve_path(paths["output_dir"], "paths.output_dir", project_root)
    source_run = _resolve_path(
        pipeline["source_run"],
        "pipeline.source_run",
        project_root,
        nullable=True,
    )
    assert input_dir is not None and output_dir is not None

    enforce_shape = _boolean(
        coordinates["enforce_expected_shape"],
        "coordinates.enforce_expected_shape",
    )
    expected_shape = (
        _int_tuple(
            coordinates["expected_shape_zyx"],
            "coordinates.expected_shape_zyx",
            length=3,
        )
        if enforce_shape
        else None
    )
    visualizations_enabled = _boolean(
        visualization["enabled"], "visualization.enabled"
    )
    figure2a_enabled = visualizations_enabled and _boolean(
        figure2a["enabled"], "visualization.figure2a.enabled"
    )
    show_gui = _boolean(
        figure2a["show_gui"], "visualization.figure2a.show_gui"
    )
    sampling_enabled = _boolean(sampling["enabled"], "sampling.enabled")

    rodent = RodentVasculatureConfig(
        input_dir=input_dir,
        output_root=output_dir,
        stage=str(_string(pipeline["stage"], "pipeline.stage")),
        source_run=source_run,
        cohort=str(_string(pipeline["cohort"], "pipeline.cohort")),
        sample_id=_string(pipeline["sample_id"], "pipeline.sample_id", nullable=True),
        parent_group_id=_string(
            pipeline["parent_group_id"], "pipeline.parent_group_id", nullable=True
        ),
        split=_string(pipeline["split"], "pipeline.split", nullable=True),
        max_samples=_integer(
            pipeline["max_samples"], "pipeline.max_samples", nullable=True
        ),
        spacing_xyz_um=_float_tuple(
            coordinates["spacing_xyz_um"], "coordinates.spacing_xyz_um", length=3
        ),  # type: ignore[arg-type]
        expected_shape_zyx=expected_shape,  # type: ignore[arg-type]
        smoothing_enabled=_boolean(
            preprocessing["smooth_centerlines"], "preprocessing.smooth_centerlines"
        ),
        smoothing_window_points=int(
            _integer(
                preprocessing["smoothing_window_points"],
                "preprocessing.smoothing_window_points",
            )
        ),
        resample_step_um=_number(
            preprocessing["resample_step_um"], "preprocessing.resample_step_um"
        ),
        strict_nonpositive_radius=_boolean(
            preprocessing["strict_nonpositive_radius"],
            "preprocessing.strict_nonpositive_radius",
        ),
        analysis_component_id=_integer(
            preprocessing["analysis_component_id"],
            "preprocessing.analysis_component_id",
            nullable=True,
        ),
        visualizations_enabled=visualizations_enabled,
        max_visualization_samples=int(
            _integer(visualization["max_samples"], "visualization.max_samples")
        ),
        max_direction_arrows=int(
            _integer(
                visualization["max_direction_arrows"],
                "visualization.max_direction_arrows",
            )
        ),
        figure2a_enabled=figure2a_enabled,
        figure2a_volume_opacity=_number(
            figure2a["volume_opacity"], "visualization.figure2a.volume_opacity"
        ),
        figure2a_window_size=_int_tuple(
            figure2a["window_size"], "visualization.figure2a.window_size", length=2
        ),  # type: ignore[arg-type]
        save_graphml=_boolean(exports["save_graphml"], "exports.save_graphml"),
        save_vtp=_boolean(exports["save_vtp"], "exports.save_vtp"),
        save_npz=_boolean(exports["save_npz"], "exports.save_npz"),
    )

    sampling_config = SamplingConfig(
        output_root=output_dir,
        seed=int(_integer(sampling["seed"], "sampling.seed")),
        anchor_mode=str(_string(sampling_anchor["mode"], "sampling.anchor.mode")),
        min_anchor_distance_um=_number(
            sampling_anchor["min_distance_um"], "sampling.anchor.min_distance_um"
        ),
        max_candidate_anchors=int(
            _integer(
                sampling_anchor["max_candidates"], "sampling.anchor.max_candidates"
            )
        ),
        roi_size_um=_float_tuple(
            sampling_roi["size_um"], "sampling.roi.size_um", length=3
        ),  # type: ignore[arg-type]
        min_branch_count=int(
            _integer(
                sampling_roi["min_branch_count"], "sampling.roi.min_branch_count"
            )
        ),
        max_cut_ports=_integer(
            sampling_roi["max_cut_ports"],
            "sampling.roi.max_cut_ports",
            nullable=True,
        ),
        feature_mode=str(
            _string(sampling_features["mode"], "sampling.features.mode")
        ),
        radius_quantiles=_float_tuple(
            sampling_features["radius_quantiles"],
            "sampling.features.radius_quantiles",
        ),
        scaler=str(
            _string(sampling_features["scaler"], "sampling.features.scaler")
        ),
        radius_feature_weight=_number(
            sampling_features["radius_weight"], "sampling.features.radius_weight"
        ),
        structure_feature_weight=_number(
            sampling_features["structure_weight"],
            "sampling.features.structure_weight",
        ),
        clustering_method=str(
            _string(
                sampling_clustering["method"], "sampling.clustering.method"
            )
        ),
        n_clusters=int(
            _integer(
                sampling_clustering["n_clusters"],
                "sampling.clustering.n_clusters",
            )
        ),
        exploratory_k=_int_tuple(
            sampling_clustering["exploratory_k"],
            "sampling.clustering.exploratory_k",
        ),
        kmeans_max_iter=int(
            _integer(
                sampling_clustering["kmeans_max_iter"],
                "sampling.clustering.kmeans_max_iter",
            )
        ),
        selection_mode=str(
            _string(sampling_selection["mode"], "sampling.selection.mode")
        ),
        target_selected_count=int(
            _integer(
                sampling_selection["target_count"],
                "sampling.selection.target_count",
            )
        ),
        representatives_per_cluster=int(
            _integer(
                sampling_selection["representatives_per_cluster"],
                "sampling.selection.representatives_per_cluster",
            )
        ),
        max_selected_overlap=_number(
            sampling_selection["max_overlap"], "sampling.selection.max_overlap"
        ),
        min_representative_distance_um=_number(
            sampling_selection["min_representative_distance_um"],
            "sampling.selection.min_representative_distance_um",
        ),
        compare_feature_modes=_boolean(
            sampling_comparison["compare_feature_modes"],
            "sampling.comparison.compare_feature_modes",
        ),
        max_roi_previews=int(
            _integer(
                sampling_visualization["max_roi_previews"],
                "sampling.visualization.max_roi_previews",
            )
        ),
    )

    rodent.validate()
    sampling_config.validate()
    graph_stage = rodent.stage in {"all", "hierarchical-graph"}
    if graph_stage and figure2a_enabled and (show_gui or sampling_enabled) and not rodent.save_npz:
        raise ValueError(
            "exports.save_npz must be true when the saved Figure 2(a) view loads "
            "the graph for GUI or sampling-layer rendering"
        )
    return SWCROIRunConfig(
        source_path=source_path,
        rodent=rodent,
        sampling=sampling_config,
        sampling_enabled=sampling_enabled,
        show_gui=show_gui,
        verbose=_boolean(runtime["verbose"], "runtime.verbose"),
    )
