"""Explicit operational validity limits; thresholds are policies, not universal laws."""

from dataclasses import dataclass, asdict

import numpy as np


class PhysicsValidityError(RuntimeError):
    pass


@dataclass(frozen=True)
class PairInteractionValidity:
    min_distance_over_radius: float = 10.0
    min_distance_over_sum: float = 5.0
    max_radial_fraction: float = 0.2
    max_kd: float = 0.3

    def check(self, positions, radii, radial_fraction, frequency, sound_speed):
        x, r = np.asarray(positions), np.asarray(radii)
        i, j = np.triu_indices(len(r), 1)
        d = np.linalg.norm(x[i] - x[j], axis=1)
        ri, rj = d / r[i], d / r[j]
        rs = d / (r[i] + r[j])
        amplitude = np.asarray(radial_fraction)
        near = (
            (ri < self.min_distance_over_radius)
            | (rj < self.min_distance_over_radius)
            | (rs < self.min_distance_over_sum)
        )
        high = (amplitude[i] > self.max_radial_fraction) | (amplitude[j] > self.max_radial_fraction)
        retardation = 2 * np.pi * frequency * d / sound_speed > self.max_kd
        valid = ~near & ~high & ~retardation
        return dict(
            pair_count=len(d),
            valid_pair_count=int(valid.sum()),
            far_field_fraction=float(valid.mean()) if len(d) else 1.0,
            near_field_pairs=int(near.sum()),
            amplitude_invalid_pairs=int(high.sum()),
            retardation_invalid_pairs=int(retardation.sum()),
            contact_count=int((rs <= 1).sum()),
            minimum_pair_separation_m=float(d.min()) if len(d) else None,
            min_d_over_Ri=float(ri.min()) if len(d) else None,
            min_d_over_Rj=float(rj.min()) if len(d) else None,
            min_d_over_sum=float(rs.min()) if len(d) else None,
        )


@dataclass(frozen=True)
class PairPhysicsRegime:
    """Separation warnings, informed by Stage 2.5; not empirical accuracy certificates."""
    far_field_over_sum: float = 5.0
    transition_over_sum: float = 2.0
    contact_relative_tolerance: float = 1e-5

    def __post_init__(self):
        if not 1 < self.transition_over_sum < self.far_field_over_sum:
            raise ValueError("Require 1 < transition < far-field separation")
        if not 0 <= self.contact_relative_tolerance < 0.01:
            raise ValueError("Invalid contact classification tolerance")

    def population(self, positions, radii):
        i, j = np.triu_indices(len(radii), 1)
        ratio = np.linalg.norm(positions[j]-positions[i], axis=1)/(radii[i]+radii[j])
        labels = np.full(len(i), "FAR_FIELD", dtype="U16")
        labels[ratio < self.far_field_over_sum] = "TRANSITION"
        labels[ratio < self.transition_over_sum] = "NEAR_CONTACT"
        labels[ratio <= 1+self.contact_relative_tolerance] = "CONTACT"
        result = dict(pair_count=len(i), evidence="Stage 2.5 separation evidence informs configurable warning thresholds; no accuracy guarantee")
        for key, label in [("far_field", "FAR_FIELD"), ("transition", "TRANSITION"), ("near_field", "NEAR_CONTACT"), ("contact", "CONTACT")]:
            result[key+"_fraction"] = float(np.mean(labels == label)) if len(i) else float(label == "FAR_FIELD")
        result["outside_validated_fraction"] = 0.0
        return result


class PhysicsValidityMonitor:
    def __init__(
        self,
        radial_model,
        *,
        policy="strict",
        pair=None,
        min_wall_distance_over_radius=5.0,
        regime_map_path=None,
        max_invalid_pair_fraction=0.0,
        max_transition_pair_fraction=0.05,
        regime=None,
    ):
        if policy not in ("strict", "legacy_demo", "pair_reference", "simulator"):
            raise ValueError("Unknown near_field_policy")
        self.radial_model, self.policy = radial_model, policy
        self.pair = pair or PairInteractionValidity()
        self.wall_ratio = min_wall_distance_over_radius
        self.history = []
        self.last = {}
        self.accepted_duration_s = 0.0
        self.invalid_duration_s = 0.0
        self.pair_time = 0.0
        self.valid_pair_time = 0.0
        self.response_context = {}
        self.regime_map = None
        if regime_map_path:
            from mevo_mb.validation.regime_maps import PairPhysicsRegimeMap
            self.regime_map = PairPhysicsRegimeMap(regime_map_path)
        self.simple_regime = PairPhysicsRegime(**(regime or {})) if policy == "simulator" else None
        self.invalid_limit = max_invalid_pair_fraction
        self.transition_limit = max_transition_pair_fraction
        if not 0 <= self.invalid_limit <= 1 or not 0 <= self.transition_limit <= 1:
            raise ValueError("Pair warning thresholds must be fractions")

    def record_accepted_interval(self, dt_s):
        self.accepted_duration_s += dt_s
        self.invalid_duration_s += dt_s * bool(self.last.get("invalid_reason"))
        self.pair_time += dt_s * self.last.get("pair_count", 0)
        self.valid_pair_time += dt_s * self.last.get("valid_pair_count", 0)

    def observe(
        self,
        positions,
        radii,
        time_s,
        *,
        radial_fraction,
        frequency,
        sound_speed,
        wall_distance_m=None,
        response_valid=True,
        harmonic_resolved=True,
        secondary_enabled=True,
        pressure_amplitude=None,
        pressure_phase_rad=None,
        shell_model=None,
    ):
        check = self.pair.check(positions, radii, radial_fraction, frequency, sound_speed)
        regime = dict(
            far_field_fraction=0.0,
            transition_fraction=0.0,
            near_field_fraction=0.0,
            contact_fraction=check["contact_count"] / max(check["pair_count"], 1),
            outside_validated_fraction=1.0 - check["contact_count"] / max(check["pair_count"], 1),
            pair_count=check["pair_count"],
            evidence="No matching empirical map supplied",
        )
        if self.regime_map is not None and pressure_amplitude is not None:
            phases = np.asarray(pressure_phase_rad if pressure_phase_rad is not None else 0.0)
            if np.max(np.abs(np.exp(1j * phases) - np.exp(1j * phases.flat[0]))) < 1e-10:
                regime = self.regime_map.population(
                    positions, radii, frequency, pressure_amplitude, shell_model, sound_speed
                )
        if self.simple_regime is not None:
            regime = self.simple_regime.population(np.asarray(positions), np.asarray(radii))
        if not check["pair_count"]:
            regime.update(far_field_fraction=1.0, outside_validated_fraction=0.0)
        invalid_fraction = (
            regime["near_field_fraction"] + regime["contact_fraction"] + regime["outside_validated_fraction"]
        )
        regime["invalid_pair_fraction"] = invalid_fraction
        warning = secondary_enabled and (
            invalid_fraction > self.invalid_limit or regime["transition_fraction"] > self.transition_limit
        )
        walls = (
            0
            if wall_distance_m is None
            else int(np.sum(np.asarray(wall_distance_m) / radii < self.wall_ratio))
        )
        reasons = []
        if not response_valid:
            reasons.append("response table bounds/cell/medium invalid")
        if not harmonic_resolved and secondary_enabled:
            reasons.append("secondary-force harmonic truncation unresolved")
        if secondary_enabled and check["near_field_pairs"]:
            reasons.append("pair separation outside far-field policy")
        if secondary_enabled and check["retardation_invalid_pairs"]:
            reasons.append("quasistatic acoustic retardation limit exceeded")
        if np.max(radial_fraction) > self.pair.max_radial_fraction:
            reasons.append("oscillation amplitude exceeds configured small-motion limit")
        if walls:
            reasons.append("near-wall hydrodynamics unvalidated")
        warning = warning or bool(reasons)
        self.last = dict(
            time_s=float(time_s),
            radial_model=self.radial_model,
            near_field_policy=self.policy,
            within_response_table=bool(response_valid),
            max_radial_fraction=float(np.max(radial_fraction)),
            wall_events=walls,
            invalid_reason=reasons,
            pair_regime=regime,
            validity_status="PHYSICS_VALIDITY_WARNING" if warning else "WITHIN_CONFIGURED_POLICY",
            **check,
        )
        self.history.append(self.last.copy())
        if reasons and self.policy == "strict":
            raise PhysicsValidityError("; ".join(reasons))
        return self.last

    def summary(self):
        counts = sum(x["pair_count"] for x in self.history)
        valid = sum(x["valid_pair_count"] for x in self.history)
        regime_keys = (
            "far_field_fraction",
            "transition_fraction",
            "near_field_fraction",
            "contact_fraction",
            "outside_validated_fraction",
            "invalid_pair_fraction",
        )
        fractions = {
            k: sum(x["pair_regime"][k] * x["pair_count"] for x in self.history) / counts
            if counts
            else (1.0 if k == "far_field_fraction" else 0.0)
            for k in regime_keys
        }
        result = dict(
            **self.last,
            response_context=self.response_context,
            policy_thresholds=asdict(self.pair),
            checks=len(self.history),
            sample_weighted_far_field_fraction=valid / counts if counts else 1.0,
            near_field_events=sum(x["near_field_pairs"] for x in self.history),
            total_wall_events=sum(x["wall_events"] for x in self.history),
            accepted_duration_s=self.accepted_duration_s,
            invalid_accepted_duration_s=self.invalid_duration_s,
            time_weighted_far_field_fraction=self.valid_pair_time / self.pair_time if self.pair_time else 1.0,
            physics_status="demo only"
            if self.policy == "legacy_demo"
            else (
                "outside validity"
                if any(x["invalid_reason"] for x in self.history)
                else "partially validated"
            ),
            note="Sample fraction includes RK stages; time fraction uses accepted interval endpoint flags. Passing policies is not experimental validation",
        )
        result["NEAR_CONTACT_PHYSICS_APPROXIMATED"] = any(x["pair_regime"]["near_field_fraction"]+x["pair_regime"]["contact_fraction"] > 0 for x in self.history)
        result["pair_regime"] = dict(fractions, weighting="evaluation-pair weighted, includes RK stages")
        result["validity_status"] = (
            "PHYSICS_VALIDITY_WARNING"
            if any(x["validity_status"] == "PHYSICS_VALIDITY_WARNING" for x in self.history)
            else "WITHIN_CONFIGURED_POLICY"
        )
        return result
