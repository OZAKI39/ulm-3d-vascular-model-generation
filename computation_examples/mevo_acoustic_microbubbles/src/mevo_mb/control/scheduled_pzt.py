"""Simulation-only ideal, phase-coherent four-disc schedules; no instrument I/O."""

from dataclasses import asdict, dataclass
import math

ORDER = ("Y-", "Y+", "X-", "X+")
MODEL = "HOMOGENEOUS_WATER_CONTROL_PROXY"
MODES = ("ONE_HOT_PZT", "AXIS_PAIR_CONTROL", "FOUR_CHANNEL_INDEPENDENT")
PLANNED_MAPPING = {"instrument_A": {"CH1": "X-", "CH2": "X+"},
                   "instrument_B": {"CH1": "Y-", "CH2": "Y+"}}


@dataclass(frozen=True)
class Segment:
    segment_id: int
    start_s: float
    duration_s: float
    mode: str
    amplitude_Pa_equivalent: tuple
    phase_rad: tuple
    enabled: tuple
    frequency_Hz: float = 2.5e6
    field_model: str = MODEL
    amplitude_unit: str = "Pa-equivalent"
    timing_model: str = "IDEAL_SCHEDULED_SEGMENTS"

    def __post_init__(self):
        for field in ("amplitude_Pa_equivalent", "phase_rad", "enabled"):
            object.__setattr__(self, field, tuple(getattr(self, field)))
        if not isinstance(self.segment_id, int) or self.segment_id < 0:
            raise ValueError("Nonnegative integer segment ID required")
        if not all(math.isfinite(v) for v in (self.start_s, self.duration_s)) or self.start_s < 0 or self.duration_s <= 0:
            raise ValueError("Finite nonnegative start and positive duration required")
        if self.mode not in MODES or self.frequency_Hz != 2.5e6 or self.field_model != MODEL:
            raise ValueError("Only registered modes and the existing 2.5 MHz water basis are admitted")
        if self.amplitude_unit != "Pa-equivalent" or self.timing_model != "IDEAL_SCHEDULED_SEGMENTS":
            raise ValueError("No voltage conversion or measured hardware timing is available")
        if any(len(v) != 4 for v in (self.amplitude_Pa_equivalent, self.phase_rad, self.enabled)):
            raise ValueError("Channel order must be Y-, Y+, X-, X+")
        if any(type(v) is not bool for v in self.enabled):
            raise ValueError("Boolean enables required")
        if not all(math.isfinite(v) and 0 <= v <= 1500 for v in self.amplitude_Pa_equivalent):
            raise ValueError("Source limit 0–1500 Pa-equivalent from F10Q; no clipping")
        if not all(math.isfinite(v) for v in self.phase_rad):
            raise ValueError("Finite phases required")
        if self.mode == "ONE_HOT_PZT" and sum(self.enabled) > 1:
            raise ValueError("ONE_HOT_PZT permits at most one enabled channel")
        if self.mode == "AXIS_PAIR_CONTROL":
            for i, j in [(0, 1), (2, 3)]:
                if self.enabled[i] != self.enabled[j] or self.amplitude_Pa_equivalent[i] != self.amplitude_Pa_equivalent[j]:
                    raise ValueError("Pair members share amplitude and enable; relative phase is free")

    @classmethod
    def from_legacy(cls, event, *, segment_id=0, start_s=0.0, duration_s=None, mode="FOUR_CHANNEL_INDEPENDENT"):
        if event.get("flow_enabled", False):
            raise ValueError("Strict zero background flow")
        return cls(segment_id, start_s, event["duration_s"] if duration_s is None else duration_s,
                   mode, tuple(event["amplitude_Pa"]), tuple(event["phase_rad"]), tuple(event["enabled"]),
                   frequency_Hz=event.get("frequency_Hz", 2.5e6), field_model=event.get("field_model", MODEL),
                   amplitude_unit=event.get("amplitude_unit", "Pa-equivalent"))

    def event(self):
        return dict(time_s=self.start_s, duration_s=self.duration_s,
                    amplitude_Pa=list(self.amplitude_Pa_equivalent), phase_rad=list(self.phase_rad),
                    enabled=list(self.enabled), flow_enabled=False)

    def record(self):
        return asdict(self)


def axis_pair(*, x_amplitude=0.0, y_amplitude=0.0, x_delta_deg=0.0, y_delta_deg=0.0,
              x_enabled=True, y_enabled=False, duration_s=0.002, start_s=0.0, segment_id=0):
    # delta = phase(plus) - phase(minus); both pair reference phases are zero.
    return Segment(segment_id, start_s, duration_s, "AXIS_PAIR_CONTROL",
                   (y_amplitude, y_amplitude, x_amplitude, x_amplitude),
                   (0.0, math.radians(y_delta_deg), 0.0, math.radians(x_delta_deg)),
                   (y_enabled, y_enabled, x_enabled, x_enabled))


def one_hot(channel, amplitude=1500.0, *, duration_s=0.002, start_s=0.0, segment_id=0):
    if channel not in ORDER:
        raise ValueError("Unknown channel")
    i = ORDER.index(channel)
    return Segment(segment_id, start_s, duration_s, "ONE_HOT_PZT",
                   tuple(amplitude if j == i else 0.0 for j in range(4)), (0.0,) * 4,
                   tuple(j == i for j in range(4)))


@dataclass(frozen=True)
class Schedule:
    segments: tuple

    def __post_init__(self):
        object.__setattr__(self, "segments", tuple(self.segments))
        if not self.segments or len(self.segments) > 8:
            raise ValueError("One to eight registered segments required")
        end = 0.0
        ids = set()
        for s in self.segments:
            if not isinstance(s, Segment) or abs(s.start_s - end) > 1e-15 or s.segment_id in ids:
                raise ValueError("Unique IDs and contiguous segments beginning at zero required")
            ids.add(s.segment_id)
            end = s.start_s + s.duration_s

    @property
    def duration_s(self):
        s = self.segments[-1]
        return s.start_s + s.duration_s

    def steps(self, dt_s):
        if not math.isfinite(dt_s) or dt_s <= 0:
            raise ValueError("Positive finite nominal time step required")
        for s in self.segments:
            count = int(math.floor(s.duration_s / dt_s + 1e-10))
            for k in range(count):
                yield s.start_s + k * dt_s, dt_s, s
            remainder = s.duration_s - count * dt_s
            if remainder > 1e-15:
                yield s.start_s + count * dt_s, remainder, s

    def record(self):
        return dict(channel_order=list(ORDER), segments=[s.record() for s in self.segments],
                    planned_mapping=PLANNED_MAPPING, hardware_executed=False,
                    calibration="ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED",
                    phase_reference="Re[P exp(+i omega t)]; common continuous simulation clock; phase offsets do not reset at switches",
                    timing="IDEAL_SCHEDULED_SEGMENTS; split the integration step exactly at each event boundary")

    @classmethod
    def from_record(cls, data):
        if data["channel_order"] != list(ORDER):
            raise ValueError("Historical channel order cannot change")
        return cls(tuple(Segment(**{k: tuple(v) if k in {"enabled", "phase_rad", "amplitude_Pa_equivalent"} else v
                                   for k, v in item.items()}) for item in data["segments"]))
