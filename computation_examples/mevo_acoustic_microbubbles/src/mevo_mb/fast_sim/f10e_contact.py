"""Passive bounded observer around F10B tracing and F10C full-particle execution.

No force, contact, wall, event, or acceptance kernel is copied or changed.
Checking after each physical step permits stopping at the first original error.
The observer is enabled only in a registered window; no intermediate restart.
"""

from collections import deque
import gzip
import json
import numpy as np
from . import kernels as k
from .f10b_contact import ContactTrace, distances
from .f10c_engine import FullParticleBatch


class ContactAuditBatch(ContactTrace, FullParticleBatch):
    def __init__(self, *args, audit_start_s=0.0, numerical_mode="M0_DT10_PBD12", **kwargs):
        self.audit_enabled = False
        self.audit_start_s = audit_start_s
        self.numerical_mode = numerical_mode
        self.physical_steps = 0
        self.first_rejection = None
        self.traces = deque(maxlen=32)
        super().__init__(*args, **kwargs)
        if self.batch != 1:
            raise ValueError("Step diagnosis is registered only for one original member")

    def step(self, dt):
        if self.audit_enabled:
            return ContactTrace.step(self, dt)
        self.trace = None
        return FullParticleBatch.step(self, dt)

    def pairs(self):
        result = super().pairs()
        if self.trace is not None:
            self.trace.update(
                pair_force_N=self.pair_buffer.numpy().copy(),
                harmonics_complex_volume_m3=self.harmonics.numpy().copy(),
                pressure_complex_Pa=self.pressure.numpy().copy(),
                nearfield_cumulative_diag=self.near_diag.numpy().copy(),
            )
        return result

    def launch(self, kernel, dim, args):
        if self.trace is not None and kernel is k.accept:
            self.trace["error_bits_before_accept"] = int(self.error.numpy()[0])
        return super().launch(kernel, dim, args)

    def advance(self, dt, count, graphs=True):
        for _ in range(count):
            start = self.physical_steps * dt
            self.audit_enabled = start >= self.audit_start_s - 1e-14
            failure = None
            try:
                # Inherited arithmetic and original F6 contact-budget test.
                FullParticleBatch.advance(self, dt, 1, graphs=graphs and not self.audit_enabled)
            except RuntimeError as exc:
                failure = str(exc)
            q = self.quality()
            bits = int(self.error.numpy()[0])
            if bits:
                failure = "GPU_ERROR_" + str(bits)
            event = dict(
                step_index_zero_based=self.physical_steps,
                step_number_one_based=self.physical_steps + 1,
                start_s=start,
                attempted_end_s=start + dt,
                error_bits=bits,
                reason=failure,
            )
            if self.trace is not None:
                trace = self.trace
                trace.update(
                    event,
                    quality=q,
                    accepted=failure is None,
                    secondary_force_N=trace["total_force_N"] - trace["primary_force_N"],
                    prediction=distances(trace["predicted_positions_m"], self.p.radius_m),
                    final=distances(trace["proposal_positions_m"], self.p.radius_m),
                )
                trace["balanced_velocity_m_s"] = (
                    trace["predicted_positions_m"] - trace["initial_positions_m"]
                ) / dt
                trace["positions_after_step_m"] = self.x.numpy().copy()
                d = np.linalg.norm(
                    trace["proposal_positions_m"][:, None] - trace["proposal_positions_m"][None], axis=-1
                )
                ii, jj = np.where(np.triu(d < 1.99 * self.p.radius_m, 1))
                trace["residual_violating_pairs"] = np.c_[ii, jj]
                self.traces.append(trace)
            self.physical_steps += 1
            if failure:
                self.first_rejection = event
                raise RuntimeError(failure)

    def export_observer(self, out):
        from scripts.fast_pzt_f9.common import clean, write

        metadata = dict(
            numerical_mode=self.numerical_mode,
            dt_scope="Per-record actual dt",
            radius_m=self.p.radius_m,
            physics=self.p.metadata(),
            chi_zero=[self.chi_zero.real, self.chi_zero.imag],
            response_pressure_knots_Pa=self.response.knots,
            first_rejection_step=self.first_rejection,
            traced_steps=len(self.traces),
            attempted_steps=self.physical_steps,
            audit_start_s=self.audit_start_s,
            scope="Passive readbacks; original kernels. Last 32 audited steps retained. Diagnostic cost is not production performance.",
            original_quality=self.quality(),
            training_ready=False,
        )
        write(out / "OBSERVER.json", metadata)
        with gzip.open(out / "CONTACT_TRACE.json.gz", "wt") as stream:
            json.dump(
                clean(dict(metadata, traces=list(self.traces))),
                stream,
                separators=(",", ":"),
                allow_nan=False,
            )
