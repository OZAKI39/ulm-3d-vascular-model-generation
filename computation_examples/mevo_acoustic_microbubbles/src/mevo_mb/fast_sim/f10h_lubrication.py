"""F10H explicit lubrication-count diagnostic; inherited physical kernels unchanged.

Graph replay is one physical step followed by the original host quality check.
The same 178nm pass budget is checked before acceptance, also in CUDA graphs.
Only B=1 is admitted until a separate batch gate passes.
"""

from collections import deque
from dataclasses import asdict, replace
import gzip
import json
import numpy as np
import warp as wp
from . import kernels as k
from .engine import FastSimulator
from .nearfield import F2Physics
from .f10b_contact import ContactTrace
from .f10c_engine import FullParticleBatch

MODES = {"M0_LUB2_PBD12_DT10": (1e-5, 2), "M0_LUB4_PBD12_DT10": (1e-5, 4), "M0_LUB4_PBD12_DT5": (5e-6, 4)}


@wp.kernel
def preaccept_pass_budget(diag: wp.array2d(dtype=wp.float64), error: wp.array(dtype=int), radius: wp.float64):
    i = wp.tid()
    if diag[i, 6] > wp.float64(0.1) * radius:
        wp.atomic_or(error, 0, 2048)


@wp.kernel
def accepted_clock(clock: wp.array(dtype=int), error: wp.array(dtype=int)):
    if error[0] == 0:
        clock[0] = clock[0] + 1


class LubricationDiagnostic(ContactTrace, FullParticleBatch):
    mode_registry = MODES
    registered_physical_steps = (1e-5, 5e-6, 2.5e-6)

    def __init__(self, *args, numerical_mode="M0_LUB2_PBD12_DT10", diagnostic=True, **kwargs):
        if numerical_mode not in self.mode_registry:
            raise ValueError("Only explicitly registered F10H modes")
        self.numerical_mode = numerical_mode
        self.diagnostic = bool(diagnostic)
        self.audit_enabled = False
        self.physical_steps = 0
        self.accepted_steps = 0
        self.first_rejection = None
        self.traces = deque(maxlen=32)
        self.ledger = []
        self.all_steps = []
        self.executed_controls = []
        self.control_plan = []
        self.control_index = 0
        kwargs["graph_chunk"] = 1
        super().__init__(*args, **kwargs)
        if self.batch != 1:
            raise ValueError(
                "F10H batch not admitted; one original member only; no batch production admission"
            )
        if asdict(self.p) != asdict(F2Physics()):
            raise ValueError(
                "Supply unchanged F2 physics; opt-in iteration count is applied after allocation"
            )
        self.p = replace(self.p, mobility_iterations=self.mode_registry[numerical_mode][1])
        self.final_proposal = self.tmp1
        self.previous_diag = self.diag.numpy().copy()
        self.previous_wall_stats = self.wall_stats.numpy().copy()
        wp.load_module(module=preaccept_pass_budget.module, device=self.device)

    def configure_schedule(self, events, dt):
        if dt != self.mode_registry[self.numerical_mode][0]:
            raise ValueError("Mode and dt mismatch")
        if not events or events[0]["time_s"] != 0:
            raise ValueError("Control must start at zero")
        ticks = [round(e["time_s"] / dt) for e in events]
        if any(abs(i * dt - e["time_s"]) > 1e-14 for i, e in zip(ticks, events)) or any(
            b <= a for a, b in zip(ticks, ticks[1:])
        ):
            raise ValueError("Controls must align exactly to this physical substep grid")
        self.control_plan = list(zip(ticks, events))

    def control(self, event):
        record = super().control(event)
        if hasattr(self, "executed_controls"):
            self.executed_controls.append(record)
        return record

    def graph(self, dt, count):
        if count != 1:
            raise ValueError("F10H graph must stop at every physical substep")
        old = self.audit_enabled
        self.audit_enabled = False
        try:
            return super().graph(dt, count)
        finally:
            self.audit_enabled = old

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
            )
        return result

    def launch(self, kernel, dim, args):
        if kernel is k.accept:
            self.final_proposal = args[0]  # correct for odd AND even round counts
            FastSimulator.launch(
                self, preaccept_pass_budget, self.total, [self.diag, self.error, self.p.radius_m]
            )
        if kernel is k.advance_clock:
            return FastSimulator.launch(self, accepted_clock, 1, [self.step_count, self.error])
        return super().launch(kernel, dim, args)

    def advance(self, dt, count, graphs=True):
        if dt != self.mode_registry[self.numerical_mode][0] or count < 1:
            raise ValueError("Registered dt and positive step count required")
        if self.first_rejection is not None or int(self.error.numpy()[0]):
            raise RuntimeError("F10H_FAIL_STOP: cannot advance a rejected state")
        for _ in range(count):
            start = self.accepted_steps * dt
            if self.control_index + 1 < len(self.control_plan):
                tick, event = self.control_plan[self.control_index + 1]
                if tick < self.accepted_steps:
                    raise RuntimeError("Missed control event")
                if tick == self.accepted_steps:
                    self.control(event)
                    self.control_index += 1
            self.audit_enabled = self.diagnostic
            self.previous_diag = self.diag.numpy().copy()
            self.previous_wall_stats = self.wall_stats.numpy().copy()
            before = self.x.numpy().copy()
            before_active = self.active.numpy().copy()
            failure = None
            try:
                FullParticleBatch.advance(self, dt, 1, graphs=graphs and not self.diagnostic)
            except RuntimeError as exc:
                failure = str(exc)
            q = self.quality()
            bits = int(self.error.numpy()[0])
            if bits:
                failure = "GPU_ERROR_" + str(bits)
            row = dict(
                step_number_one_based=self.physical_steps + 1,
                start_s=start,
                attempted_end_s=(self.physical_steps + 1) * dt,
                dt_s=dt,
                accepted=failure is None,
                error_bits=bits,
                reason=failure,
                control_index=self.control_index,
                quality=q,
            )
            if self.trace is not None:
                tr = self.trace
                pred = tr["predicted_positions_m"]
                pp = self.final_proposal.numpy().copy()
                tr["proposal_positions_m"] = pp
                tr.update(
                    row,
                    secondary_force_N=tr["total_force_N"] - tr["primary_force_N"],
                    positions_after_step_m=self.x.numpy().copy(),
                )
                d = np.linalg.norm(pp[:, None] - pp[None], axis=-1)
                d[np.diag_indices(self.n)] = np.inf
                pd = np.linalg.norm(pred[:, None] - pred[None], axis=-1)
                pd[np.diag_indices(self.n)] = np.inf
                row.update(
                    predicted_overlap_m=float(max(0, 2 * self.p.radius_m - pd.min())),
                    final_overlap_m=float(max(0, 2 * self.p.radius_m - d.min())),
                    round_residual_m=[r["max_overlap_m"] for r in tr["rounds"]],
                    round_max_correction_m=[r["max_correction_m"] for r in tr["rounds"]],
                    primary_max_N=float(np.linalg.norm(tr["primary_force_N"], axis=1).max()),
                    secondary_max_N=float(np.linalg.norm(tr["secondary_force_N"], axis=1).max()),
                    local_velocity_max_m_s=float(np.linalg.norm(tr["local_velocity_m_s"], axis=1).max()),
                    wall_correction_max_m=max(
                        float(np.linalg.norm(w["correction_m"], axis=1).max()) for w in tr["wall_passes"]
                    ),
                )
                tr.update(
                    prediction={"max_overlap_m": row["predicted_overlap_m"]},
                    final={"max_overlap_m": row["final_overlap_m"]},
                )
                self.all_steps.append(
                    {
                        key: tr[key]
                        for key in [
                            "initial_positions_m",
                            "positions_after_step_m",
                            "proposal_positions_m",
                            "predicted_positions_m",
                            "primary_force_N",
                            "total_force_N",
                            "free_velocity_m_s",
                            "flow_velocity_m_s",
                            "local_velocity_m_s",
                            "lubrication",
                            "start_s",
                            "attempted_end_s",
                            "dt_s",
                            "accepted",
                            "error_bits",
                        ]
                    }
                )
                self.traces.append(tr)
            self.ledger.append(row)
            self.physical_steps += 1
            if failure:
                self.first_rejection = dict(row, sample_id=self.sample_ids[0])
                self.last_legal = dict(
                    positions_m=before,
                    active=before_active,
                    diag=self.previous_diag,
                    wall_stats=self.previous_wall_stats,
                )
                self.last_legal_time_s = start
                # Original acceptance must not commit the rejected proposal.
                np.testing.assert_array_equal(self.x.numpy(), before)
                np.testing.assert_array_equal(self.active.numpy(), before_active)
                assert int(self.step_count.numpy()[0]) == self.accepted_steps
                raise RuntimeError(failure)
            self.accepted_steps += 1
            assert int(self.step_count.numpy()[0]) == self.accepted_steps

    def export_step_evidence(self, h):
        from scripts.fast_pzt_f9.common import clean

        h.attrs["actual_graph_chunk"] = 1
        h["substeps/accepted_dt_s"] = np.full(self.accepted_steps, self.mode_registry[self.numerical_mode][0])
        h["substeps/accepted_end_s"] = (
            np.arange(1, self.accepted_steps + 1) * self.mode_registry[self.numerical_mode][0]
        )
        if self.first_rejection is not None:
            for key, value in self.last_legal.items():
                h["last_legal_state/" + key] = value.reshape(1, self.n, *value.shape[1:])
            h["last_legal_state"].attrs["time_s"] = self.last_legal_time_s
            h["rejected_interval"].attrs["start_s"] = self.last_legal_time_s
            h["rejected_interval"].attrs["attempted_end_s"] = self.first_rejection["attempted_end_s"]
            h["rejected_interval"].attrs["first_rejection_json"] = json.dumps(clean(self.first_rejection))

    def export_observer(self, out):
        from scripts.fast_pzt_f9.common import write, clean

        meta = dict(
            numerical_mode=self.numerical_mode,
            radius_m=self.p.radius_m,
            physics=self.p.metadata(),
            first_rejection_step=self.first_rejection,
            attempted_steps=self.physical_steps,
            accepted_steps=self.accepted_steps,
            traced_steps=len(self.traces),
            graph_chunk=1,
            diagnostic=self.diagnostic,
            original_quality=self.quality(),
            training_ready=False,
        )
        write(out / "OBSERVER.json", meta)
        for name, data in [
            ("CONTACT_TRACE", dict(meta, traces=list(self.traces))),
            ("SUBSTEP_LEDGER", dict(meta, steps=self.ledger)),
            ("ALL_STEPS", self.all_steps),
        ]:
            with gzip.open(out / (name + ".json.gz"), "wt") as stream:
                json.dump(clean(data), stream, separators=(",", ":"), allow_nan=False)
