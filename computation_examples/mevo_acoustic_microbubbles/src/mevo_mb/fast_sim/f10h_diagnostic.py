"""B=1 CPU-contact/GPU-field diagnostic adapter; NOT a production GPU solver.

Original acoustics, pair forces, CFD, final acceptance and fail-stop are reused.
Only wall-free steps are eligible: any wall projection requests a stop. No PBD
is applied to a soft-contact proposal. CUDA capture is intentionally disabled.
"""

from time import perf_counter
import gzip
import json
import numpy as np
from . import kernels as k
from .f4_engine import F4VesselSimulator
from .f10h_lubrication import LubricationDiagnostic
from .f10g_contact import resistance, overlaps
from .f10h_contact import geometric_contact


class GeometricContactDiagnostic(LubricationDiagnostic):
    mode_registry = {"B2_GEOM_K2_DT10": (1e-5, 2), "B2_GEOM_K2_DT5": (5e-6, 2)}
    registered_physical_steps = (1e-5, 5e-6)

    def __init__(self, *args, kc_N_m, **kwargs):
        kwargs.setdefault("numerical_mode", "B2_GEOM_K2_DT10")
        kwargs["diagnostic"] = True
        self.kc_N_m = float(kc_N_m)
        self.soft_steps = []
        super().__init__(*args, **kwargs)

    def graph(self, dt, count):
        # Existing recorder preflights graph() even when graphs=False. Returning
        # None performs no numerical work; advance below forbids graph replay.
        if count != 1 or dt not in self.registered_physical_steps:
            raise ValueError("Only one registered 10/5us diagnostic step")
        return None

    def advance(self, dt, count, graphs=False):
        if graphs:
            raise ValueError("F10H CPU reference cannot be CUDA captured")
        return super().advance(dt, count, graphs=False)

    def step(self, dt):
        self.trace = None
        self._current_dt = dt
        self._pair_contact_pass = 0
        begin = perf_counter()
        x = self.x.numpy().copy()
        self.fields()
        primary = self.force.numpy().copy()
        self.pairs()
        self.launch(
            k.free_velocity,
            self.total,
            [self.force, self.free, 6 * np.pi * self.p.viscosity_Pa_s * self.p.radius_m],
        )
        free = self.free.numpy()
        force = self.force.numpy()
        flow = self.flow_query.velocity.numpy()
        field_and_read_s = perf_counter() - begin
        t = perf_counter()
        A = resistance(x, self.p.radius_m)["A"]
        assembly_s = perf_counter() - t
        result = geometric_contact(
            x, free, self.p.radius_m, 6 * np.pi * self.p.viscosity_Pa_s * self.p.radius_m, dt, self.kc_N_m, A
        )
        walls = []
        # Keep original wall kernels, but stop if they would modify this interior
        # reference. Do not silently combine wall projection with the spring.
        for label, proposal in [
            ("prediction", result["prediction_m"]),
            ("soft_contact", result["positions_m"]),
        ]:
            self.tmp1.assign(proposal)
            F4VesselSimulator.wall_pass(self, self.tmp1)
            correction = self.tmp1.numpy() - proposal
            walls.append(dict(stage=label, after_pair_round=0, correction_m=correction))
        before_commit = perf_counter()
        diag = self.diag.numpy()
        pred = result["prediction_m"]
        distances = np.linalg.norm(pred[:, None] - pred[None, :], axis=-1)
        np.fill_diagonal(distances, np.inf)
        diag[:, 2] = np.maximum(diag[:, 2], np.maximum(0, 2 * self.p.radius_m - distances.min(axis=1)))
        physical = np.linalg.norm(pred - x, axis=1)
        diag[:, 5] += physical
        diag[:, 7] = np.maximum(diag[:, 7], physical / dt)
        diag[:, 9] = np.maximum(diag[:, 9], physical)
        rounds = []
        iterate = pred.copy()
        for row in result["iterations"]:
            dx = row["correction_m"]
            iterate += dx
            diag[:, 4] += np.linalg.norm(dx, axis=1)
            diag[:, 6] = np.maximum(diag[:, 6], np.linalg.norm(dx, axis=1))
            rounds.append(
                dict(
                    round=row["iteration"],
                    positions_m=iterate.copy(),
                    correction_m=dx,
                    max_correction_m=row["max_correction_m"],
                    max_overlap_m=float(overlaps(iterate, self.p.radius_m).max()),
                    interpretation="Implicit Newton iterate, not PBD",
                )
            )
        self.diag.assign(diag)
        extra_bits = 0
        if not result["converged"] or any(np.any(w["correction_m"] != 0) for w in walls):
            extra_bits = 1048576
        if (
            result["geometric_overlap_m"].max() > 0.01 * self.p.radius_m
            or result["max_single_iteration_m"] > 0.1 * self.p.radius_m
        ):
            extra_bits |= 2048
        if extra_bits:
            self.error.assign(np.array([int(self.error.numpy()[0]) | extra_bits], dtype=np.int32))
        self.final_proposal = self.tmp1
        self.launch(
            k.accept,
            self.total,
            [
                self.tmp1,
                self.x,
                self.n,
                self.p.radius_m,
                self.diag,
                self.error,
                self.adjacency,
                self.event_count,
                self.event_ids,
                self.event_time,
                self.step_count,
                self.time_origin,
                dt,
            ],
        )
        self.launch(k.advance_clock, 1, [self.step_count])
        final = self.tmp1.numpy().copy()
        self.trace = dict(
            initial_positions_m=x,
            dt_s=dt,
            primary_force_N=primary,
            total_force_N=force,
            free_velocity_m_s=free,
            flow_velocity_m_s=flow,
            local_velocity_m_s=result["mobility_free_m_s"],
            predicted_positions_m=pred,
            proposal_positions_m=final,
            rounds=rounds,
            wall_passes=walls,
            lubrication=[],
            soft_contact=result,
            algorithm="B2 current geometric normal implicit spring with exact discrete lubrication; CPU diagnostic",
        )
        self.soft_steps.append(
            dict(
                step=len(self.soft_steps) + 1,
                initial_positions_m=x,
                result=result,
                field_and_read_s=field_and_read_s,
                resistance_assembly_s=assembly_s,
                commit_and_read_s=perf_counter() - before_commit,
                full_step_s=perf_counter() - begin,
                error_bits=int(self.error.numpy()[0]),
                accepted=not bool(int(self.error.numpy()[0])),
                actual_accepted_positions_m=self.x.numpy().copy(),
                cfd_velocity_m_s=flow,
                acoustic_free_velocity_m_s=force / (6 * np.pi * self.p.viscosity_Pa_s * self.p.radius_m),
            )
        )

    def export_observer(self, out):
        from scripts.fast_pzt_f9.common import clean, write

        super().export_observer(out)
        p = out / "OBSERVER.json"
        meta = json.loads(p.read_text())
        meta.update(
            contact_model="F10H_B2_GEOMETRIC_IMPLICIT",
            kc_N_m=self.kc_N_m,
            production=False,
            numerical_mode=self.numerical_mode,
            graph_enabled=False,
            wall_scope="Any nonzero wall correction rejects this interior-only prototype",
        )
        write(p, meta)
        with gzip.open(out / "SOFT_STEPS.json.gz", "wt") as f:
            json.dump(clean(self.soft_steps), f, separators=(",", ":"), allow_nan=False)
