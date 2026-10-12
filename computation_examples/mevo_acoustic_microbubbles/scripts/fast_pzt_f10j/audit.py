"""Read-only F10I reproduction and per-iteration instrumentation, no solver changes."""

import argparse
import inspect
import numpy as np
import warp as wp
import warp._src.optim.linear as linear
from warp.optim.linear import LinearOperator
from mevo_mb.fast_sim.f10i_gpu import FixedStateB2, product, vec
from .common import ROOT, BASE, js, write, identity, sha


@wp.kernel
def trace_linear(
    x: wp.array(dtype=wp.float64),
    b: wp.array(dtype=wp.float64),
    rec: wp.array(dtype=wp.float64),
    true: wp.array(dtype=wp.float64),
    inv: wp.array(dtype=wp.mat33d),
    p: wp.array(dtype=wp.float64),
    ap: wp.array(dtype=wp.float64),
    counter: wp.array(dtype=wp.int32),
    trace: wp.array2d(dtype=wp.float64),
    solutions: wp.array2d(dtype=wp.float64),
):
    k = counter[0]
    rr = wp.float64(0.0)
    tt = wp.float64(0.0)
    bb = wp.float64(0.0)
    pr = wp.float64(0.0)
    pt = wp.float64(0.0)
    curvature = wp.float64(0.0)
    drift = wp.float64(0.0)
    for i in range(inv.shape[0]):
        r = vec(rec, i)
        t = vec(true, i)
        rr += wp.dot(r, r)
        tt += wp.dot(t, t)
        bb += wp.dot(vec(b, i), vec(b, i))
        pr += wp.dot(r, inv[i] * r)
        pt += wp.dot(t, inv[i] * t)
        curvature += wp.dot(vec(p, i), vec(ap, i))
        drift += wp.dot(r - t, r - t)
    trace[k, 0] = wp.sqrt(rr)
    trace[k, 1] = wp.sqrt(tt)
    trace[k, 2] = wp.sqrt(bb)
    trace[k, 3] = wp.sqrt(wp.max(pr, wp.float64(0.0)))
    trace[k, 4] = wp.sqrt(wp.max(pt, wp.float64(0.0)))
    trace[k, 5] = curvature
    trace[k, 6] = wp.sqrt(drift)
    for j in range(x.shape[0]):
        solutions[k, j] = x[j]
    counter[0] = k + 1


class ObservedB2(FixedStateB2):
    """Inherit the entire original solve; add only independent trace buffers."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.evaluations = []
        self.traces = {}
        for kind, solver, inv, contact in [("A", self.cgA, self.invA, 0), ("H", self.cgH, self.invH, 1)]:
            count = 129 if kind == "A" else 25 * 129
            trace = wp.zeros((count, 7), dtype=wp.float64, device=self.device)
            solutions = wp.zeros((count, 3 * self.n), dtype=wp.float64, device=self.device)
            counter = wp.zeros(1, dtype=wp.int32, device=self.device)
            true = wp.zeros(3 * self.n, dtype=wp.float64, device=self.device)
            self.traces[kind] = (counter, trace, solutions)
            self.instrument(solver, inv, contact, true, counter, trace, solutions)
            with wp.ScopedCapture(device=self.device) as cap:
                result = solver()
            setattr(self, "linear" + kind, result)
            setattr(self, "graph" + kind, cap.graph)

    def instrument(self, solver, inv, contact, true, counter, trace, solutions):
        original = solver._M.matvec

        def mv(x, y, z, alpha, beta):
            original(x, y, z, alpha, beta)
            self.launch(
                product,
                self.n,
                [solver._x, solver._b, true, self.c, self.normal, self.a, self.current, contact, -1.0, 1.0],
            )
            self.launch(
                trace_linear,
                1,
                [
                    solver._x,
                    solver._b,
                    x,
                    true,
                    inv,
                    solver._p_and_Ap[0],
                    solver._p_and_Ap[1],
                    counter,
                    trace,
                    solutions,
                ],
            )

        solver._M = LinearOperator(solver._M.shape, wp.float64, self.device, mv)

    def evaluate(self, q):
        result = super().evaluate(q)
        # Caller variables are read solely to label the original algorithm's decisions.
        loc = inspect.currentframe().f_back.f_locals
        self.evaluations.append(
            dict(
                outer=loc.get("outer"),
                line=loc.get("line") if q is self.trial else None,
                rate=loc.get("rate") if q is self.trial else None,
                trial=q is self.trial,
                stats=result.copy(),
                q=q.numpy(),
                pred=self.pred.numpy(),
                direction=self.direction.numpy(),
                rhs=self.rhs.numpy(),
                force=self.force.numpy(),
                active=self.a.numpy(),
                normals=self.current.numpy(),
            )
        )
        return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    out = BASE / args.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    path = ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json"
    row = js(path)["rows"][-1]
    x, free = np.array(row["initial_positions_m"]), np.array(row["free_velocity_m_s"])
    wp.init()
    original = FixedStateB2(32)
    plain = original.solve(x, free)
    write(
        out / "PLAIN.json",
        dict(
            result=plain,
            input_sha256=sha(path),
            allocation_s=original.allocation_s,
            graph_build_s=original.graph_build_s,
        ),
    )
    observed = ObservedB2(32)
    result = observed.solve(x, free)
    traces = {}
    for kind, (counter, trace, solutions) in observed.traces.items():
        count = int(counter.numpy()[0])
        traces[kind] = dict(
            columns=[
                "recursive_abs",
                "true_abs",
                "rhs_norm",
                "recursive_preconditioned",
                "true_preconditioned",
                "previous_pAp",
                "recursive_true_drift",
            ],
            values=trace.numpy()[:count],
            solutions=solutions.numpy()[:count],
        )
    write(
        out / "AUDIT.json",
        dict(
            result=result,
            evaluations=observed.evaluations,
            linear=traces,
            instrumentation_position_delta_m=float(np.max(abs(result["positions_m"] - plain["positions_m"]))),
            original_saved_position_delta_m=float(
                np.max(abs(plain["positions_m"] - np.array(row["gpu"]["positions_m"])))
            ),
            instrumented_timing_is_not_performance=True,
        ),
    )
    source = inspect.getfile(linear)
    write(
        out / "WARP_CAPABILITIES.json",
        dict(
            version=wp.__version__,
            linear_source_sha256=sha(source),
            cg_signature=str(inspect.signature(linear.cg)),
            restart_supported=True,
            capture_while_supported=hasattr(wp, "capture_while"),
            actual_cg_source=inspect.getsource(linear.CG),
            actual_loop_source=inspect.getsource(linear._run_capturable_loop),
        ),
    )
    print("Original/observed", plain["reason"], result["reason"], flush=True)


if __name__ == "__main__":
    main()
