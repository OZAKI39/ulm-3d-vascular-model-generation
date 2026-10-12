"""Bounded failure attribution, real saved-state scope, cache invalidation and memory."""

import argparse
import ctypes
import gc
import h5py
import numpy as np
import warp as wp
from .common import BASE, ROOT, case, objects, identity, write
from mevo_mb.fast_sim.f10c_engine import FullParticleBatch
from mevo_mb.fast_sim.f10c_multirate import MultirateBatch
from mevo_mb.fast_sim.f10c_admission import admit_cache


def memory_used():
    class Memory(ctypes.Structure):
        _fields_ = [("total", ctypes.c_ulonglong), ("free", ctypes.c_ulonglong), ("used", ctypes.c_ulonglong)]

    lib = ctypes.CDLL("libnvidia-ml.so.1")
    assert lib.nvmlInit_v2() == 0
    handle = ctypes.c_void_p()
    assert lib.nvmlDeviceGetHandleByIndex_v2(0, ctypes.byref(handle)) == 0
    m = Memory()
    assert lib.nvmlDeviceGetMemoryInfo(handle, ctypes.byref(m)) == 0
    return m.used / 2**20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    obj = objects()
    p, b, r, g, f = obj
    cfg = case("O_1000")
    x = np.array(cfg["environment"]["positions_m"])
    event = cfg["events"][0]
    saved = []
    cases = [("interior", x, event), ("sound_off", x, case("OFF")["events"][0])]
    for name, rel, frame in [
        ("F9_nearwall", "f9/y1_10_v1/run", 0),
        ("F5_high_shear", "f5/science_v1/main", 600),
        ("F5_dispersed", "f5/science_v1/main", -1),
    ]:
        with h5py.File(ROOT / "results/fast_pzt" / rel / "trajectory.h5") as h:
            xx = h["positions_m"][frame, 0]
        cases.append((name, xx, event))
    for name, xx, ev in cases:
        row = admit_cache(xx, ev, obj)
        row["case"] = name
        saved.append(row)
    assert (
        saved[0]["status"] == "ADMITTED_FOR_NUMERICAL_TEST"
        and saved[1]["status"] == "ADMITTED_FOR_NUMERICAL_TEST"
    )
    assert all(row["status"] == "MULTIRATE_NOT_ADMITTED" for row in saved[2:])
    # An injected CFD-domain fault in member 2 must stop all members, with ID attribution.
    s = FullParticleBatch(np.stack([x] * 4), p, r, b, flow=f, vessel=g, sample_ids=["A", "B", "bad_C", "D"])
    s.control(case("OFF")["events"][0])
    bad = s.x.numpy()
    bad[2 * 32] = [1.0, 1.0, 1.0]
    s.x.assign(bad)
    s.advance(1e-5, 1, graphs=False)
    assert int(s.error.numpy()[0]) & 16384
    failure = s.failure_members()
    assert failure[2]["status"] == "FAILED_MEMBER"
    assert all(row["status"] != "PASS" for row in failure)
    np.testing.assert_array_equal(s.x.numpy(), bad)
    del s
    gc.collect()
    # Exercise actual stale-cache guards, not only CPU preflight labels.
    s = MultirateBatch(x, p, r, b, flow=f, vessel=g, field_stride=2, flow_stride=2)
    s.control(event)
    s.fields()
    s.build_phase = 1
    s.grid.fill_(wp.vec2d(1e9, 0.0))
    s.fields()
    table_error = int(s.error.numpy()[0])
    assert table_error & 2
    s.reset()
    s.control(event)
    s.fields()
    s.build_phase = 1
    bad = x.copy()
    bad[0, 0] = b.origin[0] - 1e-4
    s.x.assign(bad)
    s.fields()
    field_error = int(s.error.numpy()[0])
    assert field_error & 1
    s.reset()
    s.control(event)
    s.fields()
    s.build_phase = 1
    s.flow.launch(s.x, s.active, s.pending, s.flow_query, s.error)
    # Accepted-position validity cache no longer covers this artificial stale position.
    bad = x.copy()
    bad[0, 0] += 1e-2
    s.x.assign(bad)
    from mevo_mb.fast_sim import kernels as k

    s._current_dt = 1e-5
    s.launch(k.free_velocity, s.total, [s.force, s.free, 6 * np.pi * p.viscosity_Pa_s * p.radius_m])
    cfd_error = int(s.error.numpy()[0])
    assert cfd_error & 65536
    switches = []
    for key in ["O_1000", "OFF", "S_1000"]:
        s.reset()
        before = s.x.numpy()
        record = s.control(case(key)["events"][0])
        assert s.phase == 0 and not s.graphs
        np.testing.assert_array_equal(before, s.x.numpy())
        s.advance(1e-5, 1)
        assert int(s.error.numpy()[0]) == 0
        switches.append(
            dict(
                control=key,
                cache_invalidated=record["multirate_cache_invalidated"],
                first_step_quality=s.quality(),
            )
        )
    # dt changes invalidate cached phase and graphs before the next physical step.
    before = s.capture_s
    s.advance(5e-6, 1)
    assert s.phase == 1 and s.last_dt == 5e-6
    dt_check = dict(status="PASS", new_dt_s=s.last_dt, new_capture_s=s.capture_s - before)
    del s
    gc.collect()
    memory = []
    for size in [1, 4, 16, 64]:
        before = memory_used()
        s = FullParticleBatch(np.stack([x] * size), p, r, b, flow=f, vessel=g, graph_chunk=5)
        s.control(event)
        s.graph(1e-5, 5)
        wp.synchronize()
        prepared = memory_used()
        s.advance(1e-5, 5)
        after = memory_used()
        memory.append(
            dict(
                B=size,
                before_MiB=before,
                resident_after_capture_MiB=prepared,
                after_5_steps_MiB=after,
                observed_peak_MiB=max(before, prepared, after),
                warp_pool_high_MiB=wp.get_mempool_used_mem_high("cuda:0") / 2**20,
                scope="Separate bounded memory validation, same arrays/graph as full10ms run; NVML samples, not a continuous exact peak. Shared allocator/context retention included.",
            )
        )
        del s
        gc.collect()
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            scope_checks=saved,
            injected_failure_members=failure,
            fault_injection_scope="SOFTWARE_TEST_ONLY: invalid positions/grid injected without claiming a physical trajectory",
            cache_guard_errors=dict(
                response_table=table_error, field_domain=field_error, stale_CFD=cfd_error
            ),
            controls=switches,
            dt_change=dt_check,
            memory=memory,
            training_ready=False,
        ),
    )
    print("PASS", flush=True)


if __name__ == "__main__":
    main()
