"""Explicit single-rank native LAMMPS V4 driver; physical N/J, original-ID state.

The C++ material commits history in the native dynamic coordinator. Python only
predicts event substeps, distributes external forces, and audits accepted states.
"""

from pathlib import Path
import ctypes
import types
import traceback
import numpy as np
from mpi4py import MPI
from lammps import lammps
from .config import LIB, INPUTS, BASE, DX, V, RHO, K, G, HORIZON, S0, SC
from .io import load_npz, sha, verify_native_build
from .geometry import load_reference, write_native_inputs


class NativeSolver:

    def __init__(
        self,
        case,
        restart=None,
        protocol="biphasic",
        install_properties=True,
        surface=True,
        dt=BASE,
    ):
        verify_native_build()
        case = Path(case)
        case.mkdir(parents=True, exist_ok=True)
        (case / "keys").mkdir(exist_ok=True)
        name = "native"
        restart_state = None
        query = None
        self.comm = MPI.COMM_WORLD
        assert self.comm.size == 1
        self.ref = dict(load_reference())
        self.N = len(self.ref["X"])
        self.n = int(self.ref["clot_count"])
        self.NB = len(self.ref["bonds"])
        self.dt = dt
        self.dynamic = False
        self.restart_state = restart_state
        self.query = query
        self.loaded_library = ctypes.CDLL(str(LIB), mode=ctypes.RTLD_GLOBAL)
        if not restart:
            write_native_inputs(case / "provenance", self.ref)
        self.l = lammps(
            cmdargs=["-log", str(case / (name + ".log")), "-screen", "none"],
            comm=self.comm,
        )
        self.lib = ctypes.CDLL(str(LIB))
        self.lib.lpsvar_action.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            ctypes.c_longlong,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
        ]
        self.lib.lpsvar_action.restype = ctypes.c_int
        self.lib.lpsvar_output.argtypes = [
            ctypes.c_void_p,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
        ]
        self.lib.lpsvar_output.restype = ctypes.c_int
        if restart:
            self.l.commands_string(
                f"atom_modify map array sort 0 0\nread_restart {restart}"
            )
        else:
            self.l.commands_string(
                f"units si\nboundary m m m\natom_style bond\natom_modify map array sort 0 0\nnewton on off\nspecial_bonds lj 0 1 1 coul 1 1 1\nread_data {case / 'provenance/body.data'} extra/bond/per/atom 160 extra/special/per/atom 180\n"
            )
        self.l.commands_string(
            f"lattice sc {DX:.17g}\nneighbor 0.0001 bin\nneigh_modify every 1 delay 0 check no\ncomm_modify cutoff {3 * HORIZON:.17g}\nfix props all property/atom d_vfrac d_theta ghost yes\nset group all d_vfrac {V:.17g}\n"
        )
        self.l.command("fix nodal_mass all property/atom rmass ghost yes")
        ids = self.l.numpy.extract_atom("id")[: self.N].astype(int) - 1
        self.l.numpy.extract_atom("rmass")[: self.N] = RHO * self.ref["V"][ids]
        self.l.numpy.extract_atom("d_vfrac")[: self.N] = self.ref["V"][ids]
        if not restart:
            self.l.commands_string(
                f"pair_style bpm/peri contact_distance_factor 1.0\npair_coeff * * {18 * K / (np.pi * HORIZON ** 4):.17g} {HORIZON:.17g}\nbond_style bpm/peri/variational {case / 'provenance/reference.native'}\n"
            )
            for t in [1, 2, 3]:
                self.l.command(
                    f"bond_coeff {t} {K:.17g} {G:.17g} {HORIZON:.17g} {S0[t]:.17g} {SC[t]:.17g}"
                )
        self.l.commands_string(
            f"special_bonds lj 0 1 1 coul 1 1 1\ntimestep {dt:.17g}\nthermo 100000\nthermo_style custom step time atoms bonds pe ebond epair ke\nthermo_modify lost error norm no\n"
        )
        self.maps = sorted(
            set(
                (
                    s.split()[-1]
                    for s in Path("/proc/self/maps").read_text().splitlines()
                    if "liblammps" in s
                )
            )
        )
        assert self.maps and all((sha(p) == sha(LIB) for p in self.maps))
        self.vol = self.ref["V"]
        self.mass = RHO * self.vol
        if restart is None:
            on = S0[self.ref["bond_types"]].copy()
            off = SC[self.ref["bond_types"]].copy()
            if install_properties:
                if protocol == "biphasic":
                    field = load_npz(INPUTS / "material.npz")
                    on = field["bond_s0"]
                    off = field["bond_sc"]
                ids = np.arange(1, self.NB + 1, dtype=np.int64)[::-1].copy()
                self.lib.lpsvar_set_thresholds.argtypes = [
                    ctypes.c_void_p,
                    ctypes.c_int,
                    ctypes.POINTER(ctypes.c_longlong),
                    ctypes.POINTER(ctypes.c_double),
                    ctypes.POINTER(ctypes.c_double),
                ]
                self.lib.lpsvar_set_thresholds.restype = ctypes.c_int
                assert (
                    self.lib.lpsvar_set_thresholds(
                        self.l.lmp,
                        self.NB,
                        ids.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),
                        self.ptr(on[::-1].copy()),
                        self.ptr(off[::-1].copy()),
                    )
                    == 0
                )
            d = self.ref["initial_cut"].astype(float)
            k = 0.3 * d
            assert (
                self.lib.lpsvar_action(self.l.lmp, 1, 0, self.ptr(k), self.ptr(d)) == 0
            )
        assert np.array_equal(self.field(15), self.vol)
        if surface:
            from .surface import CurrentSurface

            self.tool = CurrentSurface(self.ref, "uniform")
        self.ext = types.SimpleNamespace(error=None, last={})
        self.force_stage_dt = 0.0

    def ptr(self, a):
        return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    def field(self, i):
        size = {
            0: 3 * self.N,
            1: 3 * self.N,
            2: self.NB,
            3: self.NB,
            4: self.N,
            5: 2 * self.NB,
            6: 2 * self.NB,
            7: self.N,
            8: 10,
            9: self.NB,
            10: self.NB,
            11: 5,
            12: self.NB,
            13: self.NB,
            14: self.NB,
            15: self.N,
        }[i]
        out = np.empty(size)
        assert self.lib.lpsvar_output(self.l.lmp, i, self.ptr(out)) == size
        return out

    def gather(self, key):
        n = self.l.extract_global("nlocal")
        ids = self.l.numpy.extract_atom("id")[:n].copy()
        vals = self.l.numpy.extract_atom(key)[:n].copy()
        assert n == self.N and np.array_equal(np.sort(ids), np.arange(1, self.N + 1))
        return vals[np.argsort(ids)]

    def state(self):
        if self.dynamic and self.ext.error:
            raise RuntimeError(self.ext.error)
        d = self.field(3)
        td = self.field(6)
        failed = d == 1
        z = dict(
            step=int(self.l.get_thermo("step")),
            time=float(self.l.get_thermo("time")),
            x=self.gather("x"),
            velocity=self.gather("v"),
            force=self.gather("f"),
            material_force=self.field(0).reshape(-1, 3),
            kappa=self.field(2),
            d=d,
            scalars=self.field(8),
            counters=self.field(11),
            native_bonds=int(self.l.get_thermo("bonds")),
            Uc=float(self.l.get_thermo("epair")),
            K=float(self.l.get_thermo("ke")),
            native_U=float(self.l.get_thermo("ebond")),
            failed_tension_max_N=float(
                np.max(abs(td[np.r_[failed, failed]]), initial=0)
            ),
            atom_ids=self.gather("id").astype(np.int64),
            box_lo=np.asarray(self.l.extract_box()[0]),
            box_hi=np.asarray(self.l.extract_box()[1]),
        )
        if self.dynamic:
            z.update({k: v.copy() for k, v in self.ext.last.items()})
        return z

    def static(self, x=None):
        if x is not None:
            ids = self.l.numpy.extract_atom("id")[: self.N]
            self.l.numpy.extract_atom("x")[: self.N] = x[ids - 1]
        self.l.command("run 0")
        return self.state()

    def close(self):
        self.l.close()

    def step(self, h):
        if h != self.dt:
            self.l.command(f"timestep {h:.17g}")
            self.dt = h
        self.force_stage_dt = h
        try:
            self.l.command("run 1 pre no post no")
        finally:
            self.force_stage_dt = 0.0
        return self.state()

    def dynamics(self):
        self.l.commands_string(
            f"group clot type 1\nfix integrate clot nve\nfix coordinator all LPS/VARIATIONAL/DYNAMIC {RHO:.17g}\nfix v3forces clot external pf/callback 1 1"
        )

        def callback(caller, step, nlocal, tags, x, f):
            f[:] = 0
            try:
                ids = tags.astype(int) - 1
                assert nlocal == self.N and np.array_equal(
                    np.sort(ids), np.arange(self.N)
                )
                current = np.empty((self.N, 3))
                velocity = np.empty_like(current)
                current[ids] = x
                velocity[ids] = self.l.numpy.extract_atom("v")[: self.N]
                t = float(self.l.get_thermo("time"))
                d = self.field(3)
                drift = velocity.copy()
                iterations = 0
                if self.force_stage_dt:
                    internal = np.empty_like(current)
                    internal[ids] = self.l.numpy.extract_atom("f")[: self.N]
                    factor = self.force_stage_dt / (2 * self.mass[: self.n, None])
                    for iterations in range(1, 13):
                        z = self.tool.evaluate(current, velocity, d, t)
                        updated = drift.copy()
                        updated[: self.n] += factor * (
                            internal[: self.n]
                            + z["force"][: self.n]
                            + z["wall_force"][: self.n]
                        )
                        err = np.max(abs(updated - velocity))
                        velocity = updated
                        if err < 2e-13 * max(1.0, float(np.max(abs(velocity)))):
                            break
                    else:
                        raise RuntimeError("Follower slip half-kick did not converge")
                z = self.tool.evaluate(current, velocity, d, t)
                f[:] = (z["force"] + z["wall_force"])[ids]
                self.ext.last = dict(
                    external_force=z["force"] + z["wall_force"],
                    external_wall_force=z["wall_force"],
                    wall_potential_J=np.array(z["wall_U"]),
                    wall_clearance=z["wall_clearance"],
                    wall_damping_power_W=np.array(z["wall_damping_power"]),
                    external_drive_force=z["drive_force"],
                    external_drag_force=z["drag_force"],
                    external_accessible_force=np.asarray(
                        z["W"] @ (z["face_force"] * z["accessible"][:, None])
                    ),
                    external_L2_force=np.asarray(
                        z["W"] @ (z["face_force"] * (self.tool.method == 1)[:, None])
                    ),
                    external_time=np.array(t),
                    force_stage_velocity=velocity,
                    verlet_drift_velocity=drift,
                    slip_iterations=np.array(iterations),
                    surface_centroid=z["centroid"],
                    surface_normal=z["normal"],
                    surface_area=z["area"],
                    surface_active=z["active"],
                    surface_velocity=z["face_velocity"],
                    surface_outer_velocity=z["outer_velocity"],
                    surface_traction=z["traction"],
                    surface_force=z["face_force"],
                    surface_weight=z["weight"],
                    surface_original=z["original"],
                    surface_new_open=z["new_open"],
                    surface_accessible=z["accessible"],
                    surface_unresolved=z["unresolved"],
                    surface_aperture_m=z["aperture_m"],
                    surface_fluid_region=z["fluid_region"],
                    surface_occluded=z["occluded"],
                    surface_method=self.tool.method,
                    surface_class=z["face_class"],
                    field_amplitude_Pa=np.array(self.tool.field.amplitude(t)),
                )
            except BaseException:
                self.ext.error = traceback.format_exc()
                f[:] = np.nan

        self.callback = callback
        self.l.set_fix_external_callback("v3forces", callback)
        self.dynamic = True
        self.l.command("run 0")
        return self.state()
