"""Load immutable reference IDs/topology and export native LAMMPS input formats."""

from pathlib import Path
import numpy as np
from .config import INPUTS, RHO, V
from .io import load_npz


def load_reference():
    ref = load_npz(INPUTS / "reference.npz")
    n = len(ref["X"])
    nb = len(ref["bonds"])
    assert np.array_equal(ref["atom_ids"], np.arange(1, n + 1))
    assert np.array_equal(ref["original_bond_ids"], np.arange(1, nb + 1))
    assert (ref["V"] > 0).all()
    return ref


def write_native_inputs(path, ref):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    N = len(ref["X"])
    n = int(ref["clot_count"])
    b = ref["bonds"]
    NB = len(b)
    with (path / "reference.native").open("x") as f:
        f.write(f"LPSVAR_REFERENCE_V2\n{N} {NB} {V:.17g}\n")
        for i, (m, v) in enumerate(zip(ref["m"], ref["V"])):
            f.write(f"{i+1} {m:.17g} {v:.17g}\n")
        for k, ((i, j), r, h, bt, cut) in enumerate(
            zip(b, ref["r0"], ref["taper"], ref["bond_types"], ref["initial_cut"])
        ):
            f.write(
                f"{k+1} {i+1} {j+1} {int(bt)} {r:.17g} {h:.17g} {float(cut)*.30:.17g} {int(cut)} {int(cut)}\n"
            )
    with (path / "body.data").open("x") as f:
        f.write(
            f"Irregular Y clot V4; original IDs preserved\n\n{N} atoms\n{NB} bonds\n\n2 atom types\n3 bond types\n\n"
        )
        for axis in "xyz":
            f.write(f"-0.03 0.03 {axis}lo {axis}hi\n")
        f.write(f"\nMasses\n\n1 {RHO*V:.17g}\n2 {RHO*V:.17g}\n\nAtoms # bond\n\n")
        for k, x in enumerate(ref["X"]):
            f.write(
                f"{k+1} 1 {1 if k<n else 2} " + " ".join(f"{v:.17g}" for v in x) + "\n"
            )
        f.write("\nBonds\n\n")
        for k, ((i, j), bt) in enumerate(zip(b, ref["bond_types"])):
            f.write(f"{k+1} {int(bt)} {i+1} {j+1}\n")
    return path
