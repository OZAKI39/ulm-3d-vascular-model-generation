"""New geometry-identified k-Wave basis; direct complex superposition verification."""
import json
import h5py
import numpy as np
import yaml
from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.acoustics.kwave_backend import KWaveCW, precompute
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from mevo_mb.utils.config import load_config
from scripts.stage4.common import ROOT, new_run, require, sha, write, finish


def configuration():
    config = load_config(ROOT / "configs/stage4/acoustic.yaml")
    geometry = require("geometry")
    config["acoustics"]["geometry_hash"] = geometry["geometry_hash"]
    return config


def main():
    config = configuration()
    out = new_run("acoustic",config,["simvascular/import"])
    (out / "acoustic").mkdir()
    path = precompute(config,out)
    # The new basis itself carries the Y identity; existing cylinder files are untouched.
    with h5py.File(path,"r+") as h:
        metadata = json.loads(h.attrs["metadata_json"])
        metadata["geometry_hash"] = config["acoustics"]["geometry_hash"]
        metadata["medium_scope"] = "Homogeneous water; lumen geometry identifies control domain, no wall acoustic interfaces"
        h.attrs["metadata_json"] = json.dumps(metadata)
    db = AcousticBasisDatabase.load_basis(path)
    check = out / "superposition"
    check.mkdir()
    array = config["acoustics"]
    u = array["basis_source_amplitude_Pa"]*np.asarray(array["superposition_amplitudes"])*np.exp(1j*np.asarray(array["superposition_phases_rad"]))
    direct,gradient,meta = KWaveCW(config,check).run(u,"direct")
    combined = db.combine(u)
    geometry = BifurcationGeometry(yaml.safe_load((ROOT / "configs/stage4/bifurcation.yaml").read_text())["geometry"])
    xyz = np.stack(np.meshgrid(*db.axes_m,indexing="ij"),axis=-1)
    mask = geometry.contains(xyz.reshape(-1,3)).reshape(xyz.shape[:-1])
    p_error = float(np.linalg.norm((combined.pressure_Pa-direct)[mask])/np.linalg.norm(direct[mask]))
    g_error = float(np.linalg.norm((combined.gradient_Pa_m-gradient)[mask])/np.linalg.norm(gradient[mask]))
    significant = mask & (np.abs(direct) > np.max(np.abs(direct[mask]))*.01)
    phase = np.angle(combined.pressure_Pa[significant]*np.conj(direct[significant]))
    checks = dict(pressure_superposition=p_error<array["validation_relative_L2_max"],
                  gradient_superposition=g_error<array["validation_relative_L2_max"],
                  finite=bool(np.isfinite(db.pressure_basis).all() and np.isfinite(db.gradient_basis).all()),
                  periodic=max(m["recording_half_window_relative_L2"] for m in db.metadata["simulations"]+[meta])<.01)
    write(out / "SUPERPOSITION.json",dict(checks=checks,pressure_relative_L2=p_error,gradient_relative_L2=g_error,
                                           masked_phase_RMS_rad=float(np.sqrt(np.mean(phase**2))),count=int(mask.sum()),direct_metadata=meta))
    finish(out,dict(status="PASS_WITH_SCOPE" if all(checks.values()) else "FAIL",checks=checks,
                    basis=str(path.relative_to(ROOT)),basis_sha256=sha(path),geometry_hash=config["acoustics"]["geometry_hash"],
                    acoustic_identity=db.metadata["acoustic_identity"],pressure_superposition_relative_L2=p_error,
                    gradient_superposition_relative_L2=g_error,
                    scope="Linear homogeneous water acoustics, 4 ideal disc sources; independent direct check, no skull or electromechanical calibration"))


if __name__ == "__main__":
    main()
