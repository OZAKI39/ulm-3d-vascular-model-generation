"""Small field-only control screen; never optimize on particle outcomes."""
import itertools
import numpy as np
import yaml
from scipy.interpolate import RegularGridInterpolator
from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
from mevo_mb.bubbles.forces import primary_bjerknes
from mevo_mb.control.pzt import PZTArray
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4.acoustic import configuration
from scripts.stage4.common import ROOT, new_run, require, sha, write, finish


def main():
    source = require("acoustic")
    config = configuration()
    geometry = BifurcationGeometry(yaml.safe_load((ROOT / "configs/stage4/bifurcation.yaml").read_text())["geometry"])
    table = config["bubble_response"]["path"]
    extension = ROOT / "results/stage4/acoustic/response_extension/STATE.json"
    if extension.exists():
        table = require("acoustic/response_extension")["table"]
    model = AdaptiveResponseTableBubbleModel(ROOT / table)
    bounds = np.asarray(model.metadata["bounds"])
    out = new_run("phase_screening",dict(config=config,response_table=table,
          phase_candidates_rad=[0.,np.pi/2,np.pi,3*np.pi/2],source_amplitude_ratios=[1.,2.,4.],
          protocol="No trajectories screened; maximize signed y-primary-force ROI average within validated table domain; mirror by swapping y sources",
          formal_eligibility="All four elements active and nonzero relative phase between y sources; single-source amplitude-only result retained as diagnostic",
          maximum_ROI_pressure_Pa=39000.),["acoustic"])
    db = AcousticBasisDatabase.load_basis(ROOT / source["basis"])
    assert db.metadata["geometry_hash"] == geometry.identity
    xyz = np.stack(np.meshgrid(np.linspace(-.0018,.001,15),np.linspace(-.001,.001,13),
                               np.linspace(-.0008,.0008,7),indexing="ij"),axis=-1).reshape(-1,3)
    xyz = xyz[geometry.contains(xyz,1.78e-6)]
    h = np.array([RegularGridInterpolator(db.axes_m,p)(xyz) for p in db.pressure_basis])
    dh = np.array([RegularGridInterpolator(db.axes_m,g)(xyz) for g in db.gradient_basis])
    phases = np.arange(4)*np.pi/2
    patterns = []
    for ratio,p0,pz in itertools.product([1.,2.,4.],phases,phases):
        patterns.append(np.array([np.exp(1j*p0),ratio,np.exp(1j*pz),np.exp(1j*pz)]))
    patterns.append(np.array([0.,1.,0.,0.],complex))
    mirror = np.array([1,0,2,3])
    rows,controls = [],[]
    for index,u0 in enumerate(patterns):
        pm = np.einsum("i,ixyz->xyz",u0,db.pressure_basis)
        pmm = np.einsum("i,ixyz->xyz",u0[mirror],db.pressure_basis)
        # Bound the whole interpolation ROI, including all cell corners, at 39 kPa.
        amplitude = 39000./max(float(abs(pm).max()),float(abs(pmm).max()))
        u = u0*amplitude
        p = np.einsum("i,in->n",u,h)
        grad = np.einsum("i,inj->nj",u,dh)
        valid = (abs(p)>=bounds[2,0]) & (abs(p)<=bounds[2,1])
        force = primary_bjerknes(model.response(1.78e-6,2.5e6,p[valid],config["medium"]),grad[valid])
        # Invalid points are recorded and excluded from force evaluation, never extrapolated.
        score = float(force[:,1].sum()/len(xyz)) if len(force) else -1e99
        rows.append(dict(index=index,valid_force_point_fraction=float(valid.mean()),mean_transverse_force_N=score,
                         formal_eligible=bool(np.count_nonzero(u0)==4 and abs(np.angle(u0[0]*np.conj(u0[1])))>1e-10),
                         positive_transverse_fraction=float(np.mean(force[:,1]>0)) if len(force) else None,
                         sampled_pressure_range_Pa=[float(abs(p).min()),float(abs(p).max())],
                         whole_ROI_peak_pressure_Pa=39000.,amplitudes_Pa=abs(u).tolist(),phases_rad=np.angle(u).tolist()))
        controls.append(u)
    rank = sorted(rows,key=lambda r:(r["formal_eligible"],r["valid_force_point_fraction"]>=.75,r["mean_transverse_force_N"]),reverse=True)
    chosen = rank[0]
    if chosen["mean_transverse_force_N"]<=0:
        finish(out,dict(status="FAIL",candidates=rows,error="No candidate with positive mean y steering force"))
        raise RuntimeError("Phase field screen has no positive candidate")
    selected = controls[chosen["index"]]
    symmetric = np.ones(4,complex)
    symmetric *= 39000./float(np.abs(np.einsum("i,ixyz->xyz",symmetric,db.pressure_basis)).max())
    case_controls = dict(B0=np.zeros(4,complex),B1=symmetric,B2=selected,B3=selected[mirror])
    # Full lumen audit: identify response-domain shortfall before any formal trajectory.
    points = np.stack(np.meshgrid(*db.axes_m,indexing="ij"),axis=-1)
    lumen = geometry.contains(points.reshape(-1,3),-np.sqrt(3)*config["grid"]["spacing_m"][0]).reshape(points.shape[:-1])
    cases = {}
    for label,u in case_controls.items():
        array = PZTArray.from_config(config)
        array.set_amplitudes(abs(u))
        array.set_phases(np.angle(u))
        np.testing.assert_allclose(array.get_control_vector(),u,atol=1e-10,rtol=1e-14)
        field = db.combine(array.get_control_vector())
        path = out / (label+"_field.h5")
        field.save(path)
        pressure = abs(field.pressure_Pa[lumen])
        cases[label] = dict(field=str(path.relative_to(ROOT)),field_sha256=sha(path),
                           amplitudes_Pa=abs(u).tolist(),phases_rad=np.angle(u).tolist(),enabled=(abs(u)>0).tolist(),
                           lumen_pressure_range_Pa=[float(pressure.min()),float(pressure.max())],
                           lumen_fraction_below_table=float(np.mean(pressure<bounds[2,0])) if label!="B0" else 0.,
                           lumen_fraction_above_table=float(np.mean(pressure>bounds[2,1])) if label!="B0" else 0.)
    mirrored_xyz = xyz.copy()
    mirrored_xyz[:,1] *= -1
    a = db.combine(selected).interpolate_complex_pressure(xyz)
    b = db.combine(selected[mirror]).interpolate_complex_pressure(mirrored_xyz)
    mirror_error = float(np.linalg.norm(a-b)/np.linalg.norm(a))
    all_covered = all(v["lumen_fraction_below_table"]==0 and v["lumen_fraction_above_table"]==0 for v in cases.values())
    status = "PASS_WITH_SCOPE" if all_covered else "WARNING"
    write(out / "SCREEN_POINTS.npz.json",dict(count=len(xyz),definition="Fixed grid in distal parent and junction, finite-radius lumen filtered"))
    np.savez_compressed(out / "SCREEN_POINTS.npz",positions_m=xyz)
    finish(out,dict(status=status,candidates=rows,ranking=[r["index"] for r in rank],chosen=chosen,cases=cases,
                   response_table=table,response_table_sha256=sha(ROOT / table),basis_sha256=source["basis_sha256"],
                   mirror_pressure_relative_L2=mirror_error,geometry_hash=source["geometry_hash"],
                   scope="Field-only finite candidate screen; no global optimization or particle-outcome selection. Lumen guard audit includes adjacent interpolation corners.",
                   next_action="Ready for guarded trajectories" if all_covered else "Resolve measured table-domain shortfall with independently validated offline APECSS coverage; do not clip pressure"))


if __name__ == "__main__":
    main()
