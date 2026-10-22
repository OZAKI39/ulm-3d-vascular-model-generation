"""Field-only tests and frozen historical-state audit; no bubble integration."""
import copy
import time
import numpy as np
import h5py
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from mevo_mb.control.scheduled_pzt import Schedule, Segment, one_hot, axis_pair, ORDER
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.fast_sim.nearfield import pair_reference
from dataclasses import replace
from .common import BASE, OLD, Q, CONFIG, ROOT, context, original, identity, read, write, sha


def primary(c, event, x, basis=None):
    basis = c["basis"] if basis is None else basis
    pb, gb = basis.query_basis(x)
    u = basis.vector(event["amplitude_Pa"], event["phase_rad"], event["enabled"])
    pressure = u @ pb
    gradient = np.einsum("i,inj->nj", u, gb)
    response = c["response"].query(pressure)  # strict actual pressure admission
    force = -.5 * np.real(response[:, 0, None] * gradient.conj())
    return pressure, gradient, force, response


def main():
    begin = time.perf_counter()
    assert (BASE/"PREREGISTRATION.json").exists()
    out = BASE/"diagnostics"; out.mkdir(exist_ok=False)
    c=context(); p=read(CONFIG); center=c["x0"].mean(0)[None]
    pb, gb = c["discs"].evaluate(center)
    cached_p, cached_g = c["basis"].query_basis(center)
    rows=[]
    for i,name in enumerate(ORDER):
        e=one_hot(name).event(); press,grad,force,_=primary(c,e,center)
        rows.append(dict(channel=name, source_Pa_equivalent=1500,
                         pressure_complex_Pa=press[0], pressure_magnitude_Pa=float(abs(press[0])),
                         gradient_complex_Pa_m=grad[0], primary_force_N=force[0],
                         force_direction=force[0]/np.linalg.norm(force[0]),
                         analytic_pressure_complex_Pa=pb[i,0]*1500,
                         analytic_gradient_complex_Pa_m=gb[i,0]*1500))
    lam=1500/2.5e6; a=.0035; zero=(a*a-(2*lam)**2)/(4*lam)
    ratio=float(abs(pb[0,0])/abs(pb[2,0]))
    specific=bool(ratio<.02 and abs(zero-.0045)<1e-5)
    # On-axis FNM reduction, and finite-difference derivative of the same field.
    axis_reference=[]
    for element in c["discs"].elements:
        z=float(np.linalg.norm(element["center_local_m"]))
        axis_reference.append(np.exp(-1j*c["discs"].k*z)-np.exp(-1j*c["discs"].k*np.sqrt(z*z+a*a)))
    axis_error=float(abs(pb[:,0]-axis_reference).max())
    fd=[]
    for axis in range(3):
        dx=np.zeros(3); dx[axis]=1e-8
        fd.append((c["discs"].evaluate(center+dx)[0]-c["discs"].evaluate(center-dx)[0])/(2e-8))
    fd=np.stack(fd,axis=-1)
    fd_error=float(abs(fd-gb).max()/abs(gb).max())
    assert axis_error<1e-12 and fd_error<1e-5
    write(out/"SINGLE_PZT.json",dict(identity=identity(), frequency_Hz=2.5e6, center_m=center[0],
          rows=rows, axis_nearfield_zero_m=zero, center_distance_Y_m=.0045,
          amplitude_efficiency_Y_over_X=ratio, intensity_efficiency_ratio=ratio**2,
          Y_scope="HOMOGENEOUS_WATER_PROXY_SPECIFIC" if specific else "NEARFIELD_HYPOTHESIS_NOT_CONFIRMED",
          analytic_onaxis_pressure_error=axis_error, analytic_gradient_FD_relative=fd_error,
          cache_center_pressure_error=float(abs(cached_p-pb).max()),
          cache_center_gradient_relative_error=float(abs(cached_g-gb).max()/abs(gb).max()),
          response_intervals_Pa=c["response"].metadata["intervals_Pa"],
          note="Computed water-proxy field and single-bubble primary force; no particle trajectory. No PDMS inference."))
    phase=[]; f0=primary(c,original("X"),center)[2][0,0]
    for shift in p["phase_perturbation_deg"]:
        e=copy.deepcopy(original("X"));e["phase_rad"][3]+=np.deg2rad(shift)
        pr,g,f,_=primary(c,e,center)
        phase.append(dict(delta_Xplus_deg=shift, pressure_Pa=float(abs(pr[0])), primary_force_N=f[0],
                          Fx_ratio=float(f[0,0]/f0), predicted_F10Q_scaled_endpoint_um=.200896324452*f[0,0]/f0,
                          status="FIELD_ONLY_NOT_TRAJECTORY"))
    slope=(next(v["Fx_ratio"] for v in phase if v["delta_Xplus_deg"]==1)-next(v["Fx_ratio"] for v in phase if v["delta_Xplus_deg"]==-1))/2
    write(out/"X_PHASE_SENSITIVITY.json",dict(rows=phase, center_Fx0_N=float(f0), fractional_slope_per_degree=slope,
          prediction=.06, hypothesis_relative_error_to_abs_prediction=float(abs(abs(slope)-.06)/.06)))
    # Independent field-only axial scan, extending beyond the small motion cache.
    peaks={}; arrays={}
    grid=np.linspace(-p["pair_peak_scan_halfwidth_m"],p["pair_peak_scan_halfwidth_m"],2401)
    for axis,ids in [("X",[2,3]),("Y",[0,1])]:
        points=np.zeros((len(grid),3)); ai="XYZ".index(axis);points[:,ai]=grid
        bp,_=c["discs"].evaluate(points)
        spline=CubicSpline(grid,bp[ids],axis=1)
        byphase={}; previous=0.0
        sequence=[0]+list(range(10,61,10))+list(range(-10,-61,-10))
        for degree in sequence:
            if degree==-10: previous=byphase[0]["position_m"]
            u=np.array([1,np.exp(1j*np.deg2rad(degree))])
            def derivative(x):
                return float(2*np.real(np.conj(u@spline(x))*(u@spline(x,1))))
            d=2*np.real(np.conj(u@bp[ids])*(u@spline(grid,1)))
            candidates=[brentq(derivative,grid[j],grid[j+1],xtol=1e-15) for j in range(len(grid)-1) if d[j]>0 and d[j+1]<=0]
            assert candidates, (axis,degree)
            peak=min(candidates,key=lambda x:(abs(x-previous),-x));previous=peak
            at=np.zeros((1,3));at[0,ai]=peak
            actual,actualg=c["discs"].evaluate(at)
            directP=u@actual[ids,0]; directG=u@actualg[ids,0,ai]
            residual=float(abs(2*np.real(directP.conj()*directG)))
            byphase[degree]=dict(position_m=peak, normalized_pressure_max=float(abs(directP)),
                                direct_intensity_derivative=residual,
                                all_local_maxima_m=candidates)
        records=[dict(delta_plus_minus_deg=degree,**byphase[degree]) for degree in p["pair_peak_phases_deg"]]
        fit=np.polyfit([v["delta_plus_minus_deg"] for v in records],[v["position_m"]*1e6 for v in records],1)
        peaks[axis]=dict(rows=records, linear_fit_um_per_degree=float(fit[0]),
                        scope="Continuous local pressure-maximum branch near center; scan +/-300um; no long-time particle simulation; not demonstrated force-equilibrium trap",
                        Y_scope="HOMOGENEOUS_WATER_PROXY_SPECIFIC" if axis=="Y" and specific else None)
        arrays[axis+"_axis_m"]=grid
        for degree in [-60,0,60]:arrays[axis+f"_Pabs_delta{degree}"]=abs(np.array([1,np.exp(1j*np.deg2rad(degree))])@bp[ids])
    np.savez_compressed(out/"AXIS_PAIR_FIELDS.npz",**arrays)
    write(out/"AXIS_PAIR_PRESSURE_MAXIMA.json",peaks)
    # Bounded proposal table. Initial COM-force prediction is never navigation evidence.
    trials=[]; targetforce=6*np.pi*.001*1.78e-6*2e-7/.002
    for amp in p["pair_static_search"]["amplitude_Pa"]:
        for degree in p["pair_static_search"]["x_delta_deg"]:
            s=axis_pair(x_amplitude=amp,x_delta_deg=degree)
            pr,g,f,v=primary(c,s.event(),c["x0"])
            row=dict(amplitude_Pa=amp, delta_deg=degree, pressure_max_Pa=float(abs(pr).max()), mean_primary_N=f.mean(0), schedule=Schedule((s,)).record())
            row["eligible"]=bool(row["pressure_max_Pa"]<=1200 and f.mean(0)[0]>0)
            row["force_error_N"]=float(abs(f.mean(0)[0]-targetforce))
            trials.append(row)
    candidates=[v for v in trials if v["eligible"]]
    assert candidates, "No admissible pair candidate; preserve field table and stop selection"
    best=min(candidates,key=lambda v:v["force_error_N"])
    controls={"X_PAIR":best["schedule"]}
    yone=max(rows[:2],key=lambda v:v["primary_force_N"][1])["channel"]
    controls["Y_ONE_HOT"]=Schedule((one_hot(yone),)).record()
    ypair=[]
    for degree in [-90,90]:
        s=axis_pair(y_amplitude=1500,y_delta_deg=degree,x_enabled=False,y_enabled=True)
        pr,g,f,v=primary(c,s.event(),c["x0"]);ypair.append((f.mean(0)[1],s))
    controls["Y_PAIR_90"]=Schedule((max(ypair,key=lambda v:v[0])[1],)).record()
    for key,scale in [("Y_FOUR_X025",.25),("Y_FOUR_X050",.5)]:
        e=copy.deepcopy(original("Y"));e["amplitude_Pa"][2]*=scale;e["amplitude_Pa"][3]*=scale
        controls[key]=Schedule((Segment.from_legacy(e),)).record()
    crows={}
    for name,schedule in controls.items():
        pr,g,f,v=primary(c,Schedule.from_record(schedule).segments[0].event(),c["x0"])
        crows[name]=dict(pressure_min_max_Pa=[float(abs(pr).min()),float(abs(pr).max())],mean_primary_N=f.mean(0),
                         COM_free_estimate_m=f.mean(0)*.002/(6*np.pi*.001*1.78e-6))
    write(BASE/"FIELD_SELECTED_SCHEDULES.json",controls)
    write(out/"PROPOSALS.json",dict(pair_trials=trials, selected_pair=best, candidate_fields=crows,
                                   no_trajectory_or_target_claim=True))
    # Historical Y endpoint sensitivity at common frozen positions separates field and geometry effects.
    fine=CachedWaterBasis(OLD/"cache/fine.h5"); hist={}; final={}
    physics=c["physics"]
    # Exactly the existing F2Simulator fallback susceptibility; never extrapolate.
    chi=c["response"].query(np.array([c["response"].knots[0]+0j]))[0,0]/c["response"].knots[0]
    for suffix in ["DT1","DT05","FINE"]:
        folder=Q/f"b_navigation_v1/runs/NAV_Y_{suffix}"
        with h5py.File(folder/"trajectory.h5") as h:
            final[suffix]=h["positions_m"][-1]
            sample={k:h["force_samples/"+k][-1] for k in h["force_samples"]}
            hist[suffix]=dict(last_force_sample=sample, last_step_contact_pair_count=int(np.count_nonzero(h["pair_force_N"][-1]>0)),
                              step_pressure_end_Pa=h["extra_metrics"][-1,3:5], time_s=float(h["time_s"][-1]))
    frozen={}
    for name,x in final.items():
        frozen[name]={}
        for bname,b in [("coarse",c["basis"]),("fine",fine)]:
            pr,g,f,v=primary(c,original("Y"),x,b)
            sec=pair_reference(x,v,pr,physics,chi)
            far=pair_reference(x,v,pr,replace(physics,force_model="f1_smooth_near_field"),chi)
            frozen[name][bname]=dict(pressure_Pa=pr,primary_N=f,secondary_N=sec,F2_increment_N=sec-far)
    write(out/"HISTORICAL_Y_FROZEN_AUDIT.json",dict(scope="Endpoint geometry versus cache sensitivity, plus last saved force sample with its true time; no causal exclusivity assumed", runs=hist, frozen=frozen,
          Y_scope="HOMOGENEOUS_WATER_PROXY_SPECIFIC" if specific else "UNCONFIRMED"))
    write(out/"EXECUTION.json",dict(identity=identity(),wall_s=time.perf_counter()-begin, particle_trajectories_run=0,
          input_sha256=read(BASE/"PREREGISTRATION.json")["input_sha256"],Y_scope_specific=specific))
    print("Y/X amplitude efficiency",ratio,"zero mm",zero*1e3,"X fractional phase slope",slope)
    print("Pair slopes um/deg",{k:v["linear_fit_um_per_degree"] for k,v in peaks.items()})
    print("Selected X pair",best["amplitude_Pa"],best["delta_deg"],"predicted force",best["mean_primary_N"])


if __name__ == "__main__":
    main()
