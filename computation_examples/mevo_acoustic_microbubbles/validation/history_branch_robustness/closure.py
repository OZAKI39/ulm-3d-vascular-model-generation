"""Offline, independently checked local periodic mobility for the history reference."""
import json
from pathlib import Path
import numpy as np
from mevo_mb.bubbles.adaptive_table import ResponseNodeCache
from mevo_mb.bubbles.apecss import RadialCase
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution
from validation.history_branch_robustness.periodic import periodic_translation


def build_closure(out,root):
    out,root = Path(out),Path(root)
    parameters = SonoVueParameterDistribution(root / "configs/stage2/sonovue_parameters.json").parameter_set(1.78e-6)
    cache = ResponseNodeCache(root / "results/stage2_5/response_table/node_cache",root / "third_party/build/apecss/ultrasound_apecss")
    nodes = {}

    def evaluate(pressure):
        key = float(pressure)
        if key in nodes:
            return nodes[key]
        case = RadialCase(radius_m=1.78e-6,frequency_Hz=2.5e6,pressure_Pa=key,
                         elasticity_N_m=parameters.shell_elasticity_N_m.value,
                         shell_viscosity_kg_s=parameters.shell_viscosity_kg_s.value,sigma_initial_N_m=.036,
                         provenance="Stage4 validation-only periodic history closure; frozen Tu2009/APECSS parameter manifold")
        node = cache.get(case)
        if node["status"] != "COMPLETE" or not node["metrics"]["valid"] or not node["metrics"]["harmonic_energy_resolved"]:
            raise ValueError("Invalid radial node for history closure")
        with np.load(cache.root / node["key"] / node["payload"]) as data:
            radius = data["radius_m"][-case.samples_per_cycle:]
        values,changes = {},{}
        for label,enabled in [("B",False),("C",True)]:
            coarse = periodic_translation(radius,case.frequency_Hz,history=enabled,n=129)
            fine = periodic_translation(radius,case.frequency_Hz,history=enabled,n=257)
            values[label] = fine["mean_velocity_per_gradient"].tolist()
            changes[label] = float(np.linalg.norm(fine["mean_velocity_per_gradient"]-coarse["mean_velocity_per_gradient"])
                                    /max(np.linalg.norm(fine["mean_velocity_per_gradient"]),1e-30))
            if changes[label]>.02 or fine["relative_linear_residual"]>1e-9:
                raise ValueError("Periodic history reference failed spectral refinement/residual gate")
        fundamental = complex(node["harmonics_real"][0],node["harmonics_imag"][0])
        expected_force = -.5*np.array([fundamental.real,fundamental.imag])
        force_phase_error = float(np.linalg.norm(fine["mean_force_per_gradient"]-expected_force)/max(np.linalg.norm(expected_force),1e-30))
        if force_phase_error > .02:
            raise ValueError("Periodic radial/gradient phase does not match the existing APECSS primary-force convention")
        result = dict(pressure_Pa=key,coefficients=values,resolution_changes=changes,cache_key=node["key"],
                      primary_force_phase_convention_error=force_phase_error,
                      radial_metrics=node["metrics"],radius_raw_payload_sha256=node["payload_sha256"])
        nodes[key] = result
        return result

    knots = [0.1,1.,10.,100.,500.,1000.,2500.,5000.,10000.,15000.,20000.,25000.,30000.,32500.,35000.,37500.,40000.]
    for p in knots:
        evaluate(p)
    # Deterministic midpoint refinement; holdouts below use different locations.
    for _ in range(8):
        add = []
        for lower,upper in zip(knots[:-1],knots[1:]):
            midpoint = (lower+upper)/2
            actual = evaluate(midpoint)
            errors = []
            for label in ["B","C"]:
                predicted = (np.array(nodes[lower]["coefficients"][label])+nodes[upper]["coefficients"][label])/2
                exact = np.array(actual["coefficients"][label])
                errors.append(np.linalg.norm(predicted-exact)/max(np.linalg.norm(exact),1e-30))
            if max(errors)>.01:
                add.append(midpoint)
        if not add:
            break
        knots = sorted(knots+add)
    else:
        raise ValueError("Periodic closure interpolation did not converge")
    coefficients = {label:np.array([nodes[p]["coefficients"][label] for p in knots]) for label in ["B","C"]}
    holdouts = []
    rng = np.random.default_rng(20261008)
    for lower,upper in zip(knots[:-1],knots[1:]):
        pressure = float(rng.uniform(lower,upper))
        actual = evaluate(pressure)
        errors = {}
        for label in ["B","C"]:
            predicted = np.array([np.interp(pressure,knots,coefficients[label][:,j]) for j in range(2)])
            exact = np.array(actual["coefficients"][label])
            errors[label] = float(np.linalg.norm(predicted-exact)/max(np.linalg.norm(exact),1e-30))
        holdouts.append(dict(pressure_Pa=pressure,relative_coefficient_error=errors))
    passed = all(max(row["relative_coefficient_error"].values())<.02 for row in holdouts)
    result = dict(status="PASS_WITH_SCOPE" if passed else "FAIL",pressure_knots_Pa=knots,coefficients={k:v.tolist() for k,v in coefficients.items()},
                  nodes=list(nodes.values()),holdouts=holdouts,spectral_sizes=[129,257],
                  scope="Adiabatic infinitely periodic local variable-radius reference; not experimental calibration or fully resolved moving history")
    (out / "PERIODIC_CLOSURE.json").write_text(json.dumps(result,indent=2)+"\n")
    if not passed:
        raise ValueError("Independent periodic closure holdout failure")
    return result


def mean_slip(closure,label,pressure,gradient):
    p = np.asarray(pressure)
    amp = abs(p)
    knots = closure["pressure_knots_Pa"]
    if np.any(amp<knots[0]) or np.any(amp>knots[-1]):
        raise ValueError("History closure out of bounds; no pressure clipping/extrapolation")
    coefficients = np.array(closure["coefficients"][label])
    a,b = [np.interp(amp,knots,coefficients[:,j]) for j in range(2)]
    relative_gradient = gradient*np.exp(-1j*np.angle(p))[:,None]
    return a[:,None]*relative_gradient.real+b[:,None]*relative_gradient.imag
