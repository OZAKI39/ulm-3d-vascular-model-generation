"""One-way frozen design + extension GEOMETRY -> fixed FEM cap pressures.

No CFD fields, optimizer calls, or trajectory inputs. Legacy section geometry
and pressure/gauge helpers are reused with explicit configured viscosity.
"""
from pathlib import Path
import json
import numpy as np
from .audit import sha256
from .balance_design import load_frozen_design
from .hydraulic_resistance import linear_radius_resistance
from .extension_transfer import section_area
from .boundary_conditions import cap_pressures,shift_pressure_gauge

PORTS=('O1','O2','O3')


def measure_extension_geometry(case,ports,*,mu_pa_s,station_count=40):
    import pyvista as pv
    if not np.isfinite(mu_pa_s) or mu_pa_s<=0:raise ValueError('Explicit positive SI viscosity required')
    if not isinstance(station_count,int) or station_count<2:raise ValueError('At least two stations required')
    surface_path=Path(case)/'SV_MESH/mesh-complete.exterior.vtp'
    surface=pv.read(surface_path)
    rows=[]
    for name in PORTS:
        port=next(p for p in ports if p['name']==name)
        origin=np.asarray(port['real_cut_xyz_um'])*1e-6
        normal=np.asarray(port['outward_normal'],float);normal/=np.linalg.norm(normal)
        length=port['fem_cap']['real_cut_to_cap_axial_um']*1e-6
        if not np.isfinite(length) or length<=0:raise ValueError('Positive extension length required')
        fractions=(np.arange(station_count)+.5)/station_count
        measurements=[section_area(surface,origin+t*length*normal,normal,port['radius_um']*1e-6) for t in fractions]
        areas=np.array([m[0] for m in measurements]);radii=np.sqrt(areas/np.pi)
        stations=np.r_[0.,fractions,1.];rfull=np.r_[radii[0],radii,radii[-1]]
        resistance=float(linear_radius_resistance(length*np.diff(stations),rfull[:-1],rfull[1:],mu=mu_pa_s).sum())
        rows.append(dict(port=name,R_extension_Pa_s_m3=resistance,mu_pa_s=float(mu_pa_s),
            length_m=float(length),station_count=station_count,station_fractions=fractions,
            cross_section_area_m2=areas,equivalent_radius_m=radii,polygon_vertices=[m[1] for m in measurements],
            surface_sha256=sha256(surface_path),
            method='Closed polygon plane sections; analytical linear equivalent-radius integral; constant half-bin endpoints',
            limitation='Equivalent circular, steady Poiseuille resistance; no junction or 3D entrance correction'))
    return rows


def transfer_prediction(prediction,extensions):
    """Analytical transformation; also usable for the pure-0D H0 regression."""
    rows=[]
    if sorted(r['port'] for r in extensions)!=list(PORTS):raise ValueError('Exactly one extension for each outlet required')
    for name in PORTS:
        ext=next(e for e in extensions if e['port']==name)
        p=prediction['real_cut_pressure_pa'][name];q=prediction['port_flow_m3_s'][name]
        R=ext['R_extension_Pa_s_m3']
        rows.append(dict(port=name,pressure_realcut_Pa=p,signed_Q_m3s=q,R_extension_Pa_s_m3=R,
            extension_deltaP_Pa=R*q,pressure_cap_raw_Pa=float(cap_pressures(p,q,R))))
    raw=np.array([r['pressure_cap_raw_Pa'] for r in rows]);shifted=shift_pressure_gauge(raw)
    for r,p in zip(rows,shifted):r['pressure_cap_shifted_Pa']=float(p)
    before=raw[:,None]-raw;after=shifted[:,None]-shifted
    np.testing.assert_allclose(before,after,rtol=1e-12,atol=1e-12)
    return dict(status='PASS',ports=rows,gauge_subtracted_constant_Pa=float(raw.min()),
        gauge_added_constant_Pa=float(-raw.min()),pairwise_raw_pressure_differences_Pa=before,
        pairwise_shifted_pressure_differences_Pa=after,
        max_pressure_difference_preservation_error_Pa=float(abs(before-after).max()),
        reference_pressure_definition='Common distal gauge from FROZEN_DESIGN; uniform cap shift preserves all differences',
        pressure_location='Physical real cut -> artificial FEM extension cap; pressure traction, not outlet flow BC')


def export_frozen_design(frozen_path,case,ports_path):
    frozen=load_frozen_design(frozen_path)
    mapping=json.loads(Path(ports_path).read_text())
    if mapping['status']!='PASS':raise ValueError('Verified port mapping required')
    if sha256(ports_path)!=frozen['input_hashes']['port_mapping_sha256']:raise ValueError('Port identity changed')
    extensions=measure_extension_geometry(case,mapping['ports'],mu_pa_s=frozen['dynamic_viscosity_pa_s'])
    result=transfer_prediction(frozen['optimized']['prediction'],extensions)
    result.update(source_model=frozen['model_name'],source_parameter_status='FROZEN_DESIGN',
        purpose='BEST_FEASIBLE_BALANCE_FORWARD_PREDICTION',design_bounds=frozen['design_bounds'],
        parameters=frozen['optimized']['parameters'],fractions_0D=frozen['optimized']['prediction']['outlet_flow_fraction'],
        mu_pa_s=frozen['dynamic_viscosity_pa_s'],source_geometry_sha256=frozen['source_geometry_sha256'],
        source_frozen_design_sha256=sha256(frozen_path),source_frozen_content_sha256=frozen['artifact_sha256'],
        FEM_surface_sha256=extensions[0]['surface_sha256'],port_mapping_sha256=sha256(ports_path),
        extension_measurements=extensions,CFD_feedback_used=False,
        implementation_sha256=sha256(__file__))
    return result
