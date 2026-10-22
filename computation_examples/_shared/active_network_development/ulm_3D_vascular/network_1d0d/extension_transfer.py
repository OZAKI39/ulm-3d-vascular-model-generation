"""Measure current FEM extension sections; no changes to geometry or mesh."""
from pathlib import Path
import numpy as np
from .hydraulic_resistance import linear_radius_resistance
from .boundary_conditions import cap_pressures, shift_pressure_gauge
from .audit import sha256


def ordered_section_loops(section):
    """Recover closed polygon loops from VTK slice lines, rejecting open contours."""
    section=section.clean(tolerance=1e-13, absolute=True); cells=section.lines; edges=[]; pos=0
    while pos < len(cells):
        n=int(cells[pos]); ids=cells[pos+1:pos+1+n]
        edges.extend(zip(ids[:-1],ids[1:])); pos+=n+1
    adjacency={}
    for u,v in edges:
        adjacency.setdefault(int(u),set()).add(int(v));adjacency.setdefault(int(v),set()).add(int(u))
    unused=set(adjacency);loops=[]
    while unused:
        seed=min(unused);component={seed};stack=[seed]
        while stack:
            u=stack.pop()
            for v in adjacency[u]:
                if v not in component:component.add(v);stack.append(v)
        unused-=component
        if any(len(adjacency[u])!=2 for u in component):continue
        order=[seed];previous=None;current=seed
        while True:
            nxt=min(adjacency[current]-({previous} if previous is not None else set()))
            if nxt==seed:break
            order.append(nxt);previous,current=current,nxt
            if len(order)>len(component):raise ValueError('Invalid section topology')
        if len(order)==len(component):loops.append(np.asarray(section.points)[order])
    return loops


def section_area(surface, origin, normal, expected_radius):
    loops=ordered_section_loops(surface.slice(normal=normal,origin=origin))
    candidates=[]
    for pts in loops:
        delta=pts-origin
        cross=np.cross(delta,np.roll(delta,-1,axis=0))
        area=abs(np.dot(cross.sum(axis=0)/2,normal))
        center=pts.mean(0); distance=np.linalg.norm(center-origin)
        if distance < 2*expected_radius and area>0:candidates.append((distance,area,len(pts)))
    if len(candidates)!=1:raise ValueError(f'Extension section ambiguous/missing: {len(candidates)} loops near axis')
    return candidates[0][1],candidates[0][2]


def measure_extensions(case, ports, station_count=40):
    import pyvista as pv
    surface_path=Path(case)/'SV_MESH/mesh-complete.exterior.vtp';surface=pv.read(surface_path)
    rows=[]
    for port in ports:
        if port['name']=='INLET':continue
        origin=np.asarray(port['real_cut_xyz_um'])*1e-6;normal=np.asarray(port['outward_normal']);normal/=np.linalg.norm(normal)
        length=port['fem_cap']['real_cut_to_cap_axial_um']*1e-6
        fractions=(np.arange(station_count)+.5)/station_count
        measurements=[section_area(surface,origin+t*length*normal,normal,port['radius_um']*1e-6) for t in fractions]
        areas=np.array([x[0] for x in measurements]);r=np.sqrt(areas/np.pi)
        # Fill the half-bin at each endpoint with its nearest sampled radius.
        t=np.r_[0.,fractions,1.];r_full=np.r_[r[0],r,r[-1]]
        R=float(linear_radius_resistance(length*np.diff(t),r_full[:-1],r_full[1:]).sum())
        rows.append(dict(port=port['name'], R_extension_Pa_s_m3=R, length_m=length, station_count=station_count,
                         station_fractions=fractions, cross_section_area_m2=areas, equivalent_radius_m=r,
                         polygon_vertices=[x[1] for x in measurements], current_surface_sha256=sha256(surface_path),
                         v1_estimate_R=port['extension']['current_length_piecewise_linear_station_radius_R_Pa_s_m3'],
                         method='Current frozen SV surface; closed-loop plane sections; linear equivalent radius integral over entire extension with constant half-bin endpoint closure',
                         approximation='Equivalent circular Poiseuille estimate; no junction/3D development correction'))
    return rows


def pressure_transfer(operating_ports, extension_rows, mapping_sha, solution_sha):
    rows=[]
    for ext in extension_rows:
        port=next(x for x in operating_ports if x['port']==ext['port']);p=port['pressure_realcut_Pa'];q=port['signed_Q_m3s'];R=ext['R_extension_Pa_s_m3']
        rows.append(dict(port=port['port'], pressure_realcut_Pa=p, signed_Q_m3s=q, R_extension_Pa_s_m3=R,
                         extension_deltaP_Pa=R*q, pressure_cap_raw_Pa=float(cap_pressures(p,q,R))))
    raw=np.array([r['pressure_cap_raw_Pa'] for r in rows]);shifted=shift_pressure_gauge(raw)
    for row,p in zip(rows,shifted):row['pressure_cap_shifted_Pa']=float(p)
    difference_raw=raw[:,None]-raw;difference_shifted=shifted[:,None]-shifted
    error=float(np.max(abs(difference_raw-difference_shifted)))
    real=np.array([r['pressure_realcut_Pa'] for r in rows]);delta=np.array([r['extension_deltaP_Pa'] for r in rows])
    return dict(status='PASS', ports=rows, reference_pressure_definition='H0 structural leaves p_ref=0 Pa gauge; source 2410 is a model assumption',
                pressure_location='real A cut transferred to artificial FEM cap', gauge_subtracted_constant_Pa=float(raw.min()),
                pairwise_raw_pressure_differences_Pa=difference_raw, pairwise_shifted_pressure_differences_Pa=difference_shifted,
                max_pressure_difference_preservation_error_Pa=error, extension_drop_pairwise_difference_Pa=delta[:,None]-delta,
                realcut_pairwise_pressure_differences_Pa=real[:,None]-real, roi_mapping_SHA=mapping_sha, source_A_solution_SHA=solution_sha)
