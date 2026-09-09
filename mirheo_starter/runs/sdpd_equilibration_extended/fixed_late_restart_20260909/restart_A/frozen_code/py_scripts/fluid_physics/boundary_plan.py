"""Finite per-port design with pinned geometry and honest native API gaps."""
import numpy as np


def build_boundary_plan(case, c, units, eos=None):
    rc=units.L0*c['space']['rc_star'];nphys=c['space']['n_star']/units.L0**3
    entries=[]
    for port in case['ports']:
        patch=port['patch'];cap=np.array(patch['center_m']);normal=np.array(patch['outward_normal'])
        measurement=port['measurement_planes']['central'];distance=float(np.dot(cap-np.array(measurement['origin_m']),normal))
        inlet=patch['kind']=='inlet';Q=case['target_volume_flow_m3_s']
        entry={'name':port['name'],'port_id':patch['port_id'],'entity_id':patch['entity_id'],
               'original_role':patch['original_role'],'boundary_origin':patch['boundary_origin'],
               'physical_requirement':({'total_Q_m3_s':Q,'total_mass_flow_kg_s':case['target_mass_flow_kg_s'],'profile':case['inlet_profile']} if inlet else {'gauge_pressure_pa':case['outlet_gauge_pressures_pa'][port['name']]}),
               'requirement_source':case['accepted_report']+'#/runtime_contract','data_origin':case['data_origin'],
               'cap_center_m':cap.tolist(),'outward_normal':normal.tolist(),'cap_area_m2':patch['area_m2'],
               'original_physical_cut_reference':port.get('original_physical_cut_reference'),
               'flux_sign':port['flux_sign'],'measurement_plane':measurement,
               'check_planes':{k:v for k,v in port['measurement_planes'].items() if k!='central'},
               'measurement_to_cap_axial_distance_m':distance,'existing_extension':port['existing_extension'],
               'apply_at':'existing distal cap; source patch face_ids identify finite cap, not an infinite plane',
               'finite_control_region':{'category':'DESIGN_CHOICE','origin_m':cap.tolist(),'outward_axis':normal.tolist(),
                   'axial_interval_m':[-rc,0.0],'transverse_mask':'finite cap polygon projected along normal, intersect local connected lumen; use per-port ownership and disjoint masks',
                   'measured_polygon_reference':port['measurement_source'],
                   'warning':'Control slab is a design, no new volume/particles/SDF generated; its effect at retained measurement planes needs convergence test.',
                   'clearance_between_measurement_and_control_m':distance-rc},
               'existing_buffer_action':'REUSE existing extension, no second extension; pressure remains prescribed at legacy distal cap',
               'state':'REQUIRES_EXTENSION','true_vessel_boundary_implemented':False}
        if inlet:
            entry.update(mean_inward_normal_velocity_m_s=Q/patch['area_m2'],
                         mean_inward_velocity_star=units.to_star(Q/patch['area_m2'],'velocity'),
                         expected_particle_supply_per_s=Q*nphys,
                         insertion={'number_density_star':c['space']['n_star'],'particle_mass_star':c['space']['m_star'],
                                    'kBT_star':c['space']['kBT_star'],'temperature_K':case['temperature_K'],
                                    'thermal_velocity_variance_star':c['space']['kBT_star']/c['space']['m_star'],
                                    'mean_velocity_direction':(-normal).tolist(),
                                    'profile_choice':'UNVERIFIED: flat normal cap seed is a possible initial choice, Q/A alone does not specify final profile',
                                    'mass_accumulator':'expected additions n_phys*Q_command*dt_phys; fractional stochastic accumulator, record actual inserted mass'},
                         existing_apis=[{'api':'Plugins.createVelocityInlet(name,pv,implicit_surface_func,velocity_field,resolution,number_density,kBT)',
                                         'scope':'static finite-mask velocity supply prototype; geometry and velocity sampled at setup, no verified runtime setter; native rate uses abs(nA dot v), requires correctly inward velocity and zero velocity on unrelated surface; <5 additions/triangle/step asserted in velocity_inlet.cu'},
                                        {'api':'Plugins.createVelocityControl','scope':'box average velocity PID only; cannot equate it with arbitrary-port total flux control'}],
                         extension={'inputs':['finite cap patch and disjoint local mask','frozen central and check polygons','Q_target, rho, kBT, m, rc, dt','controller gains and update/block intervals to calibrate'],
                                    'measure':['signed crossing mass and aperture-integrated velocity Q at original measurement plane','local density/temperature','supplied minus removed mass and total inventory'],
                                    'actions':'Regulate supply rate or local reservoir forcing with a bounded PI controller. Preserve Maxwell fluctuations. Do not remove reverse movers to satisfy a gate.',
                                    'feedback':'error = Q_target - time-averaged measured Qin; distinguish command from measured flux; anti-windup and explicit saturation/inventory failure diagnostics'},
                         next_validation=['finite-mask insertion rate, thermal variance and global mass balance in synthetic box',
                                          'steady synthetic tube: measured Q and profile, controller interval/gain sensitivity, no-slip and leakage',
                                          'later vessel: frozen three planes, short/long physical windows, inherited 1% Qin/mass/steady gates'])
        else:
            entry.update(pressure_mapping={'formula':'p_target_star = p_reference_star + (p_gauge_j - 0 Pa gauge)/P0',
                                          'gauge_reference_pa':0.0,'reference_category':'DESIGN_CHOICE',
                                          'legacy_offset_pa':case['pressure_reference']['legacy_numerical_offset_pa'],
                                          'legacy_offset_use':'provenance only, never applied as thermodynamic absolute pressure',
                                          'P0_pa':units.scales['pressure'],'measured_EOS':eos if eos else {'status':'UNVERIFIED'}},
                         existing_apis=[{'api':'Plugins.createDensityOutlet(name,pvs,number_density,region,resolution)','scope':'negative-inside bounded region, removes excess particles only; no reservoir supply or complete pressure BC'},
                                        {'api':'Plugins.createDensityControl','scope':'native smooth level-set density PI(D) forcing; reusable controller component with calibrated EOS, not a full pressure/reservoir boundary'},
                                        {'api':'Plugins.createVirialPressurePlugin','scope':'local virial sum only, positive-inside predicate; add local kinetic stress, actual V, flow subtraction and pair-boundary treatment'},
                                        {'api':'Plugins.createPlaneOutlet','scope':'not selected: infinite half-space sink, neither finite branch isolation nor prescribed pressure'}],
                         extension={'inputs':['finite port/slab and local connected-volume masks','pressure reference and measured p(n,T) with uncertainty and domain limits','target gauge pressure, reservoir number density, kBT, particle mass','controller gains, saturation, density bounds, signed measurement polygons'],
                                    'measure':['time-averaged normal traction including kinetic and virial stress in reservoir/control slab','number density and temperature','signed Qin/Qout, insertion/removal inventory'],
                                    'actions':'Two-way local reservoir exchange plus smooth normal force/density feedback, Maxwell insertion and unbiased spatial removal as required for reservoir inventory. Never condition removal on velocity sign.',
                                    'limits':'No extrapolation beyond measured EOS; reject unattainable target or density ratio outside 0.9..1.1. Do not change prescribed gauge pressure.',
                                    'pressure_location':'Maintain control reference at existing distal cap. If external reservoir later needed, separately quantify added hydraulic resistance and match cap pressure.'},
                         next_validation=['isolated periodic slab reservoir: pressure/density command, thermal stress and signed mass closure',
                                          'two disconnected synthetic branches: no cross-branch deletion/forcing',
                                          'small pressure-driven tube: pressure difference and Q, cap-to-measurement resistance, buffer length sensitivity',
                                          'later multi-outlet vessel: short AND long mean Qout >= -0.05*abs(Qin); report violations without clipping'])
        entries.append(entry)
    return {'status':'REQUIRES_EXTENSION','plan_complete':True,'ports':entries,
            'wall_plan':{'status':'PLANNED_WITH_EXISTING_APIS','requirements':case['walls'],
                         'reuse':'Native stationary wall signed-distance confinement/bounce-back and frozen wall particles with DPD interaction (local tests/doc_scripts/walls.py).',
                         'next_validation':'Synthetic channel before real SDF: no-slip profile/slip length, leakage, wall density and rc/SDF-spacing convergence. No claim of validation from periodic fluid tests.'},
            'cpp_cuda_modified':False,'native_api_source':'vendor/Mirheo/src/mirheo/bindings/plugins.cpp:196,219,235,751,771,788',
            'required_next_plugins':['bounded arbitrary-polygon signed-flux and local traction measurement','runtime adjustable finite-cap flow/reservoir supply','two-way local pressure reservoir combining EOS, density force controller and exchange'],
            'current_scope':'Implementation design only. Native components read from pinned source; none of the real vessel boundaries executed.'}
