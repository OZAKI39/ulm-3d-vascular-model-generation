#!/usr/bin/python3
"""Freeze supported BC choices, derive units, and prepare audited flux quadrature.

Only writes RUN_DIR. Imports Step2's STL reader without executing its workflow.
No new port identification, lattice population reconstruction, or solver execution.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

HC = Path('/home/lzy/projects/hemocell_starter')
S2 = Path('/home/lzy/projects/compre_output/step2/20260912_225418')
sys.path.insert(0, str(HC/'test_code/step2_vascular_voxelization_smoke'))
from prepare_port_contract import read_stl, write_json, projected_inside

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def triangulate_loop(coords):
    points = vtk.vtkPoints(); points.SetDataTypeToDouble()
    for p in coords:
        points.InsertNextPoint(*p)
    lines = vtk.vtkCellArray()
    for i in range(len(coords)):
        lines.InsertNextCell(2); lines.InsertCellPoint(i); lines.InsertCellPoint((i+1)%len(coords))
    poly = vtk.vtkPolyData(); poly.SetPoints(points); poly.SetLines(lines)
    cutter = vtk.vtkContourTriangulator(); cutter.SetInputData(poly); cutter.Update()
    out = cutter.GetOutput()
    assert out.GetNumberOfPolys() == len(coords)-2, 'Flux contour triangulation is incomplete'
    pts = vtk_to_numpy(out.GetPoints().GetData())
    faces = vtk_to_numpy(out.GetPolys().GetData()).reshape(-1,4)
    assert np.all(faces[:,0]==3)
    return pts[faces[:,1:]]

def cut_loops(tri, center, normal):
    # Exact intersection with full frozen STL, then tolerance-only endpoint joining.
    distance = (tri-center)@normal
    hits = tri[(distance.min(axis=1)<0)&(distance.max(axis=1)>0)]
    segments = []
    for t in hits:
        d=(t-center)@normal; endpoints=[]
        for i,j in [(0,1),(1,2),(2,0)]:
            if d[i]*d[j]<0:
                endpoints.append(t[i]+d[i]/(d[i]-d[j])*(t[j]-t[i]))
        assert len(endpoints)==2, 'Plane coincident with STL vertex: stop measurement preparation'
        segments.append(endpoints)
    segments=np.array(segments)
    keys=np.rint(segments.reshape(-1,3)/1e-13).astype(np.int64)
    _, first, inverse=np.unique(keys,axis=0,return_index=True,return_inverse=True)
    pts=segments.reshape(-1,3)[first]; edges=inverse.reshape(-1,2)
    adj={i:[] for i in range(len(pts))}
    for a,b in edges:
        adj[a].append(b); adj[b].append(a)
    assert all(len(v)==2 for v in adj.values()), 'Open/non-manifold flux contour'
    unseen=set(adj); loops=[]
    while unseen:
        start=min(unseen); prev=None; cur=start; loop=[]
        while cur in unseen:
            unseen.remove(cur); loop.append(cur)
            nxt=next(x for x in adj[cur] if x != prev)
            prev,cur=cur,nxt
        assert cur==start
        loops.append(pts[loop])
    return loops

def quadrature(tri, dx):
    pending=list(tri); refined=[]
    while pending:
        t=pending.pop()
        lengths=[np.linalg.norm(t[a]-t[b]) for a,b in [(0,1),(1,2),(2,0)]]
        if max(lengths)<=dx:
            refined.append(t); continue
        a,b=[(0,1),(1,2),(2,0)][int(np.argmax(lengths))]; c=3-a-b; mid=(t[a]+t[b])/2
        pending.extend([np.array([t[a],mid,t[c]]),np.array([mid,t[b],t[c]])])
    t=np.array(refined); area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
    bary=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])
    return np.einsum('av,nvc->nac',bary,t).reshape(-1,3),np.repeat(area/3,3),t

def main(run):
    physical=json.loads((run/'contracts/physical_bc_contract.json').read_text())
    geom=json.loads((run/'contracts/geometry_reuse_contract.json').read_text())
    assert physical['status']==geom['status']=='PASS'
    source_contexts={
      'examples/pipeflow/pipeflow.cpp':'periodic/body-force pipe; unsuitable physical forcing here',
      'examples/pipeflow_with_preinlet/pipeflow_with_preinlet.cpp':'preinlet and Cartesian 0N pressure',
      'examples/curvedflow_with_preinlet/curvedflow_with_preinlet.cpp':'hardcoded Cartesian pressure box; not arbitrary cap',
      'cases/stl_preinlet/stl_preinlet.cpp':'STL with Cartesian preinlet pressure',
      'cases/AR2/AR2.cpp':'Cartesian 0N pressure reference',
      'cases/preinlet_shear/preinlet_shear.cpp':'Cartesian 0P pressure / planar velocity',
      'palabos/src/boundaryCondition/boundaryInstantiator3D.h':'1037-1040 planar domain precondition; 1272-1312 directions 0/1/2 orientations +/-',
      'palabos/src/boundaryCondition/zouHeBoundary3D.hh':'WrappedZouHeBoundaryManager pressure is direction/orientation-specialized',
      'palabos/examples/showCases/aneurysm/aneurysm.cpp':'95-137,285-365: native Guo + plug/Poiseuille + constant-density Neumann-velocity multi-opening flow',
      'palabos/src/offLattice/triangleBoundary3D.h':'public tagDomain, getTag, BoundaryProfiles3D.defineProfile',
      'palabos/src/offLattice/triangleBoundary3D.hh':'724+: tagDomain only triangles whose three vertices pass; TriangleFlowShape exact triangle intersections and profiles',
      'palabos/src/offLattice/offLatticeBoundaryProfiles3D.hh':'VelocityPlugProfile; 527-556 DensityNeumann fixed density, NOT zero-gradient pressure',
      'palabos/src/offLattice/guoOffLatticeModel3D.hh':'114-190 dry-node stencil; 485+ density override; 513+ normal momentum; existing population reconstruction',
      'palabos/src/offLattice/offLatticeBoundaryCondition3D.hh':'insert() no packed-field arguments supported; computePressure subtracts MEAN density, unsuitable for our fixed gauge reference',
      'mechanics/constantConversion.cpp':'36-54: nuLU=(tau-.5)/3, dt=nuLU*dx²/nuPhys, dm=rhoPhys*dx³, df=dm*dx/dt²',
      'palabos/src/latticeBoltzmann/nearestNeighborLattices3D.hh':'D3Q19 weights, velocities, cs²=1/3',
      'palabos/src/latticeBoltzmann/dynamicsTemplates.h':'BGK second-order equilibrium; second moment = cs²*rho*I + rho*u*u',
      'palabos/src/basicDynamics/isoThermalDynamics.hh':'360-365 BGKdynamics::computeEquilibrium delegates to bgk_ma2_equilibrium',
    }
    evidence=[]
    for rel,context in source_contexts.items():
        p=HC/rel; text=p.read_text()
        evidence.append({'path':str(p),'sha256':sha(p),'context':context,'lines':len(text.splitlines())})
    ports=geom['ports']
    feasibility={'status':'PASS','kind':'existing public API feasibility; runtime mapping/stencil gate still mandatory',
      'method':'native GuoOffLatticeModel3D regularized second order, BGK D3Q19',
      'custom_population_reconstruction':False,'evidence':evidence,'ports':{}}
    md=['# BC API feasibility — PASS (public API; Stage A runtime verification required)',
      'HemoCell current bundled Palabos provides arbitrary-triangle Guo off-lattice boundary reconstruction. The official aneurysm example uses velocity inlets and constant-density pressure outlets. No new population reconstruction is introduced.',
      'Cartesian addPressureBoundary0P/0N/1P/1N/2P/2N specialize missing populations by coordinate direction and require a planar Cartesian Box3D. Small boxes do not remove that assumption. They are not used here.',
      'Only frozen Step2 cap triangle IDs receive port tags; all other surface triangles remain stationary no-slip. No hole discovery, synthetic circles, new cap, or box-face port is used.',
      'Step2 closed lumen defines physical volume. Step2 opened mask and all 497 port labels are preserved as verified port diagnostics. Native Guo completes outerBorder ghost nodes (official aneurysm arrangement); they are numerical support, not added physical lumen. Runtime checks must prove unchanged closed flags and cap IDs. A width-2 envelope changes storage only.']
    for name,p in ports.items():
        method='Q-normalized native VelocityPlugProfile3D' if name=='inlet' else 'native DensityNeumannBoundaryProfile3D: fixed rho, normal momentum extrapolation'
        feasibility['ports'][name]={'target':physical['inlet']['condition_type'] if name=='inlet' else 'pressure','normal':p['normal'],'triangle_ids':p['triangle_ids'],'supported':True,'method':method,'limitations':'cap edge stencils mix neighboring wall and cap data; finite-grid error must be measured; role ASSUMED'}
        md.append(f"## {name}\nNormal: {p['normal']}; frozen label {p['label']}, {len(p['triangle_ids'])} cap triangles. Target: {feasibility['ports'][name]['target']}. Public API: **YES**, {method}. Uses actual surface normals, without dominant-axis substitution. Finite-grid cap/wall junction errors and short transients remain to be checked.")
    md.append('## Local source evidence\n'+'\n'.join(f"- `{e['path']}`: {e['context']} (SHA256 {e['sha256']})." for e in evidence))
    write_json(run/'contracts/bc_api_feasibility.json',feasibility)
    (run/'contracts/bc_api_feasibility.md').write_text('\n\n'.join(md)+'\n')

    dx=geom['dx_effective_m']; rho=physical['fluid']['density_kg_m3']; nu=physical['fluid']['kinematic_viscosity_m2_s']
    cs2=1/3; tau=1.; nulu=cs2*(tau-.5); dt=nulu*dx*dx/nu
    velunit=dx/dt; punit=rho*velunit**2; qunit=dx**3/dt; q=physical['inlet']['target_volume_flow_m3_s']
    tri=read_stl(geom['input_stl']); inlettri=tri[ports['inlet']['triangle_ids']]
    areas=np.linalg.norm(np.cross(inlettri[:,1]-inlettri[:,0],inlettri[:,2]-inlettri[:,0]),axis=1)/2
    area=float(areas.sum()); ume=q/area; diameter=2*np.sqrt(area/np.pi)
    unit={'status':'PASS','dx_status':'CANDIDATE_BASELINE','tau_status':'STEP3_CANDIDATE_TAU','dx_requested_m':geom['dx_requested_m'],'dx_effective_m':dx,
      'rho_phys_kg_m3':rho,'nu_phys_m2_s':nu,'cs2':cs2,'tau':tau,'nu_lu':nulu,'dt_s':dt,
      'reference_physical_gauge_pressure_pa':0.,'reference_lattice_density':1.,'gauge_reference':physical['pressure_reference']['zero_definition'],
      'velocity_conversion':{'formula':'u_phys=u_lu*dx/dt','m_s_per_lu':velunit},
      'pressure_conversion':{'formula':'p_gauge_phys=(rho_lu-1)*cs2*rho_phys*(dx/dt)^2; p_gauge_lu=p_phys/[rho_phys*(dx/dt)^2]','pa_per_lu':punit,'note':'Fixed reference rho=1; do NOT use offLatticeBC.computePressure() which subtracts instantaneous average density.'},
      'flow_conversion':{'formula':'Q_phys=Q_lu*dx^3/dt','m3_s_per_lu':qunit},
      'mass_conversion':{'kg_per_lu':rho*dx**3},'force_conversion':{'newton_per_lu':rho*dx**4/dt**2,'used':False},
      'target_Q_phys':q,'target_Q_lu':q/qunit,'outlets':{},'sources':[e for e in evidence if any(s in e['path'] for s in ['constantConversion','nearestNeighbor','dynamicsTemplates','isoThermalDynamics','offLatticeBoundaryCondition3D.hh'])]}
    for name,p in physical['outlets'].items():
        pg=p['pressure_pa']; unit['outlets'][name]={'pressure_phys_pa':pg,'pressure_gauge_lu':pg/punit,'density_lu':1+pg/punit/cs2}
    densities=[1.]+[v['density_lu'] for v in unit['outlets'].values()]
    # Pressure acceleration bound over 1000 steps assuming the entire pressure range
    # acts across a single dx, without viscous damping. Conservative startup estimate.
    pressure_range=max(0.,*[v['pressure_pa'] for v in physical['outlets'].values()])-min(0.,*[v['pressure_pa'] for v in physical['outlets'].values()])
    pbound=ume+pressure_range/(rho*dx)*1000*dt
    pre={'status':'PASS','expected_inlet_mean_velocity_m_s':ume,'inlet_area_m2':area,'equivalent_diameter_m':diameter,'Re_inlet':ume*diameter/nu,
      'expected_inlet_u_lu':ume/velunit,'expected_inlet_mach':ume/velunit/np.sqrt(cs2),
      'conservative_1000_step_speed_bound_m_s':pbound,'max_expected_mach':pbound/velunit/np.sqrt(cs2),
      'pressure_speed_bound_assumption':'Full gauge pressure range across one dx for 1000dt without viscous damping; startup screening bound, not expected solution or mathematical stability proof.',
      'rho_lu_min_expected':min(densities),'rho_lu_max_expected':max(densities),'max_pressure_density_deviation':max(abs(np.array(densities)-1)),
      'tau':tau,'Q_lu':q/qunit,'runtime_mach_limit':.05,'runtime_density_limits':[.99,1.01],
      'density_limit_rationale':'1% compressibility ceiling, much larger than imposed pressure offsets; selected for this candidate before execution, not inherited from old CFD.'}
    assert tau>.5 and pre['max_expected_mach']<.05 and min(densities)>.99 and max(densities)<1.01
    write_json(run/'contracts/lattice_unit_contract.json',unit); write_json(run/'diagnostics/lbm_preflight.json',pre)
    origin=np.array(geom['physical_origin_m']); nx,ny,nz=geom['lattice_shape']
    native=np.memmap(S2/'diagnostics/palabos_native_flags.u8',dtype='u1',mode='r',shape=(nz,ny,nx))
    metadata=json.loads((S2/'diagnostics/cap_identification.json').read_text())['ports']
    allquad=[]; planes=[]
    for name,p in ports.items():
        normal=np.array(p['normal']); center=np.array(p['center']); plane_center=center-4*dx*normal
        loops=cut_loops(tri,plane_center,normal); candidates=[]
        for loop in loops:
            t=triangulate_loop(loop)
            if projected_inside(plane_center[None,:],t,normal)[0]:
                candidates.append((loop,t))
        assert len(candidates)==1, f'{name}: ambiguous intersected lumen contour'
        loop,t=candidates[0]
        cap_radius=np.linalg.norm(tri[p['triangle_ids']].reshape(-1,3)-center,axis=1).max()
        assert np.linalg.norm(loop-plane_center,axis=1).max()<2*cap_radius, 'Plane reaches neighboring branch'
        ext=metadata[name]['extension_length_actual_s4_um']*1e-6
        assert 4*dx < ext/4
        points,weights,refined=quadrature(t,dx); lu=(points-origin)/dx
        bases=np.floor(lu).astype(int); flags=[]
        for ox in [0,1]:
            for oy in [0,1]:
                for oz in [0,1]:
                    ix=bases+np.array([ox,oy,oz]); assert np.all(ix>=0) and np.all(ix<np.array([nx,ny,nz]))
                    flags.extend(native[ix[:,2],ix[:,1],ix[:,0]].tolist())
        assert set(flags)<={1,2,3,4}, 'Unallocated/undetermined interpolation node'
        pdata={'port':name,'label':p['label'],'outward_normal':p['normal'],'cap_center_m':p['center'],'center_m':plane_center.tolist(),
          'offset_inward_lu':4,'offset_inward_m':4*dx,'known_extension_length_m':ext,'full_plane_contour_count':len(loops),
          'selected_contour_count':1,'selection':'Unique actual STL contour containing projected frozen cap center; extent <2 cap radius; distance inside first quarter of known artificial extension',
          'area_m2':float(weights.sum()),'area_lu2':float(weights.sum()/dx**2),'quadrature_points':len(points),'refined_triangles':len(refined),
          'native_interpolation_node_flag_counts':{str(c):flags.count(c) for c in sorted(set(flags))},'max_quadrature_triangle_edge_m':dx}
        planes.append(pdata)
        for pos,w in zip(lu,weights/dx**2): allquad.append([p['label'],*pos,w,*normal])
        vtkpoints=vtk.vtkPoints(); vtkpoints.SetDataTypeToDouble(); cells=vtk.vtkCellArray()
        for face in t:
            cells.InsertNextCell(3)
            for v in face: cells.InsertCellPoint(vtkpoints.InsertNextPoint(*v))
        poly=vtk.vtkPolyData(); poly.SetPoints(vtkpoints); poly.SetPolys(cells)
        writer=vtk.vtkXMLPolyDataWriter(); writer.SetFileName(str(run/f'output/measurement_{name}.vtp')); writer.SetInputData(poly); writer.SetDataModeToBinary(); assert writer.Write()==1
    np.savetxt(run/'contracts/flux_quadrature.tsv',allquad,fmt=['%d']+['%.17g']*7)
    write_json(run/'contracts/flux_measurement_contract.json',{'status':'PASS','method':'STL cross-section triangle quadrature + trilinear live lattice velocity','planes':planes,'dx_m':dx,'dt_s':dt,'quadrature_sha256':sha(run/'contracts/flux_quadrature.tsv')})
    (run/'contracts/flux_measurement_contract.md').write_text('# Flux measurement contract\n\nIndependent live-field integration, not mean velocity times contract area. Each plane is 4 effective dx inward from the unchanged Step2 cap center along its outward normal. Full frozen STL is intersected with that plane; the unique local closed contour containing the projected center is triangulated. Neighboring contours are excluded and recorded, without changing vessel geometry. Every integration triangle has edges ≤dx and uses three degree-2 quadrature points. Each point samples 8 live lattice velocities by trilinear interpolation. Outer-border ghost velocities participate where needed for interpolation; outside NoDynamics velocity is zero. Ghosts are excluded from physical mass and ParaView fluid mask.\n\nQ_outward_LU=sum(area_LU² * u_LU dot n_out); Q_phys=Q_LU*dx³/dt. Qin=-Q_outward(inlet); outlets preserve signed outward values. Closure=abs(Qin-sum(Qout))/max(abs(Qin),Qtarget*1e-12). Target error=abs(Qin-Qtarget)/Qtarget. Effective measured areas, normals, quadrature counts, native sampling flags, and extension checks are in flux_measurement_contract.json; actual sampling coordinates/weights are in flux_quadrature.tsv.\n\nThese planes cut off short end sections. Their instantaneous flux imbalance is not an exact discrete whole-lumen mass residual; physical lumen mass and its relative drift are reported separately. At startup, pressure waves, storage and cap-to-plane delay can produce substantial closure and target error. There is no steady-state acceptance claim. Finite-grid quadrature and wall interpolation errors remain unverified by refinement.\n')
    # Machine input: no coordinates or physical constants hardcoded in executable.
    with (run/'contracts/solver_parameters.txt').open('w') as f:
        f.write(f"{geom['input_stl']}\n{S2}\n")
        f.write(' '.join(map(str,[nx,ny,nz,*origin,dx,dt,tau,rho,q,.05,.99,1.01,10,1000]))+'\n')
        for name,p in ports.items():
            dens=1 if name=='inlet' else unit['outlets'][name]['density_lu']
            f.write(' '.join(map(str,[p['label'],*p['normal'],dens,len(p['triangle_ids'])]))+'\n')
            for tid in p['triangle_ids']:
                f.write(' '.join(map(str,[tid,*tri[tid].ravel()]))+'\n')
    execution={'status':'FROZEN_BEFORE_BUILD_AND_TIMESTEPS','mpi_ranks':1,'mpi2':'UNVERIFIED','executable':str(run/'build/vascular_pure_fluid'),
      'system_PATH':'/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin','mpirun':'/usr/bin/mpirun','mpicxx':'/usr/bin/mpicxx',
      'source_head':'fd21a850a16d0ba17ef3d521123badd53864a3a2a2','hemocell_head':'5a410848bd5c57d5ae1c171112e78eab4a82e650',
      'dx':dx,'dt':dt,'tau':tau,'rho':rho,'nu':nu,'Qtarget':q,'outlet_pressures_pa':{k:v['pressure_pa'] for k,v in physical['outlets'].items()},
      'stages':[0,1,10,100,1000],'snapshot_steps':[0,10,100,1000],'diagnostic_interval':10,'safety_check_interval':1,
      'initialization':'uniform rest, rhoLU=1','NUMERICAL_STARTUP_RAMP':'linear from 0 at stage A to full physical target at step 10; both Q and gauge pressures; no tuning',
      'inlet_profile_status':'STEP3_NUMERICAL_INLET_PROFILE_ASSUMPTION','inlet_realization':'uniform normal profile normalized by actual tagged native cap projected area after official 0.001 LU inflation, integral equals Qtarget; frozen physical cap area remains unchanged',
      'inlet_profile_rationale':'Current target Q is specified without experimentally established profile. Constant native plug profile is the minimal auditable numerical assumption on the actual noncircular cap; artificial extension permits downstream profile development.',
      'wall':'native Guo off-lattice stationary no-slip','outlet':'native densityNeumann fixed individual gauge-density values, normal momentum extrapolation',
      'short_run_acceptance':{'all_stage_safety_checks':True,'all_four_tag_profiles_actually_called':True,'geometric_profile_integral_relative_error_max':1e-12,
        'final_inlet_flow_positive':True,'final_plane_inlet_target_error_max':.25,'note':'25% is a provisional short-run consistency tolerance including startup and finite-grid error; not a validated accuracy level. Failure remains failure, no retuning.',
        'relative_mass_drift_abs_max':.01,'outlet_backflow':'record and interpret pressure/transient; not automatic failure','closure':'report; no steady threshold in short startup'},
      'steady_state_acceptance':'NOT_ASSESSED','max_timesteps':1000,'cell_count':0,'automatic_parameter_tuning':False,'native_stencil_envelope':2,
      'geometry_reuse':'same Step2 closed pipeline reconstructed solely to build native TriangleHash/flow shape; exact bulk flags, transform, triangle coordinates and labels checked before any timestep'}
    # Read actual audited head instead of trusting textual constants in this script.
    import subprocess
    execution['source_head']=subprocess.check_output(['/usr/bin/git','-C','/home/lzy/projects/ulm_3D_vascular','rev-parse','HEAD'],text=True).strip()
    write_json(run/'contracts/execution_contract.json',execution)
    geom['method']=execution['geometry_reuse']; write_json(run/'contracts/geometry_reuse_contract.json',geom)
    print(json.dumps({'bc_api':'PASS','lattice_units':'PASS','preflight':pre,'flux_planes':planes},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('run',type=Path)
    main(parser.parse_args().run)
