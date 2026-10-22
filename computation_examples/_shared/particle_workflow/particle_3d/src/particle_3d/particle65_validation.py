"""Independent scalar checks and actual LAMMPS parity for V1."""
import numpy as np
from .particle65_cases import *
from .particle65_motion import assemble_v1
from .hydrodynamic_resistance import PhysicalNearField,evaluate_block,sphere_self_diagonal
from .resistance_assembly import assemble_resistance_system
from .resistance_solver import solve_resistance
from .nearfield_handoff import handoff_constraints
from .lammps_bridge import LammpsParticleBridge
from .particle6_cases import state_errors
from .particle65_checkpoint import write_checkpoint,read_checkpoint
from .particle6_checkpoint import frozen_provenance


def static_scans():
    policy=NearFieldRegularizationV1();walls=[];pairs=[]
    chis=np.unique(np.r_[np.geomspace(.00001,.1,180),.0015,.002,.01,.02,.03,.04,.05,.06])
    for ratio in [None,1,2,4]:
        for chi in chis:
            a=1e-6;h=float(chi*a);b=a if ratio is None else ratio*a
            shapes={17:Sphere([0,0,a+h],a)} if ratio is None else {17:Sphere([0,0,0],a),203:Sphere([a+b+h,0,0],b)}
            spec=PhysicalNearField(17,None if ratio is None else 203,h,[0,0,1] if ratio is None else [-1,0,0],0.)
            z,n,r=policy.evaluate(spec,shapes,MU);old,_,oldrow=evaluate_block(spec,shapes,MU)
            lead,_,_=evaluate_block(PhysicalNearField(spec.particle_i_id,spec.particle_j_id,h,n,0.,'VALIDATION_GAPS_ONLY'),shapes,MU)
            free={i:np.r_[(np.array([0,0,-1e-4]) if ratio is None else np.array([5e-5 if i==17 else -5e-5,0,0])),np.zeros(3)] for i in shapes}
            # Solve actual V1 matrix at the exact prescribed geometries; raw geometry
            # differs from the prescribed h only by float64 subtraction.
            wall=WallGeometry([plane_triangle()]) if ratio is None else None
            system=assemble_v1(shapes,free,MU,wall,policy=policy);solution=solve_resistance(system)
            row=np.r_[n,np.zeros(3)] if ratio is None else np.r_[n,np.zeros(3),-n,np.zeros(3)]
            self_inv=sum(1/(6*np.pi*MU*s.radius_m) for s in shapes.values())
            actual_z=system.blocks[0]['coefficient_kg_s']
            actual_reference=1/(1+actual_z*self_inv)
            measured=float((row@solution.velocity)/(row@system.free))
            r.update(radius_i_m=a,radius_j_m=None if ratio is None else b,radius_ratio=ratio,
                old_p5_default_zeta_kg_s=old,old_p5_leading_uncapped_zeta_kg_s=lead,
                old_p5_default_vn_ratio=1/(1+old*self_inv),v1_unconstrained_vn_ratio=measured,
                analytic_v1_unconstrained_vn_ratio=actual_reference,scalar_solution_error=abs(measured-actual_reference),
                below_lower_role='DIAGNOSTIC_STATIC_SCAN_ONLY_NOT_ACCEPTED_DYNAMICS',
                old_p5_eligibility_role=oldrow['eligibility_role'])
            (walls if ratio is None else pairs).append(r)
    return walls,pairs


def bridge_case(folder):
    shapes={17:Sphere([-1.001e-6,0,1.002e-6],1e-6),203:Sphere([1.001e-6,0,1.002e-6],1e-6),
        901:Sphere([10e-6,0,3e-6],1e-6)}
    def provider(i,s,t):return np.array([2e-5 if i==17 else -2e-5 if i==203 else 0,2e-6,-1e-5 if i!=901 else 0]),np.array([.1,.2,.3])
    wall=WallGeometry([plane_triangle()]);query=ValidationNeighborPolicy(20e-6,.5e-6,'V1_PARITY_ONLY_EXTRA_CANDIDATES_INTENTIONAL')
    particles=records_for(shapes,provider);a=Particle65Stepper(particles,query,provider,MU,wall=wall);a.matrix_audit=[]
    provenance=frozen_provenance(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    errors=[];matrix_rows=[];snapshots=[];restart_error=None
    bridge=LammpsParticleBridge(particles[::-1],query);b=Particle65Stepper(bridge.read(),query,provider,MU,wall=wall,bridge=bridge);b.matrix_audit=[]
    allforces=[];allcommands=[]
    try:
        for k in range(1,7):
            a.step_to(k*.0005);b.step_to(k*.0005)
            err=state_errors(a.read(),b.read(),steps=k);errors.append(dict(step=k,**err))
            ar,br=a.matrix_audit[-1],b.matrix_audit[-1]
            diff={name:float(np.max(np.abs(np.asarray(ar[name])-br[name]),initial=0.)) for name in ['R','b','U','J']}
            raw=sorted(bridge.rebuild());eligible=ar['eligible_pairs'];filtered=br['eligible_pairs']
            matrix_rows.append(dict(step=k,**diff,standalone_eligible_pairs=eligible,lammps_candidates=raw,
                filtered_nearfield_pairs=filtered,constraints_equal=ar['constraint_ids']==br['constraint_ids'],
                standalone=ar,bridge=br))
            snapshots.append(dict(time_s=a.time_s,standalone=[p.to_dict() for p in a.read()],bridge=[p.to_dict() for p in b.read()]))
            if k==3:
                before=b.read();write_checkpoint(folder,b,provenance,REPO);allforces+=bridge.force_audits;allcommands+=bridge.commands
                bridge.close();b=read_checkpoint(folder,provenance,provider,MU,wall=wall);bridge=b.bridge;b.matrix_audit=[]
                restart_error=state_errors(before,b.read(),steps=k)
        allforces+=bridge.force_audits;allcommands+=bridge.commands
        return dict(status='PASS' if all(e['exact_equal'] for e in errors) and all(all(r[k]==0 for k in ['R','b','U','J']) and r['constraints_equal'] and r['standalone_eligible_pairs']==r['filtered_nearfield_pairs'] for r in matrix_rows) else 'FAIL',
            errors=errors,matrix_rows=matrix_rows,snapshots=snapshots,policy=query.to_dict(),restart_error=restart_error,
            restart_at_step=3,total_steps=6,force_audits=allforces,commands=allcommands,provenance=provenance)
    finally:bridge.close()
