"""Small numerical evidence set; eight FEM locations, five 2D reductions."""
from pathlib import Path
import ast,types,sys,csv,json
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.particle9a_provenance import require_current_flow
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block,legacy_planar_coefficients_si
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle9a_motion import augment_planar_system
from particle_3d.resistance_solver import solve_resistance
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_gap import wall_gap
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle3_cases import plane_triangle
from particle_3d.particle81_simulation import environment,dump
REPORT=ROOT/'particle_3d/reports/particle9a_2mmps'


def csv_file(name,rows):
    with (REPORT/'data'/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)


def main():
    require_current_flow()
    ref=REPORT/'reference/legacy_particle_mobility.py';tree=ast.parse(ref.read_text())
    tree.body=[node for node in tree.body if not isinstance(node,ast.Try)]
    for node in ast.walk(tree):
        if isinstance(node,ast.FunctionDef):node.decorator_list=[]
    legacy=types.ModuleType('legacy');exec(compile(tree,str(ref),'exec'),legacy.__dict__)
    reciprocity=[dict(xi=xi,**legacy_planar_coefficients_si(xi,freeze_asymptote=False)) for xi in [.001,.003,.01,.03,.1,.3,1.]]
    csv_file('reciprocity_audit.csv',reciprocity)
    n=np.array([.6,0,.8]);t=np.array([-.8,0,.6]);a=1e-6;mu=.00345312;gamma=500.;G=gamma*np.outer(t,n)
    reductions=[]
    for xi in [.003,.01,.1,.5,1.]:
        v=gamma*a*(1+xi)*t;u=np.r_[v,[0,-gamma/2,0]]
        _,_,d=planar_wall_affine_block(a,mu,xi*a,n,G,u)
        old=legacy.background_hydrodynamic_velocity_scalar(v[0]*1e6,v[2]*1e6,G[0,0],G[0,2],G[2,0],G[2,2],n[0],n[2],1.,xi,mu)
        vt=np.array(d['target_tangential_velocity_xyz'])@t;omega=d['target_tangential_omega_xyz'][1]
        oldvt=(old[0]*t[0]+old[1]*t[2])*1e-6
        c=legacy_planar_coefficients_si(xi);fs,ty,_=legacy.wall_shear_load_scalar(mu,1.,xi,gamma);base=1/(6*np.pi*mu*a)
        ev=abs(vt-oldvt);eo=abs(omega-old[2]);budget_v=old[3]*base*c['projection_entry_error']*abs(ty*1e-18/a);budget_o=old[3]*base*c['projection_entry_error']*abs(fs*1e-12)/a
        reductions.append(dict(xi=xi,Vt_2D_m_s=oldvt,Vt_P9A_m_s=vt,omega_2D_s_inv=old[2],omega_P9A_s_inv=omega,
            abs_Vt_error_m_s=ev,abs_omega_error_s_inv=eo,analytic_projection_budget_Vt_m_s=budget_v,
            analytic_projection_budget_omega_s_inv=budget_o,pass_budget=bool(ev<=budget_v+1e-15 and eo<=budget_o+1e-9)))
    csv_file('planar_2d_reduction.csv',reductions)
    wall=WallGeometry([plane_triangle()]);comparison=[]
    # 81 display samples produce a smooth figure, not trajectory experiments.
    for xi in np.geomspace(.002,2.,81):
        shape=Sphere([0,0,a*(1+xi)],a);G=np.array([[0,0,500.],[0,0,0],[0,0,0]])
        free=np.array([500*a*(1+xi),0,-1e-4,0,250.,0])
        old=assemble_v1({1:shape},{1:free},mu,wall)
        new=augment_planar_system(old,{1:shape},mu,wall,lambda x:G,wall_gap)
        p65=solve_resistance(old).velocity;p9=solve_resistance(new).velocity;d=new.planar_diagnostics[0]
        comparison.append(dict(xi=xi,p65_Vt_mm_s=p65[0]*1e3,p9a_Vt_mm_s=p9[0]*1e3,
            p65_omega_y_s_inv=p65[4],p9a_omega_y_s_inv=p9[4],normal_difference_m_s=p9[2]-p65[2],
            m_tt=d['effective_m_tt'],m_rr=d['effective_m_rr'],m_cross=d['effective_m_cross']))
    csv_file('planar_comparison.csv',comparison)
    env=environment();field=env.field;centers=field.points_m[field.tetra].mean(axis=1)
    # Fixed spatial selections and one high-gradient cell. No field alteration.
    ids=np.linspace(0,len(centers)-1,100,dtype=int)
    gradient_norm=np.linalg.norm(field.gradients_s_inv,axis=(1,2))
    branch_cell=int(np.argmax(np.linalg.norm(field.velocity_nodes_m_s[field.tetra].mean(axis=1),axis=1)))
    selected=[int(ids[k]) for k in [8,22,38,55,72,88]]+[int(np.argmax(gradient_norm)),branch_cell]
    ratios=[2.,.5,.1,.03,.03,.003,.02,.1];rows=[]
    for cell,xi in zip(selected,ratios):
        center=centers[cell];_,distance=env.wall.nearest_center_triangle(center)
        radius=distance/(1+xi);shape=Sphere(center,radius);s=field.sample(center);gap=wall_gap(shape,env.wall)
        free=np.r_[s.velocity_m_s,.5*s.vorticity_s_inv]
        old=assemble_v1({1:shape},{1:free},env.mu,env.wall)
        new=augment_planar_system(old,{1:shape},env.mu,env.wall,lambda x:s.velocity_gradient_s_inv,wall_gap)
        v65=solve_resistance(old).velocity;v9=solve_resistance(new).velocity;d=new.planar_diagnostics[0]
        speed_ratio=np.linalg.norm(v9[:3])/max(np.linalg.norm(free[:3]),1e-12)
        rows.append(dict(cell_id=cell,selection='HIGH_GRADIENT_CELL' if cell==selected[-2] else 'HIGHEST_SPEED_BRANCH_VISUAL_FOCUS' if cell==branch_cell else 'FIXED_SPATIAL_CELL',
            center_m=center.tolist(),radius_m=radius,
            u_m_s=s.velocity_m_s.tolist(),gradient_s_inv=s.velocity_gradient_s_inv.tolist(),
            p65_velocity_m_s=v65[:3].tolist(),p9a_velocity_m_s=v9[:3].tolist(),p9a_omega_s_inv=v9[3:].tolist(),
            normal_velocity_difference_m_s=float((v9[:3]-v65[:3])@gap.normal_inward),
            speed_ratio_to_bulk=float(speed_ratio),finite=bool(np.isfinite(v9).all()),**d))
    dump(REPORT/'data/real_fem_smoke.json',dict(role='EIGHT_STATIC_LOCAL_SOLVES_NOT_A_PARAMETER_STUDY',
        radius_role='DIAGNOSTIC_RADIUS_MATCHED_TO_LOCAL_CLEARANCE_NOT_A_SONOVUE_POPULATION',rows=rows))
    csv_file('real_fem_smoke.csv',[{k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows])
    summary=dict(reciprocity_max_numerator=max(abs(r['reciprocity_numerator']) for r in reciprocity),
        reciprocity_rounding_bound=reciprocity[0]['rounding_numerator_bound'],
        normal_max_difference_m_s=max(abs(r['normal_difference_m_s']) for r in comparison),
        reduction_max_Vt_error_m_s=max(r['abs_Vt_error_m_s'] for r in reductions),
        reduction_max_omega_error_s_inv=max(r['abs_omega_error_s_inv'] for r in reductions),
        reduction_within_derived_rounding_budget=all(r['pass_budget'] for r in reductions),
        fem_max_abs_normal_shear_force_N=max(abs(r['wall_normal_shear_force_N']) for r in rows),
        fem_max_speed_ratio=max(r['speed_ratio_to_bulk'] for r in rows),fem_finite=all(r['finite'] for r in rows),
        fem_normal_max_difference_m_s=max(abs(r['normal_velocity_difference_m_s']) for r in rows))
    dump(REPORT/'data/numerical_validation.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
