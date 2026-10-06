"""Machine-readable audit of verified deployed source; no source mutations."""
from common import *
import xml.etree.ElementTree as ET

def main():
    root=SOURCE/'Code/Source/solver';manifest=read(HERE/'EXACT_DEPLOYED_SOURCE_MANIFEST.json');inputs=read(HERE/'INPUT_SHA256.json')
    specs={
      'assembly_dispatch':('eq_assem.cpp','fluid::construct_fluid'),
      'fluid_entry':('fluid.cpp','void construct_fluid'),
      'vms_switch':('fluid.cpp','if (lM.nFs == 1)'),
      'momentum':('fluid.cpp','void fluid_3d_m'),
      'continuity':('fluid.cpp','void fluid_3d_c'),
      'tau_momentum':('fluid.cpp','double kT = 4.0 * pow(ctM/dt,2.0)'),
      'tau_continuity':('fluid.cpp','tauM = 1.0 / (rho * sqrt( kT + kU + kS ))'),
      'time_parameters':('initialize.cpp','eq.am ='),
      'metric':('nn.cpp','void gnn('),
      'metric_order':('nn.cpp','tet4{1, 2, 3, 0}'),
      'mesh_orientation':('read_msh.cpp','void check_tet_conn'),
      'TH_reader':('read_files.cpp','use_taylor_hood_type_basis.defined()'),
      'TH_pressure':('fs.cpp','case ElementType::TET10:'),
      'TH_edge_pressure':('fs.cpp','void thood_val_rc'),
      'TH_integrator':('Integrator.cpp','fs::thood_val_rc'),
      'TH_pressure_corrector':('Integrator.cpp','void Integrator::corrector_taylor_hood'),
      'PETSc_assembly':('PetscLinearAlgebra.cpp','void PetscLinearAlgebra::assemble'),
      'VTK_TET10_write':('VtkData.cpp','return VTK_QUADRATIC_TETRA'),
      'VTK_TET10_read':('vtk_xml_parser.cpp','{VTK_QUADRATIC_TETRA, 10}'),
      'restart':('output.cpp','restart_file.write((char*)Yn.data()'),
      'arbitrary_basis':('FE/Basis/BasisFactory.cpp','std::make_unique<LagrangeBasis>(req.topology, order)'),
      'TET4_quadrature':('nn_elem_gip.h','{ElementType::TET4,'),
      'TET4_quadrature_default':('Parameters.cpp','"Quadrature_modifier_TET4"'),
    }
    refs={}
    for key,(rel,needle) in specs.items():
        p=root/rel;h=sha(p);assert h==manifest['Code/Source/solver/'+rel],rel;inputs[str(p)]=h
        lines=p.read_text().splitlines();hits=[i+1 for i,s in enumerate(lines) if needle in s];assert hits,(key,needle)
        refs[key]=dict(file=rel,sha256=h,lines=hits,needle=needle)
    xml=BASE/'cases/balanced_pressure_final/run/solver.xml';e=ET.parse(xml).getroot();g=e.find('GeneralSimulationParameters');eq=e.find('Add_equation')
    fluid=(root/'fluid.cpp').read_text();m=fluid[fluid.index('void fluid_3d_m'):];c=fluid[fluid.index('void fluid_3d_c'):fluid.index('void fluid_3d_m')]
    for s in [m,c]:assert 'ctM  = 1.0' in s and 'ctC  = 36.0' in s
    assert eq.attrib['type']=='fluid' and eq.find('Use_taylor_hood_type_basis') is None
    assert float(g.findtext('Time_step_size'))==.01 and float(eq.findtext('Density'))==1056 and float(eq.findtext('Viscosity/Value'))==.00345312
    result=dict(status='SOURCE_CONFIRMED',production_identity='PRODUCTION_SOURCE_IDENTITY.json',exact_source_manifest='EXACT_DEPLOYED_SOURCE_MANIFEST.json',source_commit=read(HERE/'PRODUCTION_SOURCE_IDENTITY.json')['source_commit'],references=refs,
      assembly_chain=['Integrator::assemble_equations','eq_assem::global_eq_assem Equation_fluid','fluid::construct_fluid: lM.nFs=1 -> vmsStab=true','fs::get_thood_fs iOpt=1, nn::gnn, fluid_3d_m','fs::get_thood_fs iOpt=2, nn::gnn, fluid_3d_c','PetscLinearAlgebra::assemble -> FsilsLinearAlgebra -> com_mod.Val -> PETSc CUDA solve'],
      active_parameters=dict(dt_s=.01,ctM=1.,ctC=36.,rho=1056.,mu=.00345312,nu=3.27e-6,rho_inf=.5,alpha_m=5/6,alpha_f=2/3,gamma=2/3),
      active_flags=dict(nFs=1,vmsFlag=True,eNoNw=4,eNoNq=4,Taylor_Hood=False,ALE=False,constant_viscosity=True,Brinkman_inverse_permeability=0.,URIS=False,body_force=False,backflow_coefficient=0.,global_omega_correction=False),
      solver_xml_sha256=sha(xml),boundary={b.attrib['name']:{x.tag:x.text for x in b} for b in eq.findall('Add_BC')},
      formulas=dict(J='[x0-x3,x1-x3,x2-x3] after reader orientation',Kxi='J^-T J^-1',kT='4*(ctM/dt)^2 + (Brinkman_inverse_permeability*mu/rho)^2 + urisFactorTotal^2; active=40000 s^-2',kU='u^T Kxi u',kS='ctC*(mu/rho)^2*(Kxi:Kxi)',tau_M='1/(rho*sqrt(kT+kU+kS))',tau_C='1/(tau_M*tr(Kxi))',up='-tau_M*[rho*(a_alpha_m - body_force + (u_alpha_f.grad)u_alpha_f) + grad(p) - div(mu*(G+G.T)) + inactive_Brinkman_URIS_terms]',
      fine_scale_extra='tau_B=rho/sqrt(up^T Kxi up), code epsilon floor when is_zero; ua=u+up; pa=p-tau_C*div(u)',
      continuity='R_p[a] += w*(Nq[a]*div(u) - up.grad(Nq[a])); PSPG pressure-pressure tangent w*alpha_f*gamma*dt*tau_M*gradNa.gradNb',
      grad_div='pa=p-tau_C div(u): momentum includes +tau_C div(u) div(test); LSIC/grad-div-like, not exact element mass constraint',
      generalized_alpha='u_alpha_f=(1-alpha_f)*u_old+alpha_f*u_new; a_alpha_m=(1-alpha_m)*a_old+alpha_m*a_new; tangent T1=alpha_f*gamma*dt, amd=alpha_m/T1',
      steady='a->0 does not remove explicit 4/dt^2 in tau; actual exported velocity used for OFFLINE evaluation, not asserted identical to last intermediate state'),
      JIA_ESMAILY_TRADITIONAL_FORM_MATCH=True,match_scope='Time-scale formula family: paper tauSUPG=rho*tauM and paper nu_C=tauC/rho; actual ctC=36 vs paper CI=3; production has additional VMS terms. Not identical full algorithms.',
      paper_counterpart=dict(tau='1/sqrt(omega^2+u.K.u+CI*nu^2*K:K)',omega='global ||partial_t u||L2 / ||u||L2; Eulerian, not convective material acceleration',initialization='first step omega=2/dt; previous accepted state for subsequent omega, once per step held during nonlinear iterations',ALE='paper Eq10 uses time derivative at fixed mesh coordinates and relative velocity in denominator',code_access='Paper points to MUPFES author site; current direct access returned Google sign-in, public implementation not obtained/verified this session',offline='omega=0 steady limiting case, keep actual ctC=36; NOT a time-consistent CFD solution'))
    put('STABILIZATION_SOURCE_AUDIT.json',result);put('INPUT_SHA256.json',inputs)
    put('DECISION_GATE.json',dict(status='STOP_C',source_identity=True,traditional_family=True,paper_formula_known=True,public_author_implementation_verified=False,material_local_dt_influence=False,already_time_consistent=False,reason='Target time shares <2.227% over all retained samples; removing kT changes tauM <=1.1322%, at 1189 A only0.29857%; no evidence for local dt-dominance mechanism.',new_solver_patch=False,new_pipe_CFD=0,new_BraVa_CFD=0,pipe_result='NOT_RUN',BraVa_result='NOT_RUN',nonlocal_dt_effect='NOT_EVALUABLE from local term decomposition alone'))
    print('SOURCE_AUDIT_COMPLETE',len(refs),'STOP_C')
if __name__=='__main__':main()
