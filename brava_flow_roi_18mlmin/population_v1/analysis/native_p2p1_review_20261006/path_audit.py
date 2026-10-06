"""Bind the inspected native call chain to exact retained-source hashes and fixture."""
from common import *

def main():
    unit=read(HERE/'NATIVE_UNIT_FIXTURE.json')
    assert unit['status']=='PASS'
    selections={
      'configuration':('read_files.cpp',1738,1778),
      'TET10_class':('nn_elem_props.h',58,66),
      'VTK_reader_and_faces':('vtk_xml_parser.cpp',38,119),
      'orientation':('read_msh.cpp',657,704),
      'face_node_selection':('read_msh.cpp',447,520),
      'node_distribution':('distribute.cpp',209,245),
      'native_basis_order':('nn.cpp',150,202),
      'mixed_quadrature':('fs.cpp',73,197),
      'pressure_space':('fs.cpp',320,383),
      'fluid_dispatch':('eq_assem.cpp',398,412),
      'mixed_fluid_assembly':('fluid.cpp',480,754),
      'momentum':('fluid.cpp',1768,2226),
      'continuity':('fluid.cpp',1443,1765),
      'strong_Dirichlet_all_face_nodes':('set_bc.cpp',995,1144),
      'edge_pressure_matrix':('fs.cpp',394,468),
      'edge_pressure_before_solve':('Integrator.cpp',128,156),
      'edge_pressure_output_correction':('Integrator.cpp',1087,1208),
      'PETSc_assembly':('PetscLinearAlgebra.cpp',44,85),
      'FSILS_to_sparse':('FsilsLinearAlgebra.cpp',30,65),
      'PETSc_local_to_global':('petsc_impl.cpp',402,536),
      'output_all_nodes':('vtk_xml.cpp',1498,1545),
      'VTK_TET10_serialization':('VtkData.cpp',713,746),
      'checkpoint_all_nodes':('output.cpp',250,347),
      'restart_reader':('initialize.cpp',90,150),
    }
    evidence={}
    for label,(name,a,b) in selections.items():
        p=SOURCE/name;lines=p.read_text().splitlines()
        evidence[label]=dict(file='Code/Source/solver/'+name,sha256=sha(p),line_start=a,line_end=b,
                             excerpt='\n'.join(f'{i+1}: {lines[i]}' for i in range(a-1,min(b,len(lines)))))
    p=SOURCE/'petsc_impl.cpp';lines=p.read_text().splitlines()
    needles=['MatSetFromOptions','VecSetFromOptions','MatSetValues','VecCreate','MatCreate','cuda','cusparse']
    evidence['GPU_sparse_route']=dict(file='Code/Source/solver/petsc_impl.cpp',sha256=sha(p),
        matches=[dict(line=i+1,text=l) for i,l in enumerate(lines) if any(n in l for n in needles)])
    result=dict(status='NATIVE_P2P1_PATH_CONFIRMED',scope='Native basis, mixed spaces, 3D local fluid assembly, pressure-edge matrix handling executed; full GPU CFD is a later gate and remains NOT_RUN.',
      source_identity=read(HERE/'SOURCE_IDENTITY.json'),native_unit_fixture=unit,call_chain=[
      'VTK_QUADRATIC_TETRA -> eNoN10, TRI6 boundary connectivity',
      'Use_taylor_hood_type_basis true -> mesh.nFs=2',
      'fs velocity TET10/10 and pressure TET4/4 at common reference quadrature points',
      'Integrator -> eq_assem -> fluid::construct_fluid -> fluid_3d_m/fluid_3d_c',
      'PetscLinearAlgebra::assemble -> FsilsLinearAlgebra -> Val/R; pressure-edge correction before solve',
      'Existing production PETSc options select aijcusparse/cuda; same 4-by-4 stored node blocks',
      'VTK all 10 velocity nodes; pressure edge values interpolated; full nodal Yn/An/Dn checkpoint'],
      velocity_reference_corners=[[1,0,0],[0,1,0],[0,0,1],[0,0,0]],
      edge_order=[[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]],
      solver_to_FE_basis=[1,2,3,0,5,9,8,4,6,7],
      ownership='Input supplies one global ID per unique shared edge; reader/partitioner distribute all eNoN nodes. Solver does not create midpoint DOFs from TET4.',
      wall_BC='All face.gN nodes in the inherited WALL TRI6 faces receive the three zero velocity components. An edge with wall endpoints but no wall face incidence is not automatically fixed.',
      storage='Independent global unknowns 3*N_TET10 + N_corner; stored slots 4*N_TET10 with pressure edge identity rows. Per-element active 34/stored40, not 40 independent Taylor-Hood DOFs.',
      native_stabilization_semantics='nFs=2 selects vmsFlag=false: continuity tauM/PSPG disabled; momentum tauC/tauB zero, ua=u, pa=p, while momentum residual tauM/up terms remain. No parameter patch; native basis switch changes weak-form terms as well as polynomial degree.',
      curved_geometry_caveat='construct_fluid uses xwl=10 nodes and xql=4 corner nodes with separate gnn calls; velocity is curved but pressure geometry is affine. Momentum weights use velocity Jacobian. In continuity, only g=0 pressure gnn overwrites Jacobian; pressure lShpF=true skips that call at later points, whose weights therefore use velocity Jacobian. Affine fixture does not validate curved geometry; see separate D2 native assembly probe.',
      actual_P2_GPU_run='NOT_RUN',production_patch=False,evidence=evidence,
      regression_inventory='NATIVE_REGRESSION_INVENTORY.json')
    put('NATIVE_P2P1_PATH_AUDIT.json',result)
    print(result['status'])

if __name__=='__main__':main()
