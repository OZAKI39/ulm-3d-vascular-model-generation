#!/usr/bin/python3
"""Export sparse native snapshots without allocating the full 68M-cell field.

VTU cells are the exact Step2 voxels (not a new fluid grid); global lattice IDs
permit exact mask/label verification. Only RUN_DIR is written.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy

def js(path,obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def main(run):
    g=json.loads((run/'contracts/geometry_reuse_contract.json').read_text())
    u=json.loads((run/'contracts/lattice_unit_contract.json').read_text())
    s2=Path('/home/lzy/projects/compre_output/step2/20260912_225418')
    nx,ny,nz=g['lattice_shape']; dx=g['dx_effective_m']; origin=np.array(g['physical_origin_m'])
    opened=np.memmap(s2/'diagnostics/opened_flag_matrix.u8',mode='r',dtype='u1')
    closed=np.memmap(s2/'diagnostics/closed_flag_matrix.u8',mode='r',dtype='u1')
    labels=np.memmap(s2/'diagnostics/port_label_field.u8',mode='r',dtype='u1')
    ids=np.flatnonzero(opened).astype(np.int64); fluid=np.array(closed[ids]); portlabels=np.array(labels[ids])
    ijk=np.column_stack([ids%nx,(ids//nx)%ny,ids//(nx*ny)])
    # VTK hexahedron connectivity; voxel centers are origin+ijk*dx.
    offsets=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    cornerijk=(ijk[:,None,:]+offsets).reshape(-1,3)
    cornerid=(cornerijk[:,2]*(ny+1)+cornerijk[:,1])*(nx+1)+cornerijk[:,0]
    unique,inv=np.unique(cornerid,return_inverse=True)
    corners=np.column_stack([unique%(nx+1),(unique//(nx+1))%(ny+1),unique//((nx+1)*(ny+1))])
    xyz=origin+(corners-.5)*dx
    points=vtk.vtkPoints(); points.SetData(numpy_to_vtk(xyz,deep=True))
    cells=vtk.vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(len(ids)+1,dtype=np.int64)*8,deep=True),numpy_to_vtkIdTypeArray(inv.astype(np.int64),deep=True))
    grid=vtk.vtkUnstructuredGrid();grid.SetPoints(points);grid.SetCells(vtk.VTK_HEXAHEDRON,cells)
    def add(name,array):
        a=numpy_to_vtk(np.ascontiguousarray(array),deep=True);a.SetName(name);grid.GetCellData().AddArray(a)
    add('GlobalLatticeId',ids);add('FluidMask',fluid);add('Step2OpenedEnvelope',np.ones(len(ids),dtype='u1'))
    add('PortLabel',portlabels);add('PhysicalFieldValid',fluid.copy())
    dtype=np.dtype([('index','<u8'),('rho','<f8'),('u','<f8',3)])
    checks=[];collection=[]
    for source in sorted((run/'diagnostics').glob('fields_*.bin'),key=lambda p:int(p.stem.split('_')[1])):
        step=int(source.stem.split('_')[1]); num=int(np.fromfile(source,dtype='<u8',count=1)[0])
        arr=np.fromfile(source,dtype=dtype,offset=8)
        assert len(arr)==num==int(fluid.sum()) and source.stat().st_size==8+40*num
        order=np.argsort(arr['index']); arr=arr[order]; pos=np.searchsorted(ids,arr['index'])
        assert np.array_equal(ids[pos],arr['index']) and np.array_equal(ids[fluid==1],arr['index'])
        rho=np.zeros(len(ids),dtype=np.float64);v=np.zeros((len(ids),3),dtype=np.float64)
        rho[pos]=arr['rho'];v[pos]=arr['u']*u['velocity_conversion']['m_s_per_lu']
        pressure=np.zeros(len(ids),dtype=np.float64);pressure[pos]=(arr['rho']-1)*u['cs2']*u['pressure_conversion']['pa_per_lu']
        assert np.isfinite(rho).all() and np.isfinite(v).all() and np.isfinite(pressure).all()
        for name,value in [('DensityLU',rho),('Density_kg_m3',rho*u['rho_phys_kg_m3']),('Velocity_m_s',v),('VelocityMagnitude_m_s',np.linalg.norm(v,axis=1)),('GaugePressure_Pa',pressure)]:
            add(name,value)
        grid.GetCellData().SetActiveVectors('Velocity_m_s');grid.GetCellData().SetActiveScalars('VelocityMagnitude_m_s')
        time=vtk.vtkDoubleArray();time.SetName('TimeValue');time.InsertNextValue(step*u['dt_s']);grid.GetFieldData().AddArray(time)
        dest=run/f'output/flow_step_{step:04d}.vtu'
        writer=vtk.vtkXMLUnstructuredGridWriter();writer.SetFileName(str(dest));writer.SetInputData(grid);writer.SetDataModeToBinary();writer.SetCompressorTypeToZLib()
        assert writer.Write()==1
        reader=vtk.vtkXMLUnstructuredGridReader();reader.SetFileName(str(dest));reader.Update(); reread=reader.GetOutput()
        assert reread.GetNumberOfCells()==len(ids)
        for name,expected in [('FluidMask',fluid),('PortLabel',portlabels),('GlobalLatticeId',ids),('Velocity_m_s',v),('DensityLU',rho),('GaugePressure_Pa',pressure)]:
            assert np.array_equal(vtk_to_numpy(reread.GetCellData().GetArray(name)),expected),f'Roundtrip failure {name}'
        center=vtk.vtkCellCenters();center.SetInputData(reread);center.Update()
        error=float(np.max(np.abs(vtk_to_numpy(center.GetOutput().GetPoints().GetData())-(origin+ijk*dx))))
        assert error<1e-16
        checks.append({'step':step,'file':str(dest),'status':'PASS','vtk_reader_roundtrip':'PASS','cell_center_max_error_m':error,'physical_fluid_cells':int(fluid.sum()),'step2_opened_cells':len(ids),'port_label_counts':{str(i):int((portlabels==i).sum()) for i in range(1,5)},'bytes':dest.stat().st_size})
        collection.append(f'    <DataSet timestep="{step*u["dt_s"]:.17g}" group="" part="0" file="{dest.name}"/>')
        print(f'Export verified: {dest.name}, {len(ids)} exact Step2 cells',flush=True)
    (run/'output/flow.pvd').write_text('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n  <Collection>\n'+'\n'.join(collection)+'\n  </Collection>\n</VTKFile>\n')
    js(run/'diagnostics/paraview_export_check.json',{'status':'PASS' if checks else 'NOT_RUN','format':'VTU subset of original voxel cells, inline binary zlib; no resampling','snapshots':checks})
    js(run/'output/field_units.json',{'coordinate':'m','time':'s','Velocity_m_s':'m/s','VelocityMagnitude_m_s':'m/s','DensityLU':'dimensionless','Density_kg_m3':'kg/m3','GaugePressure_Pa':'Pa, fixed physical gauge zero at rhoLU=1','FluidMask':'1 only for Step2 closed lumen','PortLabel':'exact unchanged Step2 labels 1 inlet, 2/3/4 outlets; all 497 label nodes are native boundary-support markers outside closed lumen','PhysicalFieldValid':'equals FluidMask; all numerical fields at marker-only cells are zero placeholders, not measured ghost values','Step2OpenedEnvelope':'original Step2 opened voxel subset; includes closed lumen + 497 labeled support voxels','GlobalLatticeId':'x fastest; (z*ny+y)*nx+x','rendering':'Threshold FluidMask=1 before viewing velocity, density or pressure. Use separate PortLabel threshold to see unchanged port markers.'})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('run',type=Path);main(parser.parse_args().run)
