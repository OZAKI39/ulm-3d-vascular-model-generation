#!/usr/bin/python3
"""Independent VTK readback, physical containment and source immutability audit."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from prepare_port_contract import write_json, EXPECTED_SHA

def snapshot(root):
    out={}
    for base,dirs,files in os.walk(root):
        for name in dirs+files:
            p=Path(base)/name; st=p.lstat()
            out[str(p.relative_to(root))]=[st.st_mode,st.st_size,st.st_mtime_ns,os.readlink(p) if p.is_symlink() else None]
    return out

def verify(run, source):
    from vtkmodules.vtkIOXML import vtkXMLImageDataReader, vtkXMLPolyDataReader
    from vtkmodules.vtkIOGeometry import vtkSTLReader
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkPolyData
    from vtkmodules.vtkFiltersModeling import vtkSelectEnclosedPoints
    from vtkmodules.vtkFiltersCore import vtkImplicitPolyDataDistance
    from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    s=json.loads((run/'diagnostics/voxelization_summary.json').read_text())
    t=json.loads((run/'diagnostics/transform_check.json').read_text())
    shape=tuple(reversed(s['lattice_shape'])); dx=s['effective_dx_m'];origin=np.array(s['physical_origin_m'])
    checks=[]
    for filename,field,raw in [('closed_flag_matrix','ClosedFluid','closed_flag_matrix'),('opened_flag_matrix','OpenedFluid','opened_flag_matrix'),('port_label_field','PortLabel','port_label_field')]:
        reader=vtkXMLImageDataReader();reader.SetFileName(str(run/(filename+'.vti')));reader.Update();data=reader.GetOutput()
        arr=vtk_to_numpy(data.GetCellData().GetArray(field)).reshape(shape)
        reference=np.memmap(run/'diagnostics'/(raw+'.u8'),dtype=np.uint8,mode='r',shape=shape)
        check=dict(file=filename+'.vti',readback_matches_raw=bool(np.array_equal(arr,reference)),
            cell_count=data.GetNumberOfCells(),dimensions=list(data.GetDimensions()),
            origin_m=list(data.GetOrigin()),spacing_m=list(data.GetSpacing()),
            label_histogram=np.bincount(arr.ravel()).tolist())
        check['pass']=bool(check['readback_matches_raw'] and np.allclose(data.GetOrigin(),origin-.5*dx,rtol=0,atol=1e-18) and
            np.allclose(data.GetSpacing(),[dx]*3,rtol=0,atol=1e-20) and np.array_equal(data.GetDimensions(),np.array(s['lattice_shape'])+1))
        checks.append(check)
        del arr,reference,data,reader
    reader=vtkXMLPolyDataReader();reader.SetFileName(str(run/'cap_triangles.vtp'));reader.Update()
    poly=reader.GetOutput(); cap_labels=vtk_to_numpy(poly.GetCellData().GetArray('PortLabel'))
    cap_ok=bool(poly.GetNumberOfCells()==sum(p['reconstructed_triangle_count'] for p in s['ports'].values()) and
        all(np.count_nonzero(cap_labels==p['label'])==p['reconstructed_triangle_count'] for p in s['ports'].values()))
    # A separate VTK enclosed-point method checks the original, uninflated physical STL.
    closed=np.memmap(run/'diagnostics/closed_flag_matrix.u8',dtype=np.uint8,mode='r',shape=shape)
    ijk=np.argwhere(closed!=0)[:,::-1]; xyz=origin+ijk*dx
    points=vtkPoints();points.SetData(numpy_to_vtk(xyz,deep=True));cloud=vtkPolyData();cloud.SetPoints(points)
    stl_reader=vtkSTLReader();stl_reader.SetFileName(s['input_stl']);stl_reader.MergingOn();stl_reader.Update()
    enclosed=vtkSelectEnclosedPoints();enclosed.SetInputData(cloud);enclosed.SetSurfaceData(stl_reader.GetOutput())
    enclosed.SetTolerance(1e-9);enclosed.CheckSurfaceOn();enclosed.Update()
    inside=vtk_to_numpy(enclosed.GetOutput().GetPointData().GetArray('SelectedPoints'))
    exceptions=xyz[inside==0]
    distances=[]
    if len(exceptions):
        distance=vtkImplicitPolyDataDistance();distance.SetInput(stl_reader.GetOutput())
        distances=[abs(distance.EvaluateFunction(point)) for point in exceptions]
    tolerance=s['inflate_lu']*dx+1e-14
    containment=dict(method='VTK vtkSelectEnclosedPoints, original frozen uninflated meter STL',
        closed_fluid_points_checked=len(xyz),points_reported_outside=len(exceptions),
        max_outside_distance_m=max(distances,default=0),expected_inflation_bound_m=tolerance,
        pass_with_recorded_inflation=bool(max(distances,default=0)<=tolerance),
        note='An inflated Palabos surface may include nodes within 0.001 LU of the original surface; no input geometry is edited')
    write_json(run/'diagnostics/vtk_artifact_verification.json',dict(images=checks,cap_vtp_pass=cap_ok,independent_containment=containment))
    print('VTK_READBACK',all(x['pass'] for x in checks),'CONTAINMENT',containment,flush=True)
    hc=source.parent.parent
    tracked=json.loads((run/'diagnostics/hemocell_tracked_before.json').read_text()); bad=[]
    for name,digest in tracked.items():
        p=hc/name;actual=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
        if actual!=digest:bad.append(name)
    comparisons={}
    roots={'hemocell':hc,'step1':Path('/home/lzy/projects/compre_output/step1/20260912_215759'),
           'old_project':Path('/home/lzy/projects/ulm_3D_vascular')}
    allowed_prefix='test_code/step2_vascular_voxelization_smoke'
    for key,root in roots.items():
        before=json.loads((run/'diagnostics'/f'{key}_before_manifest.json').read_text());after=snapshot(root)
        changed=[p for p in before if p not in after or before[p]!=after[p]]
        # git status can create/remove an index lock, updating the .git DIRECTORY
        # mtime without changing any contained file. Preserve and distinguish it.
        git_directory_metadata=[p for p in changed if key=='hemocell' and p=='.git' and
            p in after and before[p][:2]==after[p][:2] and before[p][3]==after[p][3]]
        content_or_other_changes=[p for p in changed if p not in git_directory_metadata]
        added=sorted(set(after)-set(before))
        unexpected=[p for p in added if not (key=='hemocell' and (p=='test_code' or p==allowed_prefix or p.startswith(allowed_prefix+'/')))]
        comparisons[key]=dict(preexisting_entries=len(before),changed_or_missing=changed,
            git_directory_mtime_only_changes=git_directory_metadata,
            new_entries=added,unexpected_new_entries=unexpected,pass_check=not content_or_other_changes and not unexpected)
    step1_hashes=json.loads((run/'diagnostics/step1_hashes_before.json').read_text())
    step1_bad=[name for name,digest in step1_hashes.items() if hashlib.sha256((roots['step1']/name).read_bytes()).hexdigest()!=digest]
    input_ok=hashlib.sha256(Path(s['input_stl']).read_bytes()).hexdigest()==EXPECTED_SHA
    codes=json.loads((run/'logs/return_codes.json').read_text())
    text=''
    for args in [['status','--short'],['rev-parse','HEAD'],['diff','--stat']]:
        text+='git '+' '.join(args)+'\n'+subprocess.check_output(['git','-C',str(hc),*args],text=True)
    (run/'git_status_after.txt').write_text(text)
    audit=dict(existing_tracked_hash_count=len(tracked),existing_tracked_hash_mismatches=bad,
        tree_comparisons=comparisons,step1_content_hashes_checked=len(step1_hashes),step1_content_hash_mismatches=step1_bad,
        frozen_stl_integrity=input_ok,hemocell_core_modified=bool(bad),palabos_modified=any(x.startswith('palabos/') for x in bad),
        existing_tracked_file_modified_by_step2=bool(bad))
    audit['pass']=bool(not bad and not step1_bad and input_ok and all(c['pass_check'] for c in comparisons.values()))
    write_json(run/'diagnostics/source_immutability.json',audit)
    static={'cpp_files':[str(p) for p in source.glob('*.cpp')],
        'contains_time_step_call':any('collideAndStream(' in p.read_text() for p in source.glob('*.cpp')),
        'contains_official_helper_call':any('getFlagMatrixFromSTL(' in p.read_text() for p in source.glob('*.cpp')),
        'compiled_executable_sha256':hashlib.sha256((run/'build/closed_voxelizer').read_bytes()).hexdigest()}
    write_json(run/'diagnostics/execution_scope_check.json',static)
    all_pass=bool(s['step2_auto_check']=='PASS' and t['status']=='PASS' and all(x['pass'] for x in checks) and cap_ok and
        containment['pass_with_recorded_inflation'] and audit['pass'] and
        all(codes[k]==0 for k in ['CONFIGURE_RC','BUILD_RC','SMOKE_RC']) and
        not static['contains_time_step_call'] and not static['contains_official_helper_call'])
    s.update(step2_auto_check='PASS' if all_pass else 'FAIL',step2_status='AUTO_PASS_HUMAN_PENDING' if all_pass else 'PARTIAL',
        independent_artifact_and_source_audit='PASS' if all_pass else 'FAIL')
    write_json(run/'diagnostics/voxelization_summary.json',s)
    print('FINAL_AUTOMATED_AUDIT',all_pass,flush=True)
    assert all_pass,'Artifact/source audit failed; inspect saved diagnostics'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();verify(a.run,Path(__file__).resolve().parent)
