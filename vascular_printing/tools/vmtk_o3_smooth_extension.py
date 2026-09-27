"""VMTK subprocess worker. Probe and VTP exchange only; never imported by the main runtime.

Unsupported ramp is a hard stop. This worker does not silently select linear/TPS.
"""
import argparse
import hashlib
import importlib.metadata
import inspect
import json
from pathlib import Path
import sys
import traceback


def runtime_probe():
    import vtk
    from vmtk import vtkvmtk, vmtkscripts, vmtkflowextensions, vmtksurfaceendclipper
    flow = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
    scripts = {}
    for module in (vmtkflowextensions, vmtksurfaceendclipper):
        p = Path(inspect.getfile(module))
        scripts[p.name] = {'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
    return dict(python=sys.version, executable=sys.executable,
                vmtk_version=importlib.metadata.version('vmtk'), vtk_version=vtk.vtkVersion.GetVTKVersion(),
                flow_filter=flow.GetClassName(), scripts=scripts,
                flow_methods=[n for n in dir(flow) if any(k in n for k in ('Interpolation', 'Preserve', 'Extension', 'Transition'))],
                ramp_available=hasattr(flow, 'SetInterpolationModeToRamp'),
                preserve_shape_available=hasattr(flow, 'SetPreserveCrossSectionShape'),
                capper_available=hasattr(vmtkscripts, 'vmtkSurfaceCapper'),
                endclipper_available=hasattr(vmtkscripts, 'vmtkSurfaceEndClipper'),
                default_interpolation_enum=flow.GetInterpolationMode())


def roundtrip(request):
    import vtk
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(request['input_vtp'])
    reader.Update()
    surface = reader.GetOutput()
    if not surface.GetNumberOfCells():
        raise ValueError('VMTK_BRIDGE_EMPTY_INPUT')
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(request['output_vtp'])
    writer.SetInputData(surface)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise ValueError('VMTK_BRIDGE_WRITE_FAILED')
    return dict(points=surface.GetNumberOfPoints(), cells=surface.GetNumberOfCells())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--request', required=True, type=Path)
    p.add_argument('--result', required=True, type=Path)
    args = p.parse_args()
    try:
        request = json.loads(args.request.read_text(encoding='utf-8'))
        result = {'runtime': runtime_probe(), 'status': 'PASS', 'operation': request['operation']}
        if request['operation'] == 'roundtrip':
            result['roundtrip'] = roundtrip(request)
        elif request['operation'] != 'probe':
            # No geometry production path exists until the required runtime is available.
            raise ValueError('VMTK_RAMP_UNAVAILABLE' if not result['runtime']['ramp_available']
                             else 'O3_EXTENSION_PIPELINE_NOT_IMPLEMENTED')
        code = 0
    except Exception as exc:
        result = dict(status='FAIL', error=str(exc), traceback=traceback.format_exc())
        code = 2
    args.result.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return code


if __name__ == '__main__':
    sys.exit(main())
