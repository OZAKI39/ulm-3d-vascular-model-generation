"""File-based bridge to the existing Windows VMTK; no VMTK imports in this process."""
from pathlib import Path
import json
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def wsl_path(path, windows=False):
    return subprocess.check_output(['wslpath', '-w' if windows else '-u', str(path)], text=True).strip()


def invoke(executable, output, operation='probe', input_vtp=None, timeout=60):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    temp_probe = subprocess.run([str(executable), '-c', 'import tempfile; print(tempfile.gettempdir())'],
                               capture_output=True, text=True, timeout=timeout, check=True)
    windows_temp = Path(wsl_path(temp_probe.stdout.strip().splitlines()[-1]))
    with tempfile.TemporaryDirectory(prefix='ulm_o3_vmtk_', dir=windows_temp) as temporary:
        stage = Path(temporary)
        shutil.copy2(ROOT / 'tools/vmtk_o3_smooth_extension.py', stage / 'worker.py')
        request = {'operation': operation}
        if input_vtp is not None:
            shutil.copy2(input_vtp, stage / 'input.vtp')
            request.update(input_vtp=wsl_path(stage / 'input.vtp', True), output_vtp=wsl_path(stage / 'output.vtp', True))
        (stage / 'request.json').write_text(json.dumps(request), encoding='utf-8')
        command = [str(executable), wsl_path(stage / 'worker.py', True), '--request',
                   wsl_path(stage / 'request.json', True), '--result', wsl_path(stage / 'result.json', True)]
        run = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=timeout)
        (output / 'stdout.log').write_text(run.stdout)
        (output / 'stderr.log').write_text(run.stderr)
        result_path = stage / 'result.json'
        result = json.loads(result_path.read_text()) if result_path.exists() else {'status': 'FAIL', 'error': 'NO_WORKER_RESULT'}
        result.update(command=command, request=request, returncode=run.returncode,
                      exchange_method='Windows temporary directory, mapped with wslpath; copied VTP only')
        if (stage / 'output.vtp').exists():
            shutil.copy2(stage / 'output.vtp', output / 'roundtrip_output.vtp')
        (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def require_ramp(runtime):
    if not runtime.get('ramp_available'):
        raise RuntimeError('VMTK_RAMP_UNAVAILABLE')
    if not runtime.get('preserve_shape_available'):
        raise RuntimeError('VMTK_PRESERVE_SHAPE_API_UNAVAILABLE')
    if not runtime.get('capper_available'):
        raise RuntimeError('VMTK_CAPPER_UNAVAILABLE')


def probe_slicer(executable, output, timeout=60):
    """Launch only our own Slicer probe and exit it; no extension installation or scene loading."""
    code = ('import json,importlib.util,slicer; e=slicer.app.extensionsManagerModel(); '
            'print("O3_SLICER_PROBE="+json.dumps(dict(slicer_version=slicer.app.applicationVersion, '
            'revision=slicer.app.repositoryRevision,extensions_install_path=slicer.app.extensionsInstallPath, '
            'installed_extensions=list(e.installedExtensions), '
            'vmtk_spec=str(importlib.util.find_spec("vmtk")), '
            'vtkvmtk_spec=str(importlib.util.find_spec("vtkvmtkComputationalGeometryPython")), '
            'clip_vessel_spec=str(importlib.util.find_spec("ClipVessel"))))); slicer.app.exit(0)')
    command = [str(executable), '--no-splash', '--no-main-window', '--ignore-slicerrc', '--python-code', code]
    run = subprocess.run(command, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=timeout)
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    (output / 'stdout.log').write_text(run.stdout); (output / 'stderr.log').write_text(run.stderr)
    payloads = [s.split('O3_SLICER_PROBE=', 1)[1] for s in run.stdout.splitlines() if 'O3_SLICER_PROBE=' in s]
    result = json.loads(payloads[-1]) if payloads else {'error': 'NO_SLICER_PROBE_PAYLOAD'}
    result.update(command=command, returncode=run.returncode)
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result
