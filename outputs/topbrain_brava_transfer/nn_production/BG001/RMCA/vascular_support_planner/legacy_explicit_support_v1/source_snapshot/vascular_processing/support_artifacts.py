"""Reproducibility records for local s1-8 outputs, without printer operations."""
from __future__ import annotations
import importlib.metadata
import json
import platform
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from .support_inputs import ROOT,digest,save

def write_manifests(out,summary):
    out=Path(out).resolve()
    sources=set((ROOT/'vascular_processing').glob('support_*.py'))
    sources.update((ROOT/'vascular_processing').glob('vascular_support_pipeline.py'))
    sources.update((ROOT/'tests/data/support').glob('*'))
    for name in ['s1-8_vascular_support_planner.py','tests/test_s1_8_vascular_support_planner.py',
                 'config/vascular_support_planner.yaml','tools/prepare_support_engines.py','tools/verify_s1_8_support_outputs.py',
                 'third_party/support_sources.json','open_source_support_audit.md',
                 'docs/VASCULAR_SUPPORT_PLANNER使用说明.md']:
        sources.add(ROOT/name)
    sources.update(ROOT.glob('requirements*.txt'))
    versions={}
    for name in ['numpy','scipy','trimesh','shapely','manifold3d','python-fcl','embreex','pyvista','vtk','pytest']:
        try:versions[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:versions[name]=None
    test_summary=None
    testfile=out/'test_results.xml'
    if testfile.exists():
        suites=list(ET.parse(testfile).getroot().iter('testsuite'))
        test_summary={k:sum(int(s.get(k,0)) for s in suites) for k in ['tests','failures','errors','skipped']}
    save(out/'implementation_manifest.json',dict(architecture=summary['architecture'],python=platform.python_version(),
        dependencies=versions,tests=test_summary,files=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p)) for p in sorted(sources) if p.is_file()],
        upstream_geometry_modified=False,benchmarks_role='BENCHMARK_ONLY'))
    files=set()
    for directory in ['geometry','tables','coupon','custom_candidates']:
        files.update(p for p in (out/directory).rglob('*') if p.is_file())
    visual=out/'visual_qc_manifest.json'
    if visual.exists():files.update(Path(p) for p in json.loads(visual.read_text())['main_images'])
    files.update(out.glob('custom_*_auto_support_off.3mf'))
    for name in ['vascular_support_summary.json','vascular_support_report.md','implementation_manifest.json',
                 'support_graph.json','support_detection_config.json','source_protection.json',
                 'custom_search.json','test_results.xml','test_run.log','visual_qc_manifest.json','verification.json']:
        if (out/name).is_file():files.add(out/name)
    save(out/'artifact_manifest.json',dict(status=summary['status'],production_selected=summary['production_selected'],
        diagnostic_candidate=summary['diagnostic_candidate'],scale_factor=1.,
        artifacts=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(files) if p.is_file()]))

def verify_3mf(path):
    """Read the completed local archive; a ZIP filename alone is not success."""
    path=Path(path)
    with zipfile.ZipFile(path) as archive:
        bad=archive.testzip()
        settings=json.loads(archive.read('Metadata/project_settings.config'))
        return dict(path=str(path),zip_crc_passed=bad is None,bad_member=bad,
            gcode_members=[n for n in archive.namelist() if n.endswith('.gcode')],
            actual_settings={k:settings.get(k) for k in ['printer_model','nozzle_diameter','layer_height','line_width','filament_type','enable_support']},
            sha256=digest(path))
