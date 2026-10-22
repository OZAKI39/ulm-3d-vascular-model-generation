"""Audit installed modules and read existing data with old workspaces denied.

Run with this project's Python -I. No source lookup path is added at runtime.
"""
import ast
import argparse
import importlib
import json
import os
from pathlib import Path
import sys
import traceback
from human_vessel_flow.paths import ROOT, sha256, write_json, fresh_directory

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', nargs='?', type=Path,
                        default=ROOT/'validation/independence_check')
    output = fresh_directory(parser.parse_args().output)
    denied_names = ('computation_examples', 'vascular_printing', 'new_mouse_vessel_flow',
                    'temp_storage', 'ulm_particle_formal_p9a5', 'brava_flow_roi_18mlmin',
                    'brava_calibration_equal_split', 'sonovue_size_distribution_v0',
                    'github_sync_audits', 'new_clot_simulation', 'new_vessel_printing',
                    'new_microbubble_experiment', 'new_microbubble_aggregation')
    denied = [(ROOT.parent/name).absolute() for name in denied_names]
    attempted = []
    opened = set()
    def audit(event, args):
        if event not in ('open', 'os.listdir', 'os.scandir') or not args or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        p = Path(os.path.realpath(os.fsdecode(args[0])))
        if any(p.is_relative_to(d) for d in denied):
            attempted.append(str(p))
            raise RuntimeError('Old workspace access forbidden in independence check: ' + str(p))
        opened.add(str(p))
    sys.addaudithook(audit)
    errors=[]; imported=[]; graph={}; findings=[]
    for p in sorted((ROOT/'src').rglob('*.py')):
        relative=p.relative_to(ROOT/'src').with_suffix('')
        parts=list(relative.parts)
        if parts[-1]=='__init__':parts.pop()
        name='.'.join(parts)
        tree=ast.parse(p.read_text()); imports=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):imports.extend(a.name for a in node.names)
            if isinstance(node,ast.ImportFrom):imports.append('.'*node.level+(node.module or ''))
            if isinstance(node,ast.Call):
                call=ast.unparse(node.func)
                if call in {'exec','eval','compile','importlib.util.spec_from_file_location','runpy.run_path','sys.path.insert','sys.path.append'}:
                    findings.append(dict(file=str(p.relative_to(ROOT)),line=node.lineno,call=call))
        graph[name]=sorted(set(imports))
        try:
            module=importlib.import_module(name)
            actual=Path(module.__file__).resolve()
            assert actual==p.resolve(),(name,actual,p)
            imported.append(dict(module=name,file=str(p.relative_to(ROOT)),sha256=sha256(p)))
        except Exception:
            errors.append(dict(module=name,error=traceback.format_exc()))
    # Exercise actual data loading, not only syntax/import checks.
    try:
        from human_vessel_flow.geometry import verify_geometry
        from human_vessel_flow.flow_field import load
        from human_vessel_flow.particles import load_environment, historical_summary
        from human_vessel_flow.reference_pipe.review import review
        geometry=verify_geometry();field=load();env=load_environment();population=historical_summary()
        checks=dict(geometry=geometry['status'],flow_nodes=field.mesh.n_points,
                    environment_flow_sha256=env.flow_sha256,population=population)
    except Exception:
        errors.append(dict(module='actual_data_loading',error=traceback.format_exc()));checks={}
    links=[]
    for directory, dirs, files in os.walk(ROOT, followlinks=False):
        dirs[:] = [d for d in dirs if d not in ('.venv', '.git', '__pycache__', '.pytest_cache')]
        for name in dirs + files:
            p = Path(directory)/name
            if p.is_symlink():links.append(dict(path=str(p.relative_to(ROOT)),target=os.readlink(p)))
    info=dict(status='PASS' if not errors and not attempted and not findings and not links else 'FAIL',
              python=sys.executable,isolated_flag=sys.flags.isolated,
              modules=imported,module_count=len(imported),errors=errors,
              forbidden_source_loading=findings,denied_workspace_accesses=attempted,
              project_symlinks_excluding_standard_venv=links,data_checks=checks,
              opened_project_paths=sorted(p for p in opened if Path(p).is_relative_to(ROOT)),
              source_graph='SOURCE_IMPORT_GRAPH.json',
              informational_old_path_literals=['src/particle_3d/audit.py: legacy report metadata only; used=False, copied_2d_code=False'],
              native_dependencies='Explicit configs/cfd/native.json; separate unified solver installation')
    write_json(output/'SOURCE_IMPORT_GRAPH.json',graph)
    write_json(output/'INDEPENDENCE.json',info)
    print(json.dumps({k:v for k,v in info.items() if k not in {'modules','opened_project_paths'}},indent=2))
    if info['status']!='PASS':raise SystemExit(1)

if __name__=='__main__':main()
