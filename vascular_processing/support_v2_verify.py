"""Final artifact verification, separate from scientific/physical acceptance."""
import json,zipfile
from pathlib import Path
import numpy as np
import trimesh
from .support_inputs import save,verify,digest,ROOT
from .support_v2_orientation import rigid

def run(out):
    out=Path(out);summary=json.loads((out/'vascular_support_v2_summary.json').read_text());rows=[]
    source=verify(json.loads((out/'protected_before.json').read_text()))
    required=['custom_support.stl','mold_without_support.stl','mold_with_custom_support.stl','support_roofs.stl','support_trunks.stl','support_breakaway_parts.stl']
    for identity in summary['top4']:
        folder=out/'candidates'/identity;r=json.loads((folder/'candidate.json').read_text());m=np.array(r['orientation']['print_transform']);inv=np.array(r['orientation']['casting_restore_transform']);meshes=[];archives=[]
        for name in required:
            path=folder/name;mesh=trimesh.load_mesh(path) if path.exists() else None
            meshes.append(dict(file=name,exists=path.exists(),sha256=digest(path) if path.exists() else None,closed_positive_volume=bool(mesh is not None and mesh.is_volume),zero_area_faces=int((mesh.area_faces==0).sum()) if mesh is not None else None))
        for path in folder.rglob('*.3mf'):
            with zipfile.ZipFile(path) as archive:
                corrupt=archive.testzip();settings=json.loads(archive.read('Metadata/project_settings.config')) if 'Metadata/project_settings.config' in archive.namelist() else {}
                archives.append(dict(path=str(path),crc_passed=corrupt is None,auto_support_off=str(settings.get('enable_support')) in ['0','false','False']))
        rows.append(dict(candidate_id=identity,rigid_transform=rigid(m),inverse_verified=bool(np.allclose(m@inv,np.eye(4),atol=1e-10)),mesh_exports=meshes,archives=archives,bambu_slicer_success=r.get('slice',{}).get('slicer_success',False),hard_gate_failures=r.get('hard_gate_failures',[])))
    images=[out/'QC'/f'{i+1:02}_{stem}.png' for i,stem in enumerate(__import__('vascular_processing.support_v2_review',fromlist=['STEMS']).STEMS)]
    code=[ROOT/'s1-8_vascular_support_planner.py',*sorted((ROOT/'vascular_processing').glob('support_v2_*.py'))]
    result=dict(source_protection=source,architecture=summary['architecture'],status=summary['status'],production_selected=summary['production_selected'],coarse_count=summary['coarse_count'],full_planning_count=summary['full_planning_count'],top4_count=len(rows),candidates=rows,all_20_qc_images_exist=all(p.exists() for p in images),code_sha256={str(p):digest(p) for p in code},physical_print_validated=False,printer_upload=False,printer_start=False)
    save(out/'artifact_verification.json',result);return result
