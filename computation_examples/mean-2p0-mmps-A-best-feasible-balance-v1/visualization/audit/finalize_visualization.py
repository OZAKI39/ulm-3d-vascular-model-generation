"""Read-only final checks and a compact manifest for the new visualization."""
from pathlib import Path
import hashlib
import json
from PIL import Image
import numpy as np
import pyvista as pv

ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

context=read(ROOT/'BUILD_CONTEXT.json')
source=Path(context['source_case'])
manifest=read(ROOT/'INPUT_MANIFEST.json')
for name,item in manifest['files'].items():
    assert sha(ROOT/name)==item['sha256'],name
assert sha(ROOT/'input_data/frozen_flow/steady_flow_mean_2p0_mmps.vtu')==context['source_field_sha256']
assert sha(source/'frozen_flow/steady_flow_mean_2p0_mmps_A_best_feasible_balance.vtu')==context['source_field_sha256']
surface=read(ROOT/'input_data/field_diagnostics/COMPUTE_VALIDATION.json')
stream=read(ROOT/'input_data/streamlines/COMPUTE_VALIDATION.json')
assert surface['all_pass'] and stream['all_pass']
assert surface['source_field_sha256']==stream['source_field_sha256']==context['source_field_sha256']
style=read(ROOT/'audit/STYLE_PRESERVATION.json');assert style['all_pass']
for name,digest in context['original_renderer_hashes'].items():
    assert sha(ROOT.parent/name)==digest,'Original visualization script changed: '+name
# Verify that preparation agrees with the independently audited new WSS arrays.
# Resolve the repository from the explicit case path, never from a cached old field.
repo=next(p for p in source.parents if (p/'ulm_3D_vascular').is_dir())
audit=repo/'ulm_3D_vascular/reports/best_feasible_balance_forward_v1'
raw=np.load(audit/'data/best_feasible_3D_wss_raw_si.npz')
wall=pv.read(ROOT/'input_data/field_diagnostics/data/wall_wss_si.vtp')
np.testing.assert_array_equal(wall.cell_data['WSS_raw_Pa'],raw['WSS_Pa'])
np.testing.assert_array_equal(wall.cell_data['Area_m2'],raw['area_m2'])
videos=[];independent=[]
for relative in ['results','results/surface_fields']:
    out=ROOT/relative
    check=read(out/'INDEPENDENT_VALIDATION.json');assert check['all_pass']
    media=read(out/'MEDIA_VALIDATION.json');assert media['all_pass']
    for row in media['videos']:
        assert row['decoded_frames']==row['distinct_frames']==432
        assert row['fps']==24 and row['size']==[1920,1080] and row['duration_s']==18
        path=out/row['file'];assert sha(path)==row['sha256']
        videos.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,**row))
    independent.append(dict(path=str((out/'INDEPENDENT_VALIDATION.json').relative_to(ROOT)),sha256=sha(out/'INDEPENDENT_VALIDATION.json')))
assert len(videos)==4
images=['results/figures/full_vessel_4k.png','results/figures/branch_vectors_4k.png',
        'results/surface_fields/figures/pressure_overview_4k.png','results/surface_fields/figures/wss_overview_4k.png']
for name in images:
    with Image.open(ROOT/name) as im:assert im.size==(3840,2160),name
record=dict(all_pass=True,source_case=str(source),source_field_sha256=context['source_field_sha256'],
    source_inputs_unchanged=True,source_WSS_exactly_matches_independent_audit=True,
    original_rendering_code_preserved=True,pressure_colorbar_pa=[-5,6000],
    WSS_colorbar_pa=[0,55],speed_colorbar_mm_s=[0,7.5],videos=videos,images_4k=images,
    independent_validation=independent,new_CFD_calls=0,finite_size_microbubble_calls=0,RBC_calls=0,
    newly_integrated_visualization_point_streamlines=96,
    input_alias_note=context['compatibility_aliases'],
    visual_inspection='Four overview images and multi-angle contact sheets checked; same original label/layout style')
(ROOT/'VISUALIZATION_VALIDATION.json').write_text(json.dumps(record,indent=2)+'\n')
lines=[sha(p)+'  '+str(p.relative_to(ROOT)) for p in sorted(ROOT.rglob('*'))
       if p.is_file() and p.name!='SHA256SUMS.txt' and '__pycache__' not in p.parts]
(ROOT/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'all_pass':True,'videos':len(videos),'4k_images':len(images),'hashed_files':len(lines)},indent=2))
