"""Run after the new solve/export/render/audit; missing deliverables are failures."""
from pathlib import Path
import sys,json
import numpy as np
import pyvista as pv
import pytest
from PIL import Image
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/flow_2mmps'))
from validate import CASE,Q,sha,check_mass,mesh_matches
from audit_outputs import image_check,decode_video
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv11 import linear_gate,nonlinear_gate

def read(name):return json.loads((CASE/'reports'/name).read_text())

def test_native_solve_and_steady_evidence():
    execution=read('execution.json');assert execution['status']=='PASS'
    assert linear_gate(execution) and nonlinear_gate(execution['history'])
    physics=read('physics_validation.json');assert physics['status']=='PASS'
    assert len(physics['steady']['intervals'])>=5
    for row in physics['steady']['intervals'][-5:]:
        assert row['E_u']<=1e-5 and row['E_Q']<=1e-6
    assert physics['native_checkpoint']['complete_integration_history']

def test_new_export_is_actual_native_solution_and_original_mesh():
    frozen=CASE/'frozen_flow';manifest=json.loads((frozen/'manifest.json').read_text())
    for name,record in manifest['files'].items():assert sha(frozen/name)==record['sha256']
    assert sha(CASE/manifest['raw_solver_vtu'])==sha(frozen/'steady_flow_mean_2p0_mmps.vtu')
    measure=SolutionMeasurements(CASE/'SV_MESH/mesh_arrays.npz',Q,.002)
    grid=pv.read(frozen/'steady_flow_mean_2p0_mmps.vtu');assert mesh_matches(grid,measure)
    u,p=measure.read(frozen/'steady_flow_mean_2p0_mmps.vtu');assert check_mass(measure.measure(u,p))
    with np.load(frozen/'flow_arrays_si.npz') as a:
        assert np.array_equal(a['velocity_m_s'],u) and np.array_equal(a['pressure_pa'],p)
        assert np.array_equal(a['points_m'],measure.points)
    assert measure.measure(u,p)['wall_noslip_pass']

def test_old_frozen_artifacts_unchanged():
    for name,expected in read('old_baseline_hashes.json').items():
        assert sha(ROOT/'frozen_reference'/name)==expected

def test_all_figures_are_readable_black_background_files():
    audit=read('artifact_validation.json');assert audit['status']=='PASS'
    assert len(audit['figures'])>=6
    for item in audit['figures']:
        p=CASE/item['file'];assert sha(p)==item['sha256'];assert image_check(p)['status']=='PASS'

def test_both_animations_have_full_decode_evidence():
    audit=read('artifact_validation.json');assert len(audit['animations'])==2
    for item in audit['animations']:
        assert sha(CASE/item['file'])==item['sha256']
        assert item['status']=='PASS' and item['decoded_frames']==432
        assert item['decoder_exit_code']==0 and item['decoder_stderr']==''
        assert item['metadata']['fps']==24
    assert audit['storyboard']['source']=='actual decoded MP4 frames'

def test_blank_figure_is_rejected(tmp_path):
    p=tmp_path/'blank.png';Image.new('RGB',(1920,1080),'black').save(p)
    with pytest.raises(AssertionError):image_check(p)

def test_truncated_animation_is_rejected(tmp_path):
    source=CASE/'animations/Animation_01_full_vessel.mp4'
    path=tmp_path/'truncated.mp4';path.write_bytes(source.read_bytes()[:200])
    with pytest.raises(AssertionError):decode_video(path)

def test_cached_frame_video_is_rejected():
    path=CASE/'reports/render_attempt_01/Animation_01_full_vessel.mp4'
    with pytest.raises(AssertionError,match='No visible camera rotation'):
        decode_video(path)
