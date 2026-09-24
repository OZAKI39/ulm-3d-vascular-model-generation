import json

def test_fresh_render_windows_preserve_saved_state_and_pixels(saved_scene):
    result=json.loads((saved_scene.root/'data/render_determinism.json').read_text())
    assert result['all_pass'] and result['pixel_exact'] and result['frame_receipt_exact']
    assert result['saved_state_unmodified'] and result['input_arrays_readonly']
    assert result['source_scene_sha256']==saved_scene.sha256
