import numpy as np
import imageio_ffmpeg
from particle_3d.particle82a_media import decode_video


def test_decode_all_frames_and_reject_incorrect_count(tmp_path):
    import pytest
    file=tmp_path/'moving.mp4'
    writer=imageio_ffmpeg.write_frames(str(file),(96,64),fps=5,codec='libx264',pix_fmt_in='rgb24',pix_fmt_out='yuv420p')
    writer.send(None)
    for i in range(5):
        frame=np.zeros((64,96,3),np.uint8);frame[:,i*12:i*12+16]=[80+30*i,150,200];writer.send(frame)
    writer.close()
    result=decode_video(file,5,tmp_path/'inspection')
    assert result['decoded_frames']==5 and result['distinct_frames']==5
    assert len(result['inspection_frames'])==3
    with pytest.raises(ValueError):decode_video(file,6)
