import json
import sys
from pathlib import Path

import cv2
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from ycommon import ROOT
from channel import Channel
from microscope_detector import LEGACY,profiled_channel,detect,load_profile
from live_geometry import new_geometry
from protocol import FramePacket
from observer import Observer


def geometry():
    return new_geometry([[20,20],[620,20],[620,460],[20,460]],[640,480])


def pixels():
    a=np.full((480,640),180,np.uint8)
    for y in range(60,430,36):
        for x in range(60,600,36):cv2.circle(a,(x,y),3,30,-1)
    cv2.circle(a,(10,100),3,30,-1)  # Must be excluded by the channel mask.
    return a


def profile():
    return dict(schema='MICROSCOPE_OPENCV_PROFILE_V1',id='fixture',admitted=True,
                camera_image_wh=[640,480],parameters=LEGACY)


def test_profile_legacy_parity_and_no_piv_change():
    a=pixels();original=Channel(geometry()).scan(a);updated=profiled_channel(geometry(),profile()).scan(a)
    for key in ('candidates','areas','count','selected','signal','texture','activity','excluded_candidates'):
        np.testing.assert_array_equal(original[key],updated[key])
    assert np.all(updated['candidates'][:,0]>20)
    assert updated['scan_ms']>0


def test_mask_is_still_required_and_wrong_camera_resolution_rejected(tmp_path):
    cfg=json.loads((ROOT/'configs/default.json').read_text());cfg.update(source_mode='LIVE_CAMERA',backend='B1',live_geometry=None)
    pp=tmp_path/'profile.json';pp.write_text(json.dumps(profile()));cfg['detector_profile']=str(pp)
    observer=Observer(cfg)
    packet=FramePacket('LIVE_CAMERA','software test','test',0,None,pixels(),(640,480))
    assert observer.analyze(packet)['reason']=='MASK_NOT_VERIFIED'
    gp=tmp_path/'geometry.json';gp.write_text(json.dumps(geometry()));cfg['live_geometry']=str(gp)
    assert observer.analyze(packet)['candidates'].shape[0]>0
    wrong=profile();wrong['camera_image_wh']=[1280,1024];pp.write_text(json.dumps(wrong))
    assert Observer(cfg).analyze(packet)['reason']=='DETECTOR_PROFILE_RESOLUTION_MISMATCH'
