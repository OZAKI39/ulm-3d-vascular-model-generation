import json
from pathlib import Path
import sys

import cv2
import numpy as np
import pytest

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'src'))
sys.path.insert(0,str(BASE.parent/'realtime_collective_observation_v1/src'))
from microscope_detector import detect, load_profile
from vision_evaluation import evaluate_instances


def test_small_dim_objects_and_artifacts_are_separated():
    image=np.full((100,160),180,np.uint8)
    cv2.circle(image,(25,35),4,100,-1)
    cv2.circle(image,(75,35),4,170,-1)
    cv2.rectangle(image,(110,20),(112,75),50,-1)
    _,_,objects=detect(image)
    assert len(objects)==1
    assert np.linalg.norm(np.asarray(objects[0]['xy_local'])-[25,35])<1
    _,_,sensitive=detect(image,dict(threshold=3))
    assert len(sensitive)==2


def test_explicit_profile_admission_and_parameter_validation(tmp_path):
    p=tmp_path/'profile.json'
    p.write_text(json.dumps(dict(schema='MICROSCOPE_OPENCV_PROFILE_V1',admitted=False,parameters={})))
    with pytest.raises(ValueError):load_profile(p)
    p.write_text(json.dumps(dict(schema='MICROSCOPE_OPENCV_PROFILE_V1',admitted=True,parameters={'threshold':float('nan')})))
    with pytest.raises(ValueError):load_profile(p)


def test_large_false_positive_is_counted_in_new_study():
    gt=np.zeros((100,100),np.uint32);pred=gt.copy();pred[20:50,20:50]=1
    objects=[dict(label=1,area=900,xy_local=[34.5,34.5])]
    assert evaluate_instances(gt,pred,objects,max_prediction_area=None)['center']['fp']==1
    # Historical scoring stays reproducible for old callers.
    assert evaluate_instances(gt,pred,objects)['center']['fp']==0


def test_cellpose_checkpoint_selection_does_not_overwrite_original(tmp_path):
    import hashlib
    from cellpose_adaptation import checkpoint_for
    original=tmp_path/'models/cellpose/cyto3';original.parent.mkdir(parents=True)
    original.write_bytes(b'original');(tmp_path/'data').mkdir()
    assert checkpoint_for(tmp_path)['id']=='ORIGINAL_CYTO3'
    tuned=tmp_path/'tuned';tuned.write_bytes(b'tuned')
    profile=dict(schema='MICROSCOPE_CELLPOSE_PROFILE_V1',admitted=True,id='fixture',checkpoint=str(tuned),
                 checkpoint_sha256=hashlib.sha256(b'tuned').hexdigest())
    (tmp_path/'data/cellpose_profile.json').write_text(json.dumps(profile))
    assert checkpoint_for(tmp_path)['id']=='fixture'
    assert checkpoint_for(tmp_path,force_original=True)['id']=='ORIGINAL_CYTO3'
    assert original.read_bytes()==b'original'
    tuned.write_bytes(b'tampered')
    with pytest.raises(ValueError,match='checksum'):checkpoint_for(tmp_path)
