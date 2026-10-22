"""Read-only F10U reuse. No tracker or camera/AFG object is created."""
import sys
from pathlib import Path
V1=Path(__file__).resolve().parents[2]/'collective_motion_v1'
sys.path.insert(0,str(V1/'src'))
from affine import parameters, fit_affine, empty
from distribution import CandidateDetector, distribution
from methods import piv, cpd
from state import ConsumerGate, MockConsumer, make_state, invalidate, tensor_input, CausalWindow
import synthetic as original_synthetic

generator=original_synthetic.generator
case_params=original_synthetic.case_params

def read_frames(indices, video):
    old=original_synthetic.VIDEO
    try:
        original_synthetic.VIDEO=Path(video)
        return original_synthetic.read_frames(indices)
    finally:
        original_synthetic.VIDEO=old
