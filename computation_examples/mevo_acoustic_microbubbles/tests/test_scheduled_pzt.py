import math
from dataclasses import replace
import pytest
from mevo_mb.control.scheduled_pzt import Schedule, Segment, axis_pair, one_hot


def test_legacy_roundtrip_preserves_full_precision():
    e = dict(time_s=0.0, duration_s=.002, amplitude_Pa=[638.9182444862734, 698.8296671687397, 1499.7820148465448, 1203.1590675667037],
             phase_rad=[-.6307367637178037, -.6161282040135343, -1.240646677242141, 2.1722580656754986], enabled=[True]*4, flow_enabled=False)
    s = Segment.from_legacy(e)
    assert Schedule.from_record(Schedule((s,)).record()).segments[0].event() == e


def test_nonaligned_switch_conserves_time_and_dose():
    a = one_hot("X-", 100, duration_s=.0007503)
    b = one_hot("Y+", 250, start_s=.0007503, duration_s=.0012497, segment_id=1)
    steps = list(Schedule((a,b)).steps(1e-6))
    assert len(steps) == 2001
    assert math.fsum(h for t,h,s in steps) == pytest.approx(.002, abs=1e-16)
    assert math.fsum(h for t,h,s in steps if s.segment_id == 0) == pytest.approx(a.duration_s, abs=1e-16)
    assert steps[751][0] == b.start_s
    assert all(t+h <= s.start_s+s.duration_s+1e-16 for t,h,s in steps)
    assert all(abs(steps[i][0]+steps[i][1]-steps[i+1][0])<1e-16 for i in range(len(steps)-1))


def test_aligned_replay_never_changes_original_dt():
    rows = list(Schedule((one_hot("X-"),)).steps(1e-6))
    assert len(rows)==2000 and all(h == 1e-6 for t,h,s in rows)


@pytest.mark.parametrize("update", [dict(frequency_Hz=1e6),dict(amplitude_unit="Vpp"),dict(amplitude_Pa_equivalent=(0,0,1600,0)),dict(phase_rad=(0,0,float('nan'),0)),dict(enabled=(True,True,True,True)),dict(duration_s=-1)])
def test_reject_unsupported_controls(update):
    with pytest.raises(ValueError): replace(one_hot("X-"), **update)


def test_pair_and_schedule_invariants():
    p=axis_pair(x_amplitude=100, x_delta_deg=30)
    assert p.phase_rad[3]-p.phase_rad[2] == math.radians(30)
    with pytest.raises(ValueError): replace(p, amplitude_Pa_equivalent=(0,0,100,101))
    with pytest.raises(ValueError): Schedule((one_hot("X-",start_s=.001),))
    with pytest.raises(ValueError): Schedule((one_hot("X-"),one_hot("Y-",segment_id=1,start_s=.001)))
    with pytest.raises(ValueError): Segment.from_legacy(dict(p.event(),flow_enabled=True))
    with pytest.raises(ValueError): Segment.from_legacy(dict(p.event(),frequency_Hz=1e6))
    with pytest.raises(ValueError): Segment.from_legacy(dict(p.event(),amplitude_unit="Vpp"))


def test_mutable_input_cannot_change_sealed_command():
    a=[0,0,10,0]
    s=Segment(0,0,.002,"ONE_HOT_PZT",a,[0]*4,[False,False,True,False])
    a[2]=2000
    assert s.amplitude_Pa_equivalent[2]==10
