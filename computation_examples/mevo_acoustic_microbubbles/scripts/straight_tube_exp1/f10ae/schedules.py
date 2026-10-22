from scripts.straight_tube_exp1.f10ad_multifrequency.interface import FrequencySegment
from mevo_mb.control.scheduled_pzt import Schedule
def schedule(t):
    actions=t.get('segments',[dict(amplitudes=t['amplitudes'],phases=t['phases'],duration=t['duration'])]);segments=[];start=0.
    for i,a in enumerate(actions):
        segments.append(FrequencySegment(i,start,a['duration'],'FOUR_CHANNEL_INDEPENDENT',tuple(a['amplitudes']),tuple(a['phases']),tuple(float(v)>0 for v in a['amplitudes']),frequency_Hz=t.get('frequency',2.5e6)))
        start+=a['duration']
    s=Schedule(tuple(segments))
    if abs(s.duration_s-t['duration'])>1e-15:raise ValueError('Total duration mismatch')
    return s
