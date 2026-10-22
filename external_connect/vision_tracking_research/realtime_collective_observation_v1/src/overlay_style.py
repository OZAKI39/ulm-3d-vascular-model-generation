"""Presentation only: constant brightness within the original validity window.

The frame match is reported as a delay, never as alternating brightness.
No detections, velocities, validity decisions or time limits are modified.
"""
import math

DEFAULT_LAYERS = dict(Candidates=True, PIV=True, Mask=True, Regions=False,
                      Shape=False, Validity=False, **{'Smooth PIV':True})


def overlay_opacity(age_ms, ttl_ms=250.):
    """Do not pulse brightness during an ordinary analysis update cycle."""
    if age_ms is None or not math.isfinite(age_ms) or ttl_ms <= 0:return 0.
    if age_ms < 0 or age_ms >= ttl_ms:return 0.
    return 1.


def overlay_caption(delay_ms):
    return f'Overlay delay: {max(0.,delay_ms):.0f} ms'


DISPLAY_POLICY = dict(brightness='Constant throughout the original validity window',
    stale_results='Removed at the original TTL',temporal_interpolation=False,
    piv_smoothing='Optional display-only spatial mean and causal temporal filter; raw observations unchanged')
