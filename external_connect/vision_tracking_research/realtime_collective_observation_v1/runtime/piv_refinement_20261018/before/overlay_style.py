"""Presentation only: stable overlay brightness with a bounded expiry fade.

The frame match is reported as a delay, never as alternating brightness.
No detections, velocities, validity decisions or time limits are modified.
"""
import math

DEFAULT_LAYERS = dict(Candidates=True, PIV=True, Mask=True, Regions=False,
                      Shape=False, Validity=False)
FADE_MS = 80.


def overlay_opacity(age_ms, ttl_ms=250.):
    """Full contrast until the last 80 ms of validity; zero at expiry."""
    if age_ms is None or not math.isfinite(age_ms) or ttl_ms <= 0:return 0.
    if age_ms < 0 or age_ms >= ttl_ms:return 0.
    fade_ms = min(FADE_MS, ttl_ms)
    x = min(1., max(0., (ttl_ms-age_ms)/fade_ms))
    return x*x*(3.-2.*x)


def overlay_caption(delay_ms):
    return f'Overlay delay: {max(0.,delay_ms):.0f} ms'


DISPLAY_POLICY = dict(brightness='Constant for exact and recent frames; fade near expiry',
    fade_before_expiry_ms=FADE_MS,stale_results='Removed at the original TTL',
    temporal_interpolation=False)
