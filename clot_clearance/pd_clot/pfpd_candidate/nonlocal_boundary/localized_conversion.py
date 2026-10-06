"""Frozen boundary conversion into an independently evaluated bulk operator."""
import numpy as np


def localize(f_FD, q_bulk, b_local, volume):
    V = np.asarray(volume)[:, None]
    fbulk = np.asarray(q_bulk)*V; body = np.asarray(b_local)*V
    model = (np.asarray(q_bulk)+np.asarray(b_local))*V
    localized = fbulk-np.asarray(f_FD)
    return dict(f_FD=np.asarray(f_FD), f_bulk_target=fbulk, f_BC_localized=localized,
                f_BC_naive=-np.asarray(f_FD), f_body_local=body, f_model=model,
                f_external_NL=-fbulk+localized, f_external_classical=body+localized)


def convert_localized(mechanics, x, q_bulk, b_local, volume):
    """Backend-independent entry point; the accepted decomposition is unchanged."""
    return localize(mechanics.force(x), q_bulk, b_local, volume)
