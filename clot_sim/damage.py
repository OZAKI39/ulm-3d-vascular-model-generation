"""Read-only history/softening predictor; native code commits accepted damage."""

import numpy as np


def predict(model, ref, x, k, onset, off):
    e, _ = model.geometry(x)
    kk = np.maximum(k, np.maximum(e / ref["r0"], 0))
    d = np.zeros(len(kk))
    d[kk >= off] = 1
    part = (kk > onset) & (kk < off)
    d[part] = (
        off[part] * (kk[part] - onset[part]) / (kk[part] * (off[part] - onset[part]))
    )
    return (kk, d)
