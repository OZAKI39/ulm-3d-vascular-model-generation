# Bounded AI flow comparison

Selected official SEA-RAFT-S, commit `9137517ba24e628442aec097d3afe71d03503b75`, model `MemorySlices/Tartan-C-T-TSKH-spring540x960-S`. Official weights SHA-256: `{'models/SEA-RAFT-S/model.safetensors': 'd6a75e47f2630ba6c354ce84a322e24c3d9def668a6956f340a054b9a3211908', 'models/SEA-RAFT-S/README.md': '76dfb05738b57606974a4bf5bebfcce82e294f5ad349a3da82fc2221d3e73415'}`. Code/model BSD-3-Clause. WAFT reviewed only; no second AI flow implementation deployed. Natural-image pretrained weights were not trained on these microscope images.

Actual environment 2.11.0+cu128, CUDA 12.8, NVIDIA GeForce RTX 4060 Laptop GPU. Independent copied CUDA dependency files, no modification to old annotation/Windows environments. Model consumes real two-frame inputs. 256×256 synthetic and 256×384 true recorded ROIs, gray repeated across RGB channels, native resolution, 4 iterations; official normalization/padding. One warmup, forward + reverse flow (forward/backward error <1.5 px), original texture rejection, same Huber/log-transform fit as PIV.

24 pairs: 18 exact synthetic endpoints (six cases × 3 windows, seed 7201) + six recorded pairs. All methods below use the exact same 18 synthetic input pairs and truth. Different acceptance subsets remain explicit.

| Method | Pairs | All-six fraction | MAE vx / vy | Accepted false fraction |
| --- | --- | --- | --- | --- |
| CPD | 18 | 0.5556 | 0.1799/0.1119 | 0 |
| PIV_multipass | 18 | 0.7778 | 0.6946/0.3057 | 0 |
| statistics | 18 | 0 | NOT_AVAILABLE/NOT_AVAILABLE | NOT_AVAILABLE |
| SEA_RAFT_S | 18 | 0.9444 | 6.319/4.291 | 0.3529 |

Actual final profiled p50/p95/max forward+reverse pipeline: 37.99/49.29/61.69 ms. Fit is recorded separately, not hidden inside these numbers. Peak PyTorch allocated/reserved 102.56/126.00 MiB; process RSS and device-wide GPU utilization in `resources` (desktop usage is not isolated). GPU context/desktop memory is larger than the tensor allocation. Initial run preserved; exact-input repeat added process/device resource measurements, not retuning. Both finite runs are retained; no large pressure test.

**No accuracy advantage established.** AI accepts 17/18 synthetic pairs but 6/17 accepted estimates exceed a preregistered error tolerance. PIV accepts 14/18 with 0/14 such failures in this subset. AI can follow fixed background or average tiny structures. Its smoother field is not evidence of accuracy. Retain traditional PIV as the main line; no fine-tuning or world-model training this round. Real-image AI disagreement is diagnostic only because absolute six-parameter truth is unavailable.
