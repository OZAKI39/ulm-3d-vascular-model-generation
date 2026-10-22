# Source reuse and historical protection

Public baseline: `9cb22d729fece92fe160082edf09222a18068d74`. Original independent vision repository: `0561769933dfa799dee0778f638a8c26270c3b95`. New isolated worktree contains the synchronized byte-identical F10U modules. Server physics worktree was independently audited and is an older F5 checkout, with later physics snapshots stored separately; it is not silently used as the F10V source.

`src/legacy.py` is a read-only adapter to F10U `affine.py`, `distribution.py`, `methods.py`, `state.py`, `synthetic.py`. It directly reuses CandidateDetector (original OpenCV threshold 10), distribution, single/multipass OpenPIV, pycpd, affine parameter definitions and matrix log, CausalWindow, old ConsumerGate, MockConsumer and state construction. It imports the original camera module but never constructs a camera receiver or executes device commands. Camera and AFG sources are protected.

F10V adds only project-specific triple/split extensions, heatmap fusion, official-backbone heads, differentiable geometric adaptation, scientific evaluation and a versioned M4 adapter. OpenCV warping/connected components/Gaussian filtering, SciPy expm/logm, PyTorch linalg/autograd/AMP/DataLoader, OpenPIV, pycpd, PyAV encoding/decoding, jsonschema and Matplotlib are called directly. No optical flow, CNN backbone, video codec, generic linear solver or camera protocol is reimplemented.

Protected-source inventory is `evidence/protected_inputs_before.json`. Final verification reports any change separately. Existing F10U accuracy, failures, synthetic seeds, archived results and hardware-control status remain unchanged. The original physical simulation is not rerun or relabeled.

Official references used in the implementation review:

- [TorchVision MobileNetV3-small](https://docs.pytorch.org/vision/stable/models/generated/torchvision.models.mobilenet_v3_small.html)
- [PyTorch AMP](https://docs.pytorch.org/docs/2.14/amp.html)
- [OpenPIV](https://github.com/OpenPIV/openpiv-python), [pycpd](https://github.com/siavashk/pycpd)
- [Conformal prediction introduction](https://arxiv.org/abs/2107.07511): finite-sample rank construction requires exchangeability; no cross-domain real guarantee is claimed.
- [SynthPix](https://github.com/antonioterpin/synthpix): official principal example produces image pairs; continuous three-frame identity was not assumed.

Optional NeuFlow v2 and SynthPix experiments are explicitly NOT_RUN, with reasons and future entry points recorded; no author benchmark is substituted for a microbubble test.
