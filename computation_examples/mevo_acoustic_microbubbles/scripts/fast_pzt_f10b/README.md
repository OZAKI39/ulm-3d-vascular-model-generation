# F10B reproducibility

Use the project Python with `PYTHONPATH=src:.`. Local plots use the shared WSL environment; science uses the existing server `.venv` (Warp1.18.0+cu12). Always commit executable source before `remote.py`; each label is unique and failed attempts are retained.

Executed science commands, in sequence:

```bash
python -m scripts.fast_pzt_f10b.remote contact --label contact_v1
python -m scripts.fast_pzt_f10b.remote contact --label contact_v2 --cases A2 A3
python -m scripts.fast_pzt_f10b.remote checks --label checks_v1
python -m scripts.fast_pzt_f10b.remote rollout --label rollout_v1
python -m scripts.fast_pzt_f10b.remote support --label support_v1
python -m scripts.fast_pzt_f10b.remote checks --label checks_v2
```

The first contact job intentionally remains an infrastructure-failure record after saving A0/A1: oldF3 refused24-round construction; only A2/A3 were subsequently run. Each commit is in the corresponding IDENTITY.json. Do not rerun these names or replace historical output.

Offline / read-only archive steps:

```bash
python -m scripts.fast_pzt_f10b.archive
python -m scripts.fast_pzt_f10b.analyze
python -m scripts.fast_pzt_f10b.figures
python -m scripts.fast_pzt_f10b.verify
```

Figure review and progress-document snapshot must exist before final verification; the verifier checks their actual hashes. Figures use saved data only. The scientific plots need NumPy, h5py and Matplotlib; PyMuPDF is used for PDF rendering/review. No new package was installed.

Current G1 numerical admission is limited to this initial template / OFF,S,O /10ms. The same-code full reference timer reproduces F10A10us particle states bitwise. The measured standalone speed gate fails; do not promote the integration-only6.3x number into a whole-program or production-batch claim. No production re-entry / morphology interface is provided.
