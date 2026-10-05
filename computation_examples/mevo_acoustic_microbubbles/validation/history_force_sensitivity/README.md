# Independent history-force sensitivity

Reference-only code, excluded from the installed production package. No tracker import is permitted. See `docs/HISTORY_FORCE_SENSITIVITY.md` for the equations, three SonoVue conditions, numerical convergence, separate averaging/added-mass/history comparisons and limitations.

Run from the repository root with the project environment:

```bash
python -m scripts.stage3.history
python -m scripts.stage3.history_report
```

On the current remote instance use the finite supervisor runner documented in `docs/STAGE3_QUICKSTART.md`. Full convolution is retained for this small study; this is not a scalable population solver.

Sources read: Garbin et al. 2009, complete author PDF via the web reader (local direct PDF retrieval returned HTTP 403); Cattaneo & Supponen, full arXiv v2 PDF/HTML, 7 February 2026, with final article identity confirmed at APS. DOI, equation numbering and version distinctions are in the report. No restricted dataset was accessed.
