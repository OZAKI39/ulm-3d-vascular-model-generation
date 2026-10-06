# One editable development-workflow slide

Deliverable: `development_workflow_current_progress.pptx` — exactly one 16:9 slide. The PNG and PDF are rendered from this PowerPoint. Workflow text, boxes, arrows, progress indicators and the clot schematic are native editable objects; only the three real-data insets are images.

From the `clot_clearance/` project directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/clot_slide_mpl PYTHONDONTWRITEBYTECODE=1 \
/home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B \
analysis/development_workflow_slide/generate_slide.py
```

Edit `slide_config.py` to steer the narrative and layout content. `make_insets.py` reads existing CSV/JSON results only. The generator uses python-pptx, NumPy and Matplotlib from the existing shared environment. Rendering uses LibreOffice on PATH, or the task-local extracted runtime at `/tmp/clot_slide_runtime_20261005/root`; the existing project PyMuPDF reader generates the 300 dpi preview. No system packages or scientific libraries were replaced.

Sources, the numerical statements, interpretation limits and steering instructions are recorded in `development_workflow_current_progress_notes.md` and in the editable PowerPoint speaker notes. `source_manifest.json` records only files used for this slide. `validate_slide.py` checks the one-slide structure, editable object counts, source provenance, English text and exported page bounds.

Visual review of the actual PowerPoint render: no text clipping or overlap; stage 4 highlighted; stage 5 marked NEXT; stage 6 marked FUTURE; native clot schematic attached to the wall; no default shadows. Root cause remains INCONCLUSIVE. The next 3δ experiment has not been run.
