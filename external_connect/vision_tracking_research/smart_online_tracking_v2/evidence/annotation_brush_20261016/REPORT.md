# Default brush annotation and new-bubble shortcut

The workbench now opens with SINGLE_BUBBLE selected, an unused instance ID,
and the native napari paint brush active. Frame navigation also enters paint
mode for the current class. N selects SINGLE_BUBBLE, prepares the next unused
ID, returns to the main tab and activates the brush. Multiple strokes before
N remain one instance. P/E select brush/eraser without changing IDs.

The main tab explains the brush workflow and exposes an explicit N button.
Other classes have a separate new-object brush button. Polygon drawing is
still available; an unfinished polygon blocks N/P/E without losing vertices.
Text input keeps normal letter entry. Brush size is preserved between objects.
The revised ANNOTATION_READY.md describes painting filled visible regions,
separating instances, shortcuts, correction, and the unchanged review process.

Validation: 26 actual Qt/napari interaction checks passed, including real mouse
brush marks, N/P/E/bracket keys, independent IDs, protected existing pixels,
undo/redo, text entry, SAM-to-manual switching, polygon compatibility,
frame navigation, exact save/reload and restored-layer bindings. Tests used
an independent evidence sandbox and never launched automatic training.
All 31 pre-existing annotation files retained their SHA-256 hashes.

The initial test caught napari resetting paint to pan/zoom when visibility was
assigned even on an already visible layer. Conditional visibility changes
before mode activation fixed this. The focus test also needed to activate the
test window before asserting focus under WSL. Initial failures and diagnostic
traces are retained here; ui_checks.json records the final successful run.

No formal labels, SAM weights, training logic or scientific results changed.
The screenshot contains software fixtures only, never human ground truth.
See manifest.json for baseline identity and modified-source SHA-256 hashes.
