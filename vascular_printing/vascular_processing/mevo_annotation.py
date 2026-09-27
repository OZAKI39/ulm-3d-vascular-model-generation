"""Small optional landmark picker, isolated from the existing ROI visualization."""
from datetime import datetime, timezone
import os
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree
import yaml

from .mevo_export import source_polydata
from .mevo_graph import load_exact_graph
from .roi_landmarks import ROI_DEFINITIONS, landmark_template, load_landmarks, write_landmarks


def annotate(source: Path, output: Path, *, roi_name="RMCA_M2M3", units="mm"):
    source, output = source.resolve(), output.resolve()
    if output == source:
        raise ValueError("Landmarks cannot overwrite source SWC")
    if output.exists():
        document, _ = load_landmarks(output, source)
    else:
        document = landmark_template(source, output, units=units)
        write_landmarks(output, document)
    if sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        return {"status": "NEEDS_MANUAL_REVIEW", "landmarks": str(output),
                "message": "No DISPLAY: template saved. Use inspect + ParaView to read original_swc_id, then edit YAML; no anatomical boundaries inferred."}
    import pyvista as pv
    exact = load_exact_graph(source)
    ids = list(exact.graph)
    points = np.asarray([exact.graph.nodes[n]["coords"][:3] for n in ids])
    locator = cKDTree(points)
    names = list(ROI_DEFINITIONS)
    state = {"roi": roi_name, "role": "proximal_nodes"}
    plotter = pv.Plotter(window_size=(1200, 850))
    plotter.add_mesh(source_polydata(exact), color="cyan", line_width=2, pickable=True)
    plotter.add_axes()

    def refresh():
        roi = document["rois"][state["roi"]]
        for key, color in (("proximal_nodes", "lime"), ("distal_nodes", "red"), ("exclude_subtree_roots", "orange")):
            plotter.remove_actor(key, render=False)
            if roi[key]:
                selected = np.asarray([exact.graph.nodes[n]["coords"][:3] for n in roi[key]])
                plotter.add_point_labels(selected, [str(n) for n in roi[key]], point_color=color,
                                         point_size=12, text_color=color, name=key, always_visible=True)
        text = ("ANATOMICAL LANDMARK ANNOTATION (separate inspection tool)\n"
                + " / ".join(f"{i+1}: {name}" for i, name in enumerate(names))
                + f"\nROI: {state['roi']} | pick role: {state['role']}\n"
                + "P proximal; D distal; E exclude; right-click centerline toggles nearest ORIGINAL node ID\n"
                + "B confirmed proximal boundary; T truncated; V toggle manual verification; W save YAML\n"
                + f"boundary={roi['proximal_boundary_status']} | manually_verified={roi['manually_verified']}\n"
                + f"proximal={roi['proximal_nodes']} | distal={roi['distal_nodes']} | excluded={roi['exclude_subtree_roots']}")
        plotter.add_text(text, font_size=10, name="annotation_status")
        plotter.render()

    def pick(point):
        node = ids[int(locator.query(point)[1])]
        roi = document["rois"][state["roi"]]
        values = roi[state["role"]]
        if node in values:
            values.remove(node)
        else:
            values.append(node)
        roi["enabled"] = True
        # Editing boundaries invalidates previous human confirmation; never promote to true.
        roi["manually_verified"] = False
        refresh()

    def choose_roi(name):
        state["roi"] = name
        refresh()

    def choose_role(role):
        state["role"] = role
        refresh()

    def set_boundary(status):
        roi = document["rois"][state["roi"]]
        roi["proximal_boundary_status"] = status
        roi["manually_verified"] = False
        refresh()

    def verify():
        roi = document["rois"][state["roi"]]
        if not roi["manually_verified"]:
            from .mevo_roi import extract_graph
            try:
                extract_graph(exact, roi)
            except ValueError as exc:
                print(f"Cannot verify incomplete boundaries: {exc}", flush=True)
                return
        roi["manually_verified"] = not roi["manually_verified"]
        refresh()

    def save():
        document["annotation_date"] = datetime.now(timezone.utc).isoformat()
        backup = output.with_name(output.name + "." + datetime.now().strftime("%Y%m%d_%H%M%S_%f") + ".bak")
        with backup.open("xb") as stream:
            stream.write(output.read_bytes())
        temporary = output.with_name(output.name + ".pending")
        write_landmarks(temporary, document)
        temporary.replace(output)
        print(f"Saved anatomical annotation: {output}; backup: {backup}", flush=True)

    for i, name in enumerate(names, 1):
        plotter.add_key_event(str(i), lambda name=name: choose_roi(name))
    for key, role in (("p", "proximal_nodes"), ("d", "distal_nodes"), ("e", "exclude_subtree_roots")):
        plotter.add_key_event(key, lambda role=role: choose_role(role))
    plotter.add_key_event("b", lambda: set_boundary("confirmed"))
    plotter.add_key_event("t", lambda: set_boundary("truncated_to_available_data"))
    plotter.add_key_event("v", verify)
    plotter.add_key_event("w", save)
    plotter.enable_point_picking(callback=pick, picker="cell", show_message=False, show_point=False,
                                left_clicking=False, pickable_window=False)
    refresh()
    plotter.show()
    return {"status": "PASS_WITH_WARNINGS", "landmarks": str(output),
            "message": "Picker closed; only explicit W saves were written. Run extract to validate coverage."}
