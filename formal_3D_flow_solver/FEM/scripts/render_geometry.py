#!/usr/bin/env python3
"""Render the existing source triangles and metadata-defined cap normals only."""
import json
import sys
from pathlib import Path
import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from fem3d.audit import sha256, write_json
from fem3d.source import validate_contract

c = json.loads((ROOT/"reports/stage00/source_contract.json").read_text())
validate_contract(c)
assert sha256(c["geometry_path"]) == c["geometry_sha256"]
surface = pv.read(c["geometry_path"])
ports = c["inlets"] + c["outlets"]
plot = pv.Plotter(off_screen=True, shape=(1,2), window_size=(2200,1200), border=False)
for index in range(2):
    plot.subplot(0,index)
    plot.set_background("#f7f9fc")
    wall = surface.extract_cells(np.flatnonzero(surface["CellEntityIds"] == c["wall_entity_id"])).extract_surface()
    plot.add_mesh(wall,color="#93a6b5",smooth_shading=True,specular=0.16)
    label_points = []
    names = []
    for port in ports:
        color = "#df6c2e" if port["role"] == "inlet" else "#087f8c"
        cap = surface.extract_cells(np.flatnonzero(surface["CellEntityIds"] == port["surface_entity_id"])).extract_surface()
        plot.add_mesh(cap,color=color,lighting=False)
        center = np.array(port["plane_origin_m"])*1e6
        normal = np.array(port["outward_normal"])
        plot.add_mesh(pv.Arrow(start=center,direction=normal,scale=9,tip_length=.26,tip_radius=.11,shaft_radius=.025),color=color)
        label_points.append(center + normal*13)
        names.append(port["name"])
    plot.add_point_labels(np.array(label_points),names,font_size=22,text_color="#142535",font_family="arial",
                          shape_color="#ffffff",shape_opacity=.93,show_points=False,always_visible=True,margin=7)
    plot.add_text("Tagged source surface" if index==0 else "Opening locations | XY view",position="upper_left",font_size=19,color="#142535")
    plot.add_text("Orange: inlet   Teal: outlets   Arrows: outward normals\nDisplay coordinates: micrometres (um)",
                  position="lower_left",font_size=13,color="#334a59")
    plot.add_axes(line_width=3)
    plot.enable_parallel_projection()
    if index==0:
        plot.view_isometric()
    else:
        plot.view_xy()
    plot.reset_camera()
    plot.camera.zoom(.85)
plot.screenshot(ROOT/"reports/stage00/source_geometry_overview.png")
plot.close()
write_json(ROOT/"reports/stage00/geometry_figure_metadata.json",{
    "geometry_path":c["geometry_path"],"geometry_sha256":c["geometry_sha256"],
    "labels":[p["name"] for p in ports],"triangle_count":surface.n_cells,
    "render_only":True,"volume_mesh_generated":False,"flow_solved":False,
    "image_sha256":sha256(ROOT/"reports/stage00/source_geometry_overview.png")})
print("Saved reports/stage00/source_geometry_overview.png")
