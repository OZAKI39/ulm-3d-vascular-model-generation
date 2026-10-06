"""Three native orthogonal B2 pressure sections in the established publication style."""

import argparse
import json
from pathlib import Path
import shutil

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection
import numpy as np
from PIL import Image
import pyvista as pv

from scripts.stage4 import paper_acoustic_volume as shared

ROOT, FIELD, PRIOR = shared.ROOT, shared.FIELD, shared.PRIOR
SIZE_MM = (180, 90)
STEM = "Fig1b_acoustic_field"
CMAP = plt.get_cmap("magma")
PANEL_WIDTH = 0.34
LOWER_EDGE = 0.225
XY_HEIGHT = PANEL_WIDTH * 2 * 12.4 / 13.4
UPPER_EDGE = LOWER_EDGE + XY_HEIGHT
RIGHT_GAP_MM = 13.0
YZ_HEIGHT = XZ_HEIGHT = (XY_HEIGHT - RIGHT_GAP_MM / SIZE_MM[1]) / 2
SPECS = [
    dict(
        name="XY",
        axes=(0, 1),
        normal=2,
        limits=((-3.2, 10.2), (-6.2, 6.2)),
        box=(0.083, LOWER_EDGE, PANEL_WIDTH, XY_HEIGHT),
    ),
    dict(
        name="YZ",
        axes=(1, 2),
        normal=0,
        limits=((-6.2, 6.2), (-1.6, 1.6)),
        box=(0.580, UPPER_EDGE - YZ_HEIGHT, PANEL_WIDTH, YZ_HEIGHT),
    ),
    dict(
        name="XZ",
        axes=(0, 2),
        normal=1,
        limits=((-3.2, 10.2), (-1.6, 1.6)),
        box=(0.580, LOWER_EDGE, PANEL_WIDTH, XZ_HEIGHT),
    ),
]


def load_sections():
    previous = json.loads((PRIOR / "FIGURE_PROVENANCE.json").read_text())
    mesh_path = ROOT / next(p for p in previous["protected_inputs"] if p.endswith(".vtp"))
    config = ROOT / previous["cases"]["B2"]["run"] / "config.yaml"
    protected = {str(p.relative_to(ROOT)): shared.sha(p) for p in [FIELD, mesh_path, config]}
    for path, digest in protected.items():
        assert previous["protected_inputs"][path] == digest
    helper = Path(shared.__file__)
    protected[str(helper.relative_to(ROOT))] = shared.sha(helper)
    with h5py.File(FIELD) as h:
        axes = [h["axes"][name + "_m"][:] * 1000 for name in "xyz"]
        complex_field = h["pressure_Pa"][:]
        frequency = float(h.attrs["frequency_Hz"])
    assert complex_field.shape == tuple(map(len, axes))
    for axis in axes:
        np.testing.assert_allclose(np.diff(axis), 0.1, atol=1e-12)
    mesh = pv.read(mesh_path)
    sections = []
    for spec in SPECS:
        normal = spec["normal"]
        index = int(np.argmin(np.abs(axes[normal])))
        assert abs(axes[normal][index]) < 1e-12
        pressure = np.take(complex_field, index, axis=normal)
        amplitude = abs(pressure) / 1000
        assert np.isfinite(amplitude).all() and amplitude.max() <= shared.CLIM[1]
        plane = mesh.slice(normal=np.eye(3)[normal], origin=(0, 0, 0))
        np.testing.assert_allclose(plane.points[:, normal], 0, atol=1e-12)
        outline, offset = [], 0
        while offset < len(plane.lines):
            count = plane.lines[offset]
            points = plane.points[plane.lines[offset + 1 : offset + count + 1]] * 1000
            outline.extend(np.stack([points[:-1], points[1:]], axis=1))
            offset += count + 1
        outline = np.asarray(outline)
        assert outline.shape[1:] == (2, 3)
        sections.append(
            dict(
                spec=spec,
                pressure=pressure,
                amplitude=amplitude,
                axes=[axes[i] for i in spec["axes"]],
                outline=outline,
                slice_index=index,
                slice_coordinate_mm=float(axes[normal][index]),
            )
        )
    return sections, protected, frequency


def draw(sections, frequency, output, dpi):
    fig = plt.figure(figsize=np.array(SIZE_MM) / 25.4, dpi=dpi, facecolor="white")
    fig.text(0.045, 0.96, "Acoustic pressure field", fontsize=12, va="top")
    fig.text(0.96, 0.96, f"B2  |  {frequency / 1e6:g} MHz", fontsize=9, va="top", ha="right", color="#425665")
    audit = []
    for section in sections:
        spec = section["spec"]
        ax = fig.add_axes(spec["box"])
        horizontal, vertical = section["axes"]
        image = ax.imshow(
            section["amplitude"].T,
            origin="lower",
            interpolation="nearest",
            cmap=CMAP,
            vmin=shared.CLIM[0],
            vmax=shared.CLIM[1],
            extent=[horizontal[0] - 0.05, horizontal[-1] + 0.05, vertical[0] - 0.05, vertical[-1] + 0.05],
            aspect="equal" if spec["name"] == "XY" else "auto",
        )
        np.testing.assert_array_equal(image.get_array(), section["amplitude"].T)
        ax.add_collection(
            LineCollection(section["outline"][..., list(spec["axes"])], colors="white", linewidths=0.75)
        )
        labels = ["XYZ"[i] for i in spec["axes"]]
        ax.set(
            xlim=spec["limits"][0],
            ylim=spec["limits"][1],
            xlabel=f"{labels[0]} (mm)",
            ylabel=f"{labels[1]} (mm)",
        )
        normal_name = "XYZ"[spec["normal"]]
        ax.set_title(f"{spec['name']} plane ({normal_name} = 0 mm)", fontsize=9, pad=5, loc="left")
        ax.xaxis.labelpad = 2
        ax.set_xticks([-5, 0, 5] if labels[0] == "Y" else [0, 5, 10])
        ax.set_yticks([-5, 0, 5] if labels[1] == "Y" else [-1, 0, 1])
        ax.tick_params(labelsize=8, pad=2)
        audit.append(
            dict(
                name=spec["name"],
                fixed_axis=normal_name,
                fixed_coordinate_mm=section["slice_coordinate_mm"],
                slice_index=section["slice_index"],
                axes=labels,
                native_shape=list(section["amplitude"].shape),
                displayed_limits_mm=spec["limits"],
                pressure_range_kPa=[float(section["amplitude"].min()), float(section["amplitude"].max())],
            )
        )
    cax = fig.add_axes([0.345, 0.12, 0.34, 0.023])
    cb = fig.colorbar(
        plt.cm.ScalarMappable(norm=colors.Normalize(*shared.CLIM), cmap=CMAP),
        cax=cax,
        orientation="horizontal",
    )
    cb.set_ticks([0, 10, 20, 30, 39])
    cb.set_label(r"Pressure amplitude, $|P|$ (kPa)", fontsize=8, labelpad=4)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(labelsize=8, width=0.6, length=3)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            bb = artist.get_window_extent(renderer)
            assert bb.x0 >= 0 and bb.y0 >= 0 and bb.x1 <= width + 1 and bb.y1 <= height + 1, artist.get_text()
    # The user requested taller right panels. Only the display aspect changes;
    # coordinates, field samples and the vessel intersections remain unchanged.
    scales = []
    for spec, ax in zip(SPECS, fig.axes[:3]):
        coordinates = ax.transData.transform([[0, 0], [1, 0], [0, 1]])
        x_scale = np.linalg.norm(coordinates[1] - coordinates[0])
        y_scale = np.linalg.norm(coordinates[2] - coordinates[0])
        if spec["name"] == "XY":
            np.testing.assert_allclose(x_scale, y_scale, rtol=1e-10)
        scales.append(
            dict(
                plane=spec["name"],
                horizontal=float(x_scale),
                vertical=float(y_scale),
                vertical_display_scale=float(y_scale / x_scale),
            )
        )
    xy, yz, xz = [ax.get_window_extent(renderer) for ax in fig.axes[:3]]
    alignment = dict(
        top=yz.y1 - xy.y1,
        bottom=xz.y0 - xy.y0,
        right_column_left=yz.x0 - xz.x0,
        right_column_right=yz.x1 - xz.x1,
    )
    assert max(abs(value) for value in alignment.values()) < 1e-8, alignment
    gap_mm = (yz.y0 - xz.y1) / dpi * 25.4
    np.testing.assert_allclose(gap_mm, RIGHT_GAP_MM, atol=1e-8)
    upper_xlabel = fig.axes[1].xaxis.label.get_window_extent(renderer)
    lower_title = fig.axes[2]._left_title.get_window_extent(renderer)
    assert upper_xlabel.y0 > lower_title.y1, ("Right-panel labels overlap", upper_xlabel.y0, lower_title.y1)
    if dpi == 600:
        fig.savefig(output / f"{STEM}.png", dpi=600)
        fig.savefig(
            output / f"{STEM}.pdf", dpi=600, metadata={"Title": "B2 XY, YZ and XZ acoustic pressure sections"}
        )
        fig.savefig(output / f"{STEM}.svg", dpi=600)
        with Image.open(output / f"{STEM}.png") as im:
            im.convert("RGB").save(output / f"{STEM}.tif", compression="tiff_lzw", dpi=(600, 600))
            im.convert("RGB").resize(tuple(round(v * 0.3) for v in im.size), Image.Resampling.LANCZOS).save(
                output / f"{STEM}_preview.png", dpi=(180, 180)
            )
    else:
        fig.savefig(output / f"{STEM}_preview.png", dpi=dpi)
    plt.close(fig)
    return dict(
        sections=audit,
        equal_aspect_planes=["XY"],
        pixel_per_mm=scales,
        right_panel_gap_mm=gap_mm,
        right_panel_height_mm=YZ_HEIGHT * SIZE_MM[1],
        display_aspect="YZ/XZ enlarged vertically at user request; native axes and data unchanged",
        column_alignment_residual_pixels=alignment,
        interpolation="nearest native nodes",
        colormap="magma",
        shared_clim_kPa=shared.CLIM,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    shared.style()
    sections, protected, frequency = load_sections()
    audit = draw(sections, frequency, args.output, 180 if args.preview else 600)
    assert all(shared.sha(ROOT / p) == digest for p, digest in protected.items())
    report = dict(
        status="PREVIEW" if args.preview else "PASS",
        case="B2",
        frequency_Hz=frequency,
        protected_inputs=protected,
        renderer_sha256=shared.sha(__file__),
        figure_size_mm=SIZE_MM,
        rendering=audit,
        original_inputs_unchanged=True,
        new_simulation=False,
        field_definition="Magnitude of saved complex phasor in kPa, not RMS",
        style="White background, Arial regular, original magma purple-to-yellow map, dark axes, white vessel sections",
        export="600 dpi raster pressure maps; editable vector text and vessel outlines in PDF/SVG",
    )
    (args.output / "PROVENANCE.json").write_text(json.dumps(report, indent=2) + "\n")
    if not args.preview:
        shutil.copy2(__file__, args.output / Path(__file__).name)
        shutil.copy2(shared.__file__, args.output / "paper_acoustic_volume.py")
        with h5py.File(args.output / "source_data.h5", "w") as h:
            h.attrs["source_sha256"] = shared.sha(FIELD)
            h.attrs["frequency_Hz"] = frequency
            for section in sections:
                spec = section["spec"]
                group = h.create_group(spec["name"])
                group.attrs["fixed_axis"] = "XYZ"[spec["normal"]]
                group.attrs["fixed_coordinate_mm"] = section["slice_coordinate_mm"]
                group.attrs["fixed_index"] = section["slice_index"]
                for name, axis in zip(spec["name"], section["axes"]):
                    group.create_dataset(f"{name}_mm", data=axis)
                group.create_dataset("pressure_complex_Pa", data=section["pressure"], compression="gzip")
                group.create_dataset("pressure_amplitude_kPa", data=section["amplitude"], compression="gzip")
                group.create_dataset("vessel_outline_xyz_mm", data=section["outline"], compression="gzip")
        (args.output / "CAPTION.txt").write_text(
            "Orthogonal sections of the B2 acoustic pressure amplitude at 2.5 MHz. The XY, YZ and XZ planes "
            "are extracted at Z=0, X=0 and Y=0 mm, respectively, from the existing complex acoustic field. "
            "Every section uses native 0.1 mm grid nodes without field smoothing or interpolation; the "
            "raster display uses nearest-neighbor sampling. The shared linear color scale is 0-39 kPa and "
            "represents |P|, the magnitude of the complex pressure phasor, not RMS pressure. The original magma "
            "purple-to-yellow color map is restored. White curves "
            "show exact intersections of the existing vessel mesh with each plane. XY preserves equal "
            "spatial aspect. YZ and XZ are enlarged vertically for readability, with explicit mm coordinates; "
            "their unequal horizontal/vertical display scales are recorded in PROVENANCE.json. The two "
            "right panels have identical widths and heights, with a 13 mm gap. Their combined upper/lower "
            "edges align with the XY panel. This display aspect does not change the source coordinates, "
            "pressure samples or vessel geometry. Display ranges are X=[-3.2,10.2], "
            "Y=[-6.2,6.2] and Z=[-1.6,1.6] mm; X and Y limits match the preceding right-hand plane. "
            "The homogeneous-water acoustic model does not resolve vessel-wall material interfaces. "
            "No simulation was run and no original source data were modified. The figure contains no "
            "volume rendering, isosurfaces or volume ROI locator. PDF/SVG use raster field maps and editable "
            "vector labels and vessel outlines.\n"
        )
    print(json.dumps(dict(status=report["status"], output=str(args.output), rendering=audit)), flush=True)


if __name__ == "__main__":
    main()
