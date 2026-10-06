"""Render the existing B2 pressure volume without changing the acoustic data."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image
import pyvista as pv
from scipy.interpolate import RegularGridInterpolator
from vtkmodules.vtkCommonDataModel import vtkPiecewiseFunction

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "results/stage4/publication_figures/workflow_B2_B3_20261005_v1"
REFERENCE = ROOT / "results/stage4/publication_figures/PZT_control_reference_v6"
FIELD = ROOT / "results/stage4/phase_screening/20261002T144847677887Z/B2_field.h5"
INK = "#263746"
CMAP = colors.LinearSegmentedColormap.from_list(
    "pressure_blue", ["#f7fbfd", "#cedde7", "#9bbace", "#678ea9", "#365f7e", "#15354e"], N=256
)
CLIM = (0.0, 39.0)
ROI = np.array([[-3, 3], [-2, 2], [-1.6, 1.6]])
# These are display opacity values per 0.2 mm, not physical attenuation.
OPACITY = np.array([[0, 0], [5, 0], [12, 0.008], [20, 0.035], [27, 0.15], [39, 0.55]])
ISOSURFACES = [(28.0, 0.38), (35.0, 0.82)]


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def style():
    fonts = Path("/mnt/c/Windows/Fonts")
    for path in fonts.glob("arial*.ttf"):
        font_manager.fontManager.addfont(str(path))
    plt.rcParams.update(
        {
            "font.family": "Arial" if (fonts / "arial.ttf").is_file() else "DejaVu Sans",
            "font.size": 8,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": "#647887",
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "mathtext.fontset": "stixsans",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def load_data():
    prior = json.loads((PRIOR / "FIGURE_PROVENANCE.json").read_text())
    config = ROOT / prior["cases"]["B2"]["run"] / "config.yaml"
    surface = ROOT / next(p for p in prior["protected_inputs"] if p.endswith(".vtp"))
    protected = {str(p.relative_to(ROOT)): sha(p) for p in [FIELD, config, surface]}
    for path, digest in protected.items():
        assert prior["protected_inputs"][path] == digest, path
    for path in [REFERENCE / "Fig1a_PZT_control.png", REFERENCE / "PROVENANCE.json"]:
        protected[str(path.relative_to(ROOT))] = sha(path)
    with h5py.File(FIELD) as h:
        axes = [h["axes"][name][:] * 1000 for name in ["x_m", "y_m", "z_m"]]
        complex_pressure = h["pressure_Pa"][:]
        frequency = float(h.attrs["frequency_Hz"])
    amplitude = np.abs(complex_pressure) / 1000
    assert amplitude.shape == tuple(len(a) for a in axes)
    assert np.isfinite(amplitude).all() and amplitude.max() <= CLIM[1]
    for a in axes:
        np.testing.assert_allclose(np.diff(a), 0.1, atol=1e-12)
    selected = [np.flatnonzero((a >= lo - 1e-10) & (a <= hi + 1e-10)) for a, (lo, hi) in zip(axes, ROI)]
    roi_axes = [a[i] for a, i in zip(axes, selected)]
    roi = amplitude[np.ix_(*selected)]
    assert roi.shape == (61, 41, 33)
    grid = pv.ImageData(dimensions=roi.shape, spacing=(0.1, 0.1, 0.1), origin=[a[0] for a in roi_axes])
    grid["pressure_amplitude_kPa"] = roi.ravel(order="F")
    native_points = np.stack(np.meshgrid(*roi_axes, indexing="ij"), axis=-1).reshape(-1, 3, order="F")
    np.testing.assert_allclose(grid.points, native_points, atol=1e-12)
    np.testing.assert_array_equal(grid["pressure_amplitude_kPa"].reshape(roi.shape, order="F"), roi)
    # The right-panel vessel outline comes from the existing mesh, not a schematic.
    section = pv.read(surface).slice(normal="z", origin=(0, 0, 0))
    lines, offset, outline = section.lines, 0, []
    while offset < len(lines):
        n = lines[offset]
        points = section.points[lines[offset + 1 : offset + n + 1]] * 1000
        outline.extend(np.stack([points[:-1], points[1:]], axis=1))
        offset += n + 1
    iz = int(np.argmin(abs(axes[2])))
    assert abs(axes[2][iz]) < 1e-12
    return dict(
        axes=axes,
        roi_axes=roi_axes,
        roi=roi,
        grid=grid,
        pressure=amplitude,
        complex_roi=complex_pressure[np.ix_(*selected)],
        section=amplitude[:, :, iz],
        section_complex=complex_pressure[:, :, iz],
        outline=np.asarray(outline),
        protected=protected,
        frequency_Hz=frequency,
    )


def add_lines(plotter, segments, color, width):
    points = np.asarray(segments).reshape(-1, 3)
    lines = np.array([[2, 2 * i, 2 * i + 1] for i in range(len(segments))])
    return plotter.add_mesh(pv.PolyData(points, lines=lines), color=color, line_width=width, lighting=False)


def render(data, dpi, box, figsize):
    size = np.rint(figsize * dpi * box[2:]).astype(int)
    pl = pv.Plotter(off_screen=True, window_size=size.tolist(), lighting="three lights")
    pl.set_background("white")
    pl.render_window.SetMultiSamples(0)
    pl.render_window.SetAlphaBitPlanes(1)
    pl.enable_depth_peeling(number_of_peels=100, occlusion_ratio=0.0)
    pl.renderer.SetUseDepthPeelingForVolumes(True)
    volume = pl.add_volume(
        data["grid"],
        scalars="pressure_amplitude_kPa",
        clim=CLIM,
        cmap=CMAP,
        opacity="linear",
        opacity_unit_distance=0.2,
        shade=False,
        show_scalar_bar=False,
        mapper="gpu",
        blending="composite",
    )
    # Set the native kPa transfer function explicitly, avoiding uint8 opacity quantization.
    transfer = vtkPiecewiseFunction()
    for pressure, opacity in OPACITY:
        transfer.AddPoint(float(pressure), float(opacity))
    volume.prop.SetScalarOpacity(transfer)
    volume.prop.interpolation_type = "linear"
    volume.mapper.SetAutoAdjustSampleDistances(False)
    volume.mapper.SetSampleDistance(0.025)
    volume.mapper.SetImageSampleDistance(1.0)
    np.testing.assert_array_equal(
        volume.mapper.dataset["pressure_amplitude_kPa"], data["roi"].ravel(order="F")
    )
    for pressure, opacity in OPACITY:
        assert abs(volume.prop.GetScalarOpacity().GetValue(float(pressure)) - opacity) < 1e-12
    # Two labelled, un-smoothed isosurfaces clarify depth within the continuous volume.
    surfaces, actors = [], []
    sample = RegularGridInterpolator(data["roi_axes"], data["roi"])
    for level, opacity in ISOSURFACES:
        surface = data["grid"].contour([level], scalars="pressure_amplitude_kPa")
        # VTK contour coordinates are float32; allow only their boundary rounding
        # when checking against the authoritative float64 node coordinates.
        checked_points = np.clip(surface.points, ROI[:, 0], ROI[:, 1])
        boundary_rounding = float(np.max(np.abs(checked_points - surface.points)))
        assert boundary_rounding < 1e-6
        errors = np.abs(sample(checked_points) - level)
        assert errors.max() < 1e-4
        surface = surface.compute_normals(point_normals=True, cell_normals=False, inplace=False)
        actors.append(
            pl.add_mesh(
                surface,
                color=CMAP(level / CLIM[1])[:3],
                opacity=opacity,
                smooth_shading=True,
                ambient=0.55,
                diffuse=0.45,
                specular=0.08,
                show_scalar_bar=False,
            )
        )
        surfaces.append(
            dict(
                level_kPa=level,
                opacity=opacity,
                cells=surface.n_cells,
                max_native_interpolation_error_kPa=float(errors.max()),
                maximum_boundary_rounding_mm=boundary_rounding,
            )
        )
    xmin, xmax, ymin, ymax, zmin, zmax = -3.3, 3.3, -2.3, 2.3, -1.8, 1.8
    segments = []
    for x in [xmin, -2, 0, 2, xmax]:
        if x != xmax:
            segments.extend([[[x, ymin, zmin], [x, ymax, zmin]], [[x, ymax, zmin], [x, ymax, zmax]]])
    for y in [ymin, -1, 0, 1, ymax]:
        if y != ymin:
            segments.append([[xmin, y, zmin], [xmax, y, zmin]])
        segments.append([[xmin, y, zmin], [xmin, y, zmax]])
    for z in [zmin, -1, 0, 1, zmax]:
        segments.extend([[[xmin, ymax, z], [xmax, ymax, z]], [[xmin, ymin, z], [xmin, ymax, z]]])
    actors.append(add_lines(pl, segments, "#bec9d1", max(0.6, dpi / 180 * 0.5)))
    edges = [
        [[xmin, ymin, zmin], [xmax, ymin, zmin]],
        [[xmax, ymin, zmin], [xmax, ymax, zmin]],
        [[xmax, ymax, zmin], [xmax, ymax, zmax]],
    ]
    actors.append(add_lines(pl, edges, "#647887", max(0.8, dpi / 180 * 0.7)))
    az, el = np.deg2rad([-48, 26])
    direction = -np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    right = np.cross(direction, [0, 0, 1])
    right /= np.linalg.norm(right)
    up = np.cross(right, direction)
    corners = np.array(list(itertools.product([xmin, xmax], [ymin, ymax], [zmin, zmax])))
    horizontal, vertical = corners @ right, corners @ up
    focal = (horizontal.min() + horizontal.max()) / 2 * right + (vertical.min() + vertical.max()) / 2 * up
    scale = max(np.ptp(vertical) / (2 * 0.92), np.ptp(horizontal) / (2 * 0.92 * size[0] / size[1]))
    eye = focal - direction * 30
    pl.camera_position = [eye.tolist(), focal.tolist(), up.tolist()]
    pl.enable_parallel_projection()
    pl.camera.parallel_scale = float(scale)
    pl.show(auto_close=False)
    image = pl.screenshot(return_img=True)
    assert pl.renderer.GetLastRenderingUsedDepthPeeling()
    for actor in actors:
        pl.renderer.RemoveActor(actor)
    for actor in reversed(actors):
        pl.renderer.AddActor(actor)
    pl.render()
    delta = np.abs(image.astype(np.int16) - pl.screenshot(return_img=True).astype(np.int16))
    assert delta.max() <= 2, (delta.max(), delta.mean())
    # A renderer must show the volume itself, not just the superposed surfaces and grid.
    volume.SetVisibility(False)
    pl.render()
    no_volume = pl.screenshot(return_img=True)
    volume_delta = np.abs(image.astype(np.int16) - no_volume.astype(np.int16))
    assert np.count_nonzero(volume_delta.max(axis=-1) > 3) > 0.02 * size.prod()
    pl.close()
    return image, dict(
        box=box.tolist(),
        figsize=figsize.tolist(),
        size=size.tolist(),
        camera_eye_mm=eye.tolist(),
        camera_focal_mm=focal.tolist(),
        camera_right=right.tolist(),
        camera_up=up.tolist(),
        parallel_scale_mm=float(scale),
        camera_azimuth_deg=-48,
        camera_elevation_deg=26,
        grid_bounds_mm=[xmin, xmax, ymin, ymax, zmin, zmax],
        isosurfaces=surfaces,
        opacity_knots_kPa_alpha=OPACITY.tolist(),
        opacity_unit_distance_mm=0.2,
        ray_sample_distance_mm=0.025,
        blend_mode="composite",
        interpolation="trilinear display only",
        shade_volume=False,
        depth_peeling_used=True,
        reverse_actor_max_rgb_delta=int(delta.max()),
        visible_volume_pixel_fraction=float(np.mean(volume_delta.max(axis=-1) > 3)),
    )


def project(point, audit):
    delta = np.asarray(point) - audit["camera_focal_mm"]
    scale, size = audit["parallel_scale_mm"], audit["size"]
    u = 0.5 + delta @ audit["camera_right"] / (2 * scale * size[0] / size[1])
    v = 0.5 + delta @ audit["camera_up"] / (2 * scale)
    x, y, w, h = audit["box"]
    return [x + w * u, y + h * v]


def annotate_volume(fig, ax, audit):
    _, xmax, ymin, ymax, zmin, _ = audit["grid_bounds_mm"]
    entries = [([x, ymin, zmin], str(x), (0, -10), "center") for x in [-2, 0, 2]]
    entries += [([xmax, y, zmin], str(y), (8, -5), "left") for y in [-2, 0, 2]]
    entries += [([xmax, ymax, z], str(z), (7, 0), "left") for z in [-1, 0, 1]]
    for point, text, offset, ha in entries:
        ax.annotate(
            text,
            xy=project(point, audit),
            xycoords=fig.transFigure,
            xytext=offset,
            textcoords="offset points",
            ha=ha,
            va="center",
            fontsize=8,
            color=INK,
            annotation_clip=False,
        )
    for point, text, offset in [
        ([0, ymin, zmin], "X (mm)", (0, -25)),
        ([xmax, 0, zmin], "Y (mm)", (22, -19)),
        ([xmax, ymax, 0], "Z (mm)", (22, 14)),
    ]:
        ax.annotate(
            text,
            xy=project(point, audit),
            xycoords=fig.transFigure,
            xytext=offset,
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=9,
            color=INK,
            annotation_clip=False,
        )


def make_figure(data, output, dpi):
    figsize = np.array([180, 104]) / 25.4
    box = np.array([0.014, 0.19, 0.586, 0.62])
    image, audit = render(data, dpi, box, figsize)
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
    fig.text(0.045, 0.95, "Three-dimensional acoustic field", fontsize=12, va="top")
    fig.text(0.96, 0.95, "B2  |  2.5 MHz", fontsize=9, va="top", ha="right", color="#425665")
    fig.text(0.045, 0.855, "3D pressure volume", fontsize=9)
    ax = fig.add_axes(box)
    ax.imshow(image, interpolation="none")
    ax.axis("off")
    annotate_volume(fig, ax, audit)
    fig.text(0.045, 0.09, "Isosurfaces: 28 and 35 kPa", fontsize=8, color=INK)
    # The quantitative slice remains as context; its dashed ROI maps exactly to the volume.
    ax2 = fig.add_axes([0.68, 0.36, 0.255, 0.39])
    x, y, _ = data["axes"]
    ax2.imshow(
        data["section"].T,
        origin="lower",
        interpolation="nearest",
        cmap=CMAP,
        vmin=CLIM[0],
        vmax=CLIM[1],
        extent=[x[0] - 0.05, x[-1] + 0.05, y[0] - 0.05, y[-1] + 0.05],
    )
    ax2.add_collection(LineCollection(data["outline"][..., :2], colors="#a64e59", linewidths=0.65))
    ax2.add_patch(Rectangle((-3, -2), 6, 4, fill=False, edgecolor=INK, linewidth=0.8, linestyle=(0, (3, 2))))
    ax2.set(xlim=(-3.2, 10.2), ylim=(-6.2, 6.2), xlabel="X (mm)", ylabel="Y (mm)")
    ax2.set_xticks([0, 5, 10])
    ax2.set_yticks([-5, 0, 5])
    ax2.tick_params(labelsize=8, pad=3)
    ax2.annotate(
        "ROI",
        xy=(3, 2),
        xytext=(7, 4.5),
        fontsize=7.5,
        ha="center",
        va="center",
        color=INK,
        arrowprops=dict(arrowstyle="-", color=INK, lw=0.6),
    )
    fig.text(0.68, 0.805, "Mid-plane (Z = 0 mm)", fontsize=9)
    cax = fig.add_axes([0.705, 0.23, 0.205, 0.026])
    cb = fig.colorbar(
        plt.cm.ScalarMappable(norm=colors.Normalize(*CLIM), cmap=CMAP), cax=cax, orientation="horizontal"
    )
    cb.set_ticks([0, 20, 39])
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
    stem = "Fig1b_acoustic_field"
    if dpi == 600:
        fig.savefig(output / (stem + ".png"), dpi=600)
        fig.savefig(
            output / (stem + ".pdf"), dpi=600, metadata={"Title": "B2 three-dimensional acoustic field"}
        )
        fig.savefig(output / (stem + ".svg"), dpi=600)
        with Image.open(output / (stem + ".png")) as im:
            im.convert("RGB").save(output / (stem + ".tif"), compression="tiff_lzw", dpi=(600, 600))
            im.convert("RGB").resize(tuple(round(v * 0.3) for v in im.size), Image.Resampling.LANCZOS).save(
                output / (stem + "_preview.png"), dpi=(180, 180)
            )
    else:
        fig.savefig(output / (stem + "_preview.png"), dpi=dpi)
    plt.close(fig)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    data = load_data()
    audit = make_figure(data, args.output, 180 if args.preview else 600)
    for path, digest in data["protected"].items():
        assert sha(ROOT / path) == digest
    report = dict(
        status="PREVIEW" if args.preview else "PASS",
        case="B2",
        protected_inputs=data["protected"],
        renderer_sha256=sha(__file__),
        frequency_Hz=data["frequency_Hz"],
        roi_mm=ROI.tolist(),
        grid_shape=list(data["pressure"].shape),
        roi_shape=list(data["roi"].shape),
        native_spacing_mm=0.1,
        original_pressure_range_kPa=[float(data["pressure"].min()), float(data["pressure"].max())],
        roi_pressure_range_kPa=[float(data["roi"].min()), float(data["roi"].max())],
        amplitude_definition="abs(complex pressure phasor), not RMS or instantaneous pressure",
        scalar_mapping="Shared linear 0-39 kPa scale; rendered volume colors additionally composite in depth",
        array_mapping_verified=True,
        original_inputs_unchanged=True,
        new_simulation=False,
        rendering=audit,
        reference_style="PZT_control_reference_v6: white, Arial regular, dark axes, muted blue, orthographic -48/26 degrees",
        export="600 dpi raster volume; editable vector annotations and native-value 2D section in PDF/SVG",
    )
    (args.output / "PROVENANCE.json").write_text(json.dumps(report, indent=2) + "\n")
    if not args.preview:
        shutil.copy2(__file__, args.output / Path(__file__).name)
        with h5py.File(args.output / "source_data.h5", "w") as h:
            h.attrs["source_sha256"] = sha(FIELD)
            h.attrs["frequency_Hz"] = data["frequency_Hz"]
            h.attrs["description"] = (
                "Unmodified native nodes; exact crop for volume and exact Z=0 section. Display-only opacity."
            )
            for name, axis in zip("xyz", data["roi_axes"]):
                h.create_dataset(f"roi_axes/{name}_mm", data=axis)
            h.create_dataset("roi_pressure_complex_Pa", data=data["complex_roi"], compression="gzip")
            h.create_dataset("roi_pressure_amplitude_kPa", data=data["roi"], compression="gzip")
            for name, axis in zip("xy", data["axes"][:2]):
                h.create_dataset(f"section_axes/{name}_mm", data=axis)
            h.create_dataset("z0_pressure_complex_Pa", data=data["section_complex"], compression="gzip")
            h.create_dataset("z0_vessel_outline_mm", data=data["outline"], compression="gzip")
            h.create_dataset("display_opacity_kPa_alpha", data=OPACITY)
        (args.output / "CAPTION.txt").write_text(
            "B2 acoustic pressure amplitude at 2.5 MHz, obtained from the saved k-Wave basis superposition. "
            "Left: composite ray-cast volume of |P| in the focal region X=[-3,3], Y=[-2,2], Z=[-1.6,1.6] mm, "
            "using all 61 x 41 x 33 native nodes at 0.1 mm spacing, with 28 and 35 kPa isosurfaces. "
            "The isosurfaces are linearly interpolated from these native data; no geometric or field smoothing is applied. "
            "The layered structures are features of the simulated acoustic field, not bubble surfaces. "
            "Right: native Z=0 mm amplitude section; the dashed rectangle locates the left volume, and the muted-red outline "
            "is the original vessel section. Both panels use a linear 0-39 kPa color scale. P is a complex phasor amplitude, "
            "not RMS pressure. Volume transparency is a display transfer function, not acoustic attenuation: alpha per "
            "0.2 mm at 0/5/12/20/27/39 kPa is 0/0/0.008/0.035/0.15/0.55, with linear interpolation; "
            "values up to 5 kPa are transparent in the volume and remain visible in the quantitative section. "
            "Composited volume colors depend on the view and optical overlap; the right section provides pointwise color values. "
            "Isosurface opacities are 0.38 and 0.82, respectively. Coordinates have equal spatial scale with an orthographic "
            "view at azimuth -48 degrees and elevation 26 degrees. Source data and original solver outputs are unchanged. "
            "The model uses homogeneous water and does not resolve acoustic vessel-wall material interfaces. "
            "No simulation or new physical validation was performed. White background, subdued colors and dark Arial labels "
            "follow PZT_control_reference_v6. PDF/SVG contain a 600 dpi rendered volume and editable vector text.\n"
        )
    print(json.dumps(dict(status=report["status"], output=str(args.output), rendering=audit)), flush=True)


if __name__ == "__main__":
    main()
