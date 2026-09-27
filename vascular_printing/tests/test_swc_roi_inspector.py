from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from inspect_swc_rois import (
    Inspector, Region, TUBE_SIDES, choose_rows, extract_region, geometry,
    latest_run, read_swc, source_for, write_audit,
)


@pytest.fixture
def raw_source(tmp_path):
    path = tmp_path / "sample.swc"
    path.write_text("# Source coordinates and radius in mm\n"
                    "2 3 0 0 0 1 -1\n"
                    "9 3 5 0 0 0.31 2\n"
                    "15 3 10 3 0 0.62 9\n"
                    "20 3 15 3 0 0.31 15\n"
                    "50 3 5 -5 0 0.5 9\n")
    return path, read_swc(path)


def test_region_uses_source_ids_values_and_only_complete_original_edges(raw_source):
    _, source = raw_source
    row = {"global_node_ids": "9;15;20", "anchor_id": "15", "roi_id": "roi"}
    row.update({f"bbox_min_{axis}_um": "-20000" for axis in "xyz"})
    row.update({f"bbox_max_{axis}_um": "20000" for axis in "xyz"})
    region = extract_region(row, source, "mm")
    assert [node.id for node in region.nodes] == [9, 15, 20]
    assert [node.line for node in region.nodes] == [3, 4, 5]
    assert [node.radius for node in region.nodes] == [0.31, 0.62, 0.31]
    assert region.edges == [(9, 15), (15, 20)]
    assert region.nodes[0] is source[9]
    # A stale/mismatched source is rejected rather than silently substituted.
    row["bbox_max_x_um"] = "1000"
    with pytest.raises(ValueError, match="bounds"):
        extract_region(row, source, "mm")


@pytest.mark.parametrize("bad_row", ["2 3 0 0 0 0 -1", "2 3 0 0 0 1 88", "2 3 0 0 0 1", "2.5 3 0 0 0 1 -1"])
def test_bad_source_is_rejected_without_inventing_values(tmp_path, bad_row):
    path = tmp_path / "bad.swc"
    path.write_text(bad_row)
    with pytest.raises(ValueError):
        read_swc(path)


def test_tubes_preserve_source_rings_and_have_no_added_centerline_samples(raw_source):
    _, source = raw_source
    region = Region(9, "roi", list(source.values()), [(node.parent, node.id) for node in source.values() if node.parent != -1])
    line, tube, points, spheres = geometry(region, panel=2)
    assert line.n_lines == 3  # one junction, no new centerline vertices
    assert tube.GetNumberOfPolys() == 0  # no hidden end disks
    for node_id, point, radius in zip(line["source_node_id"], line.points, line["radius"]):
        np.testing.assert_array_equal(point, source[node_id].xyz)
        assert radius == source[node_id].radius
    centers = np.asarray([source[node_id].xyz for node_id in tube["source_node_id"]])
    radii = np.asarray([source[node_id].radius for node_id in tube["source_node_id"]])
    np.testing.assert_allclose(np.linalg.norm(tube.points - centers, axis=1), radii, atol=1e-12, rtol=1e-12)
    assert tube.n_points == line.n_points * TUBE_SIDES
    for mesh in (points, spheres, tube):
        assert set(mesh["source_node_id"]) == set(source)
        assert np.all(mesh["roi_panel"] == 2)
    centers = np.asarray([source[node_id].xyz for node_id in spheres["source_node_id"]])
    np.testing.assert_allclose(np.linalg.norm(spheres.points - centers, axis=1), spheres["radius"], atol=2e-6)


def test_view_switching_and_pick_resolve_the_source_record(raw_source, tmp_path):
    path, source = raw_source
    region = Region(9, "roi", list(source.values()), [(n.parent, n.id) for n in source.values() if n.parent != -1])
    viewer = Inspector(source, [region] * 3, path, "mm", tmp_path, off_screen=True)
    try:
        for mode, expected in enumerate(((1, 0, 0, 0), (0, 1, 1, 0), (0, 0, 0, 1))):
            viewer.set_mode(mode)
            assert tuple(actor.GetVisibility() for actor in viewer.actors[1]) == expected
        tube = geometry(region, 1)[1]
        index = int(np.flatnonzero(tube["source_node_id"] == 15)[0])
        picker = SimpleNamespace(GetDataSet=lambda: tube, GetPointId=lambda: index)
        viewer.picked(tube.points[index], picker)
        readout = viewer.plotter.renderers[4].actors["source_readout"].GetInput()
        assert "node 15 / file line 4" in readout
        assert "r = 0.62 mm" in readout
        assert "ratio=2" in readout
        viewer.toggle_color()
        assert all(actors[0].mapper.GetScalarVisibility() for actors in viewer.actors)
        viewer.reset()
        viewer.plotter.screenshot(str(tmp_path / "smoke.png"))
        assert (tmp_path / "smoke.png").stat().st_size > 10000
    finally:
        viewer.plotter.close()


def test_current_human_three_rois_match_original_swc(tmp_path):
    try:
        run = latest_run()
    except ValueError:
        pytest.skip("Local human sampling outputs are not present")
    rows = choose_rows(run, None)
    path, unit = source_for(run, rows[0])
    if not path.exists():
        pytest.skip("Local BraVa source is not present")
    before = path.read_bytes()
    source = read_swc(path)
    regions = [extract_region(row, source, unit) for row in rows]
    assert len({region.anchor for region in regions}) == 3
    assert [row["roi_id"] for row in choose_rows(run, [r.anchor for r in regions])] == [r.roi_id for r in regions]
    for region in regions:
        _, tube, _, _ = geometry(region, 0)
        centers = np.asarray([source[node_id].xyz for node_id in tube["source_node_id"]])
        radii = np.asarray([source[node_id].radius for node_id in tube["source_node_id"]])
        np.testing.assert_allclose(np.linalg.norm(tube.points - centers, axis=1), radii, atol=1e-11, rtol=1e-11)
    write_audit(tmp_path, run, path, unit, regions)
    assert path.read_bytes() == before
