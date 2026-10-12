"""CPU interface and stored-data checks; no additional GPU job required."""

import h5py
import numpy as np
from .common import ROOT, BASE, case, write
from mevo_mb.fast_sim.f10c_families import control_families


def main():
    records = []
    for key in ["O_1000", "OFF", "S_1000", "O_1000"]:
        cfg = case(key)
        records.append(
            dict(
                sample_id=f"id_{len(records)}",
                positions_m=cfg["environment"]["positions_m"],
                events=cfg["events"],
                split_group="F10A_translated_original_Y32",
            )
        )
    groups = control_families(records)
    assert [len(g["sample_ids"]) for g in groups] == [2, 1, 1]
    assert groups[0]["sample_ids"] == ["id_0", "id_3"]
    assert len(control_families(records, max_batch=1)) == 4
    try:
        control_families(records + [records[0]])
    except ValueError:
        pass
    else:
        raise AssertionError("Duplicate IDs accepted")
    archived = []
    for key, path in [
        ("OFF", BASE / "multirate_v1/Q0_OFF"),
        ("S_1000", BASE / "multirate_v1/Q0_S_1000"),
        ("O_1000", BASE / "batch_v1/B1_batch"),
    ]:
        with (
            h5py.File(path / "trajectory.h5") as h,
            h5py.File(ROOT / "results/fast_pzt/f10a/extended_10_v1" / key / "trajectory.h5") as ref,
        ):
            np.testing.assert_array_equal(h["positions_m"][:], ref["positions_m"][:])
            archived.append(key)
            for group in [
                "source",
                "medium",
                "flow",
                "bubbles",
                "control",
                "trajectory",
                "morphology",
                "quality",
            ]:
                assert group in h
    write(
        BASE / "INTERFACE_CHECKS.json",
        dict(
            status="PASS",
            family_sizes=[len(g["sample_ids"]) for g in groups],
            actual_control_scope="Sequential control families; members in each family use one shared synchronous schedule.",
            sample_identity_and_split_groups="PASS",
            duplicate_ID_rejection="PASS",
            F10A_10us_bitwise=archived,
            async_per_member_control="NOT_IMPLEMENTED",
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
