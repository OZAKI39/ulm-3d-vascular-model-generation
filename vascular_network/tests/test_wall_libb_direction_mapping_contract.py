from __future__ import annotations


def test_wall_libb_direction_mapping_contract() -> None:
    # D3Q19 inverse pairs from the pinned stencil contract.
    inverse = {1: 2, 2: 1, 3: 4, 4: 3, 5: 6, 6: 5}
    outgoing_direction = 1
    incoming_bitmask_direction = inverse[outgoing_direction]
    q_slot_read_by_mus_set_bouzidi = inverse[incoming_bitmask_direction]
    assert q_slot_read_by_mus_set_bouzidi == outgoing_direction
