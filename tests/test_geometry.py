from __future__ import annotations

import numpy as np

from molpacd.geometry import filter_collisions, generate_cap_disk


def test_generate_cap_disk_is_planar_and_centered() -> None:
    positions = generate_cap_disk(
        centroid=(1.0, 2.0, 3.0),
        axis=(0.0, 0.0, 1.0),
        radius=2.0,
        spacing=1.0,
        side_sign=1.0,
    )

    assert positions.shape[1] == 3
    assert np.allclose(positions[:, 2], 3.0)
    assert np.allclose(np.mean(positions, axis=0), np.array([1.0, 2.0, 3.0]))


def test_filter_collisions_skips_positions_too_close_to_existing_atoms() -> None:
    positions = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=float)
    existing = np.array([[0.5, 0.0, 0.0]], dtype=float)

    accepted, skipped = filter_collisions(positions, existing, min_clearance=1.0)

    assert skipped == 1
    assert accepted.tolist() == [[5.0, 0.0, 0.0]]
