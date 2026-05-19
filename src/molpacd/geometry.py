from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from molpacd.models import AtomRecord, CapOptions, Opening

BACKBONE_NAMES = {"N", "CA", "C", "O"}


def coords_from_atoms(atoms: Sequence[AtomRecord]) -> NDArray[np.float64]:
    if not atoms:
        raise ValueError("no atoms selected for analysis")
    return np.array([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)


def select_atoms(atoms: Iterable[AtomRecord], options: CapOptions) -> List[AtomRecord]:
    selected: List[AtomRecord] = []

    for atom in atoms:
        if options.selection in {"ca", "backbone"} and atom.record != "ATOM":
            continue
        if options.selection == "ca" and atom.name.strip().upper() != "CA":
            continue
        if options.selection == "backbone" and atom.name.strip().upper() not in BACKBONE_NAMES:
            continue
        selected.append(atom)

    if len(selected) < 3:
        raise ValueError(
            f"selection {options.selection!r} produced {len(selected)} atoms; at least 3 required"
        )
    return selected


def normalize_vector(vector: NDArray[np.float64]) -> NDArray[np.float64]:
    norm = float(np.linalg.norm(vector))
    if norm == 0:
        raise ValueError("axis vector must not be zero")
    return vector / norm


def parse_axis(axis: str, coords: NDArray[np.float64]) -> NDArray[np.float64]:
    axis_value = axis.strip().lower()
    if axis_value == "auto":
        centered = coords - np.mean(coords, axis=0)
        covariance = np.cov(centered, rowvar=False)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        vector = eigenvectors[:, int(np.argmax(eigenvalues))]
        largest_component = int(np.argmax(np.abs(vector)))
        if vector[largest_component] < 0:
            vector = -vector
        return normalize_vector(vector)
    if axis_value == "x":
        return np.array([1.0, 0.0, 0.0])
    if axis_value == "y":
        return np.array([0.0, 1.0, 0.0])
    if axis_value == "z":
        return np.array([0.0, 0.0, 1.0])

    parts = axis.replace(",", " ").split()
    if len(parts) != 3:
        raise ValueError("axis must be auto, x, y, z, or a 3-value vector such as 0,0,1")
    try:
        vector = np.array([float(part) for part in parts], dtype=float)
    except ValueError as exc:
        raise ValueError("axis vector must contain numeric values") from exc
    return normalize_vector(vector)


def analyze_openings(
    atoms: Sequence[AtomRecord],
    options: CapOptions,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], Opening, Opening]:
    coords = coords_from_atoms(atoms)
    center = np.mean(coords, axis=0)
    axis = parse_axis(options.axis, coords)
    projections = np.dot(coords - center, axis)

    negative = _opening_from_projection(
        coords=coords,
        center=center,
        axis=axis,
        projections=projections,
        side="negative",
        window=options.window,
        min_atoms=options.min_atoms,
    )
    positive = _opening_from_projection(
        coords=coords,
        center=center,
        axis=axis,
        projections=projections,
        side="positive",
        window=options.window,
        min_atoms=options.min_atoms,
    )
    return axis, center, negative, positive


def _opening_from_projection(
    coords: NDArray[np.float64],
    center: NDArray[np.float64],
    axis: NDArray[np.float64],
    projections: NDArray[np.float64],
    side: str,
    window: float,
    min_atoms: int,
) -> Opening:
    if side == "negative":
        limit = float(np.min(projections)) + window
        indices = np.where(projections <= limit)[0]
        if indices.size < min_atoms:
            indices = np.argsort(projections)[: min(min_atoms, len(projections))]
    else:
        limit = float(np.max(projections)) - window
        indices = np.where(projections >= limit)[0]
        if indices.size < min_atoms:
            indices = np.argsort(projections)[-min(min_atoms, len(projections)) :]

    opening_coords = coords[indices]
    raw_centroid = np.mean(opening_coords, axis=0)
    projection = float(np.dot(raw_centroid - center, axis))
    centroid = center + projection * axis
    offsets = opening_coords - centroid
    axial_offsets = np.outer(np.dot(offsets, axis), axis)
    radial_offsets = offsets - axial_offsets
    distances = np.linalg.norm(radial_offsets, axis=1)
    radius = float(np.mean(distances))
    return Opening(
        side="negative" if side == "negative" else "positive",
        atom_count=int(indices.size),
        centroid=_coord_tuple(centroid),
        radius=radius,
        projection=projection,
    )


def generate_cap_disk(
    centroid: Sequence[float],
    axis: Sequence[float],
    radius: float,
    spacing: float,
    side_sign: float,
) -> NDArray[np.float64]:
    if spacing <= 0:
        raise ValueError("spacing must be greater than zero")
    if radius <= 0:
        raise ValueError("opening radius must be greater than zero")

    center = np.array(centroid, dtype=float)
    cap_direction = normalize_vector(np.array(axis, dtype=float) * side_sign)

    if abs(float(cap_direction[2])) < 0.9:
        perp1 = np.cross(cap_direction, np.array([0.0, 0.0, 1.0]))
    else:
        perp1 = np.cross(cap_direction, np.array([1.0, 0.0, 0.0]))
    perp1 = normalize_vector(perp1)
    perp2 = normalize_vector(np.cross(cap_direction, perp1))

    local_positions: List[Tuple[float, float]] = []
    row_spacing = spacing * np.sqrt(3.0) / 2.0
    y_values = np.arange(-radius, radius + row_spacing, row_spacing)
    for row_index, y_local in enumerate(y_values):
        max_x = np.sqrt(max(radius**2 - y_local**2, 0.0))
        x_offset = 0.5 * spacing if row_index % 2 else 0.0
        x_values = np.arange(-max_x, max_x + spacing, spacing) + x_offset
        for x_local in x_values:
            if x_local**2 + y_local**2 <= radius**2:
                local_positions.append((float(x_local), float(y_local)))

    if not local_positions:
        local_positions.append((0.0, 0.0))

    local = np.array(local_positions, dtype=float)
    local -= np.mean(local, axis=0)
    return np.array([center + x * perp1 + y * perp2 for x, y in local], dtype=float)


def filter_collisions(
    positions: NDArray[np.float64],
    existing_coords: NDArray[np.float64],
    min_clearance: float,
) -> Tuple[NDArray[np.float64], int]:
    if min_clearance <= 0 or existing_coords.size == 0:
        return positions, 0

    kept = []
    skipped = 0
    for position in positions:
        distances = np.linalg.norm(existing_coords - position, axis=1)
        if float(np.min(distances)) < min_clearance:
            skipped += 1
        else:
            kept.append(position)

    if not kept:
        return np.empty((0, 3), dtype=float), skipped
    return np.array(kept, dtype=float), skipped


def _coord_tuple(values: NDArray[np.float64]) -> Tuple[float, float, float]:
    return (float(values[0]), float(values[1]), float(values[2]))
