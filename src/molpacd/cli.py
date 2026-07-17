from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any, Optional, Sequence

from molpacd._version import __version__
from molpacd.capper import add_caps, analyze_structure, remove_caps
from molpacd.io import read_structure, write_structure
from molpacd.models import AnalysisResult, CapOptions, CapResult


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except Exception as exc:  # noqa: BLE001
        if getattr(args, "debug", False):
            raise
        print(f"molpacd: error: {exc}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="molpacd",
        description="Molecular Pore Aperture Cap Designer",
    )
    parser.add_argument("--version", action="version", version=f"molpacd {__version__}")
    parser.add_argument("--debug", action="store_true", help="show tracebacks for errors")
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser("analyze", help="analyze candidate aperture caps")
    analyze.add_argument("input", type=Path)
    _add_analysis_options(analyze)
    analyze.add_argument("--json", action="store_true", help="write machine-readable JSON")
    analyze.set_defaults(func=_cmd_analyze)

    add = subparsers.add_parser("add", help="add molecular aperture caps")
    add.add_argument("input", type=Path)
    add.add_argument("-o", "--output", type=Path, required=True)
    _add_analysis_options(add)
    add.add_argument("--chain", default="Z", help="chain identifier for cap atoms")
    add.add_argument("--resname", help="1-3 character cap residue name")
    add.add_argument("--atom-name", default="O", help="cap atom name")
    add.add_argument("--element", default="O", help="cap atom element")
    add.add_argument("--seed", type=int, help="seed for automatic residue-name generation")
    add.add_argument(
        "--min-clearance",
        type=float,
        default=1.4,
        help="minimum distance from existing atoms before skipping a cap atom",
    )
    add.add_argument(
        "--shared-radius",
        dest="shared_radius",
        action="store_true",
        default=True,
        help="use the larger opening radius for all requested sides",
    )
    add.add_argument(
        "--independent-radius",
        dest="shared_radius",
        action="store_false",
        help="use each opening's own inferred radius",
    )
    add.add_argument("--json", action="store_true", help="write machine-readable JSON")
    add.add_argument(
        "--dry-run", action="store_true", help="report cap design without writing output"
    )
    add.set_defaults(func=_cmd_add)

    remove = subparsers.add_parser("remove", help="remove MolPACD aperture caps")
    remove.add_argument("input", type=Path)
    remove.add_argument("-o", "--output", type=Path, required=True)
    remove.add_argument("--resname", help="cap residue name to remove")
    remove.add_argument("--chain", help="cap chain identifier to remove")
    remove.add_argument("--atom-name", help="cap atom name to remove")
    remove.add_argument(
        "--force",
        action="store_true",
        help="allow removal by residue/chain/atom match when MolPACD metadata is absent",
    )
    remove.add_argument("--json", action="store_true", help="write machine-readable JSON")
    remove.set_defaults(func=_cmd_remove)
    return parser


def _add_analysis_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--axis",
        default="auto",
        help="analysis axis: auto, x, y, z, or a 3-value vector such as 0,0,1",
    )
    parser.add_argument(
        "--selection",
        choices=["ca", "backbone", "all"],
        default="ca",
        help="atoms used to infer the pore axis and openings",
    )
    parser.add_argument(
        "--sides",
        choices=["both", "negative", "positive"],
        default="both",
        help="which opening side or sides to cap",
    )
    parser.add_argument(
        "--spacing", type=float, default=2.8, help="cap lattice spacing in Angstrom"
    )
    parser.add_argument("--radius-scale", type=float, default=1.0, help="scale opening radii")
    parser.add_argument("--window", type=float, default=6.0, help="opening slice width in Angstrom")
    parser.add_argument("--min-atoms", type=int, default=12, help="minimum atoms per opening slice")
    parser.add_argument(
        "--inversion",
        type=float,
        default=0.0,
        help="distance in Angstrom to invert caps into aperture (0 = no inversion)",
    )


def _options_from_args(args: argparse.Namespace) -> CapOptions:
    return CapOptions(
        axis=args.axis,
        selection=args.selection,
        sides=args.sides,
        spacing=args.spacing,
        radius_scale=args.radius_scale,
        window=args.window,
        min_atoms=args.min_atoms,
        chain=getattr(args, "chain", "Z"),
        resname=getattr(args, "resname", None),
        atom_name=getattr(args, "atom_name", "O"),
        element=getattr(args, "element", "O"),
        seed=getattr(args, "seed", None),
        min_clearance=getattr(args, "min_clearance", 1.4),
        shared_radius=getattr(args, "shared_radius", True),
        inversion=getattr(args, "inversion", 0.0),
    )


def _cmd_analyze(args: argparse.Namespace) -> int:
    structure = read_structure(args.input)
    analysis = analyze_structure(structure, _options_from_args(args))
    if args.json:
        print(json.dumps(_analysis_to_dict(analysis), indent=2, sort_keys=True))
        return 0
    _print_analysis(analysis)
    return 0


def _cmd_add(args: argparse.Namespace) -> int:
    structure = read_structure(args.input)
    capped, result = add_caps(structure, _options_from_args(args))
    payload = _cap_result_to_dict(result)
    payload["dry_run"] = args.dry_run
    payload["output"] = str(args.output)

    if not args.dry_run:
        write_structure(capped, args.output)
        payload["wrote"] = str(args.output)

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    for side in result.sides:
        print(
            f"{side.side}: added {side.added_count}/{side.requested_count} atoms "
            f"(skipped {side.skipped_collision_count} collisions), radius {side.radius:.2f}"
        )
    print(
        f"cap residue {result.resname} chain {result.chain}; total added atoms {result.added_count}"
    )
    if not args.dry_run:
        print(f"wrote {args.output}")
    return 0


def _cmd_remove(args: argparse.Namespace) -> int:
    structure = read_structure(args.input)
    decapped, result = remove_caps(
        structure,
        resname=args.resname,
        chain=args.chain,
        atom_name=args.atom_name,
        force=args.force,
    )
    write_structure(decapped, args.output)
    payload = asdict(result)
    payload["output"] = str(args.output)
    payload["wrote"] = str(args.output)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    print(
        f"removed {result.removed_count} atoms matching residue {result.resname}"
        f"{f' chain {result.chain}' if result.chain else ''}"
        f"{f' atom {result.atom_name}' if result.atom_name else ''}"
    )
    print(f"wrote {args.output}")
    return 0


def _print_analysis(analysis: AnalysisResult) -> None:
    print(f"axis: {_format_vector(analysis.axis)}")
    print(f"center: {_format_vector(analysis.center)}")
    for opening in [analysis.negative, analysis.positive]:
        print(
            f"{opening.side}: atoms={opening.atom_count} "
            f"centroid={_format_vector(opening.centroid)} radius={opening.radius:.2f}"
        )


def _analysis_to_dict(analysis: AnalysisResult) -> dict[str, Any]:
    return asdict(analysis)


def _cap_result_to_dict(result: CapResult) -> dict[str, Any]:
    return asdict(result)


def _format_vector(values: Sequence[float]) -> str:
    return "(" + ", ".join(f"{value:.3f}" for value in values) + ")"
