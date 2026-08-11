#!/usr/bin/env python3
"""Plot the first N minutes of current, OSB, and backup FLT ENU errors.

The current UPD files and the OSB files are compared with the matching
kinematic Fixed files under ``result/bak``.  ENU errors are computed against
the station reference ECEF coordinates used by the full-day test configuration.

Example::

    python scripts/plot_flt_enu_comparison.py \
        --minutes 30 \
        --output build/plots/flt_enu_first_30m.png
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REFERENCE_ECEF = {
    "GODN": (1130760.6931, -4831298.6759, 3994155.1990),
    "HARB": (5084657.6078, 2670325.4787, -2768480.8416),
}

PAIRS = (
    (
        "GODN · DF",
        "GODN",
        "CODEX_FULLDAY_FLT_UPD_DF_GODN.flt",
        "CODEX_FULLDAY_FLT_OSB_DF_GODN.flt",
        "GODN-PPP_kin_DF_Fixed.flt",
    ),
    (
        "GODN · FF",
        "GODN",
        "CODEX_FULLDAY_FLT_UPD_FF_GODN.flt",
        "CODEX_FULLDAY_FLT_OSB_FF_GODN.flt",
        "GODN-PPP_kin_FF_Fixed.flt",
    ),
    (
        "HARB · DF",
        "HARB",
        "CODEX_FULLDAY_FLT_UPD_DF_HARB.flt",
        "CODEX_FULLDAY_FLT_OSB_DF_HARB.flt",
        "HARB-PPP_kin_DF_Fixed.flt",
    ),
    (
        "HARB · FF",
        "HARB",
        "CODEX_FULLDAY_FLT_UPD_FF_HARB.flt",
        "CODEX_FULLDAY_FLT_OSB_FF_HARB.flt",
        "HARB-PPP_kin_FF_Fixed.flt",
    ),
)


@dataclass(frozen=True)
class Row:
    sow: float
    xyz: tuple[float, float, float]
    status: str


def read_flt(path: Path) -> list[Row]:
    """Read the position and ambiguity state fields from a .flt file."""

    rows: list[Row] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            try:
                row = Row(
                    sow=float(fields[0]),
                    xyz=(float(fields[1]), float(fields[2]), float(fields[3])),
                    status=fields[16],
                )
            except (IndexError, ValueError) as exc:
                raise ValueError(f"malformed FLT row at {path}:{line_number}") from exc
            if all(math.isfinite(value) for value in (row.sow, *row.xyz)):
                rows.append(row)
    return sorted(rows, key=lambda row: row.sow)


def ecef_to_geodetic_lat_lon(x: float, y: float, z: float) -> tuple[float, float]:
    """Return WGS-84 geodetic latitude and longitude in radians."""

    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - e2))
    for _ in range(12):
        sin_lat = math.sin(lat)
        n = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        height = p / max(math.cos(lat), 1e-15) - n
        next_lat = math.atan2(z, p * (1.0 - e2 * n / (n + height)))
        if abs(next_lat - lat) < 1e-14:
            lat = next_lat
            break
        lat = next_lat
    return lat, lon


def ecef_to_enu(
    xyz: tuple[float, float, float],
    reference: tuple[float, float, float],
    lat: float,
    lon: float,
) -> tuple[float, float, float]:
    """Rotate an ECEF error vector into the local east/north/up frame."""

    dx, dy, dz = (xyz[index] - reference[index] for index in range(3))
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    east = -sin_lon * dx + cos_lon * dy
    north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    up = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
    return east, north, up


def enu_series(
    rows: Iterable[Row],
    reference: tuple[float, float, float],
    nominal_start_sow: float,
    minutes: float,
) -> dict[float, tuple[float, float, float, str]]:
    """Return elapsed minutes -> (E, N, U, status) within the requested window."""

    lat, lon = ecef_to_geodetic_lat_lon(*reference)
    result: dict[float, tuple[float, float, float, str]] = {}
    for row in rows:
        elapsed = (row.sow - nominal_start_sow) / 60.0
        if elapsed < 0.0 or elapsed > minutes + 1e-9:
            continue
        east, north, up = ecef_to_enu(row.xyz, reference, lat, lon)
        result[round(elapsed, 6)] = (east, north, up, row.status)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("sample_data/PPPFLT_2023305/result"),
        help="result directory containing current UPD .flt files and bak/",
    )
    parser.add_argument(
        "--osb-root",
        type=Path,
        default=Path("sample_data/PPPFLT_2023305_OSB/result"),
        help="result directory containing OSB .flt files",
    )
    parser.add_argument("--minutes", type=float, default=30.0, help="plot window in minutes")
    parser.add_argument(
        "--start-sow",
        type=float,
        default=None,
        help="nominal start SOW; defaults to the first current epoch",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("build/plots/flt_enu_first_30m.png"),
        help="PNG output path",
    )
    return parser


def plot_comparison(
    root: Path,
    osb_root: Path,
    output: Path,
    minutes: float,
    start_sow: float | None,
) -> None:
    if minutes <= 0.0 or not math.isfinite(minutes):
        raise ValueError("minutes must be a positive finite number")

    loaded: list[
        tuple[
            str,
            str,
            dict[float, tuple[float, float, float, str]],
            dict[float, tuple[float, float, float, str]],
            dict[float, tuple[float, float, float, str]],
        ]
    ] = []
    inferred_start: float | None = start_sow
    for title, site, current_name, osb_name, backup_name in PAIRS:
        current_path = root / current_name
        osb_path = osb_root / osb_name
        backup_path = root / "bak" / backup_name
        if not current_path.is_file():
            raise FileNotFoundError(current_path)
        if not osb_path.is_file():
            raise FileNotFoundError(osb_path)
        if not backup_path.is_file():
            raise FileNotFoundError(backup_path)
        current_rows = read_flt(current_path)
        if inferred_start is None:
            if not current_rows:
                raise ValueError(f"no data rows in {current_path}")
            inferred_start = current_rows[0].sow
        reference = REFERENCE_ECEF[site]
        loaded.append(
            (
                title,
                site,
                enu_series(current_rows, reference, inferred_start, minutes),
                enu_series(read_flt(osb_path), reference, inferred_start, minutes),
                enu_series(read_flt(backup_path), reference, inferred_start, minutes),
            )
        )

    assert inferred_start is not None
    components = (("E", 0), ("N", 1), ("U", 2))
    limits: dict[int, tuple[float, float]] = {}
    for _, _, current, osb, backup in loaded:
        for component, index in components:
            values = [item[index] for series in (current, osb, backup) for item in series.values()]
            if not values:
                continue
            low, high = min(values), max(values)
            span = max(high - low, 0.02)
            limits[index] = (
                min(limits.get(index, (float("inf"), float("inf")))[0], low - span * 0.08),
                max(limits.get(index, (float("-inf"), float("-inf")))[1], high + span * 0.08),
            )

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})
    figure, axes = plt.subplots(4, 3, figsize=(15, 12), sharex=True)
    figure.suptitle(
        f"First {minutes:g} minutes ENU error: current FLT vs OSB FLT vs backup Fixed",
        fontsize=15,
        fontweight="bold",
        y=0.995,
    )
    colors = {"current": "#D55E00", "osb": "#7B3294", "backup": "#0072B2"}
    line_styles = {"current": ("-", "o"), "osb": ("-.", "s"), "backup": ("--", "^")}
    axes[-1, 0].set_xlabel("Elapsed time (min)")
    axes[-1, 1].set_xlabel("Elapsed time (min)")
    axes[-1, 2].set_xlabel("Elapsed time (min)")

    for row_index, (title, _, current, osb, backup) in enumerate(loaded):
        series_by_name = {"current": current, "osb": osb, "backup": backup}
        for col_index, (component, value_index) in enumerate(components):
            axis = axes[row_index, col_index]
            for label, values in series_by_name.items():
                ordered = sorted(values.items())
                times = [time for time, _ in ordered]
                errors = [item[value_index] for _, item in ordered]
                linestyle, marker = line_styles[label]
                axis.plot(
                    times,
                    errors,
                    color=colors[label],
                    linestyle=linestyle,
                    marker=marker,
                    markersize=3,
                    linewidth=1.2,
                    label=label,
                    alpha=0.95,
                )

            common_times = sorted(set(current) & set(osb) & set(backup))
            mismatch_times = [
                time
                for time in common_times
                if len({current[time][3], osb[time][3], backup[time][3]}) > 1
            ]
            if mismatch_times:
                axis.scatter(
                    mismatch_times,
                    [current[time][value_index] for time in mismatch_times],
                    color="#000000",
                    marker="x",
                    s=28,
                    linewidths=1.0,
                    zorder=5,
                    label="status disagreement" if col_index == 0 else "_nolegend_",
                )

            axis.axhline(0.0, color="#777777", linewidth=0.7, alpha=0.7)
            axis.grid(True, linewidth=0.35, alpha=0.35)
            axis.set_xlim(0.0, minutes)
            if value_index in limits:
                axis.set_ylim(*limits[value_index])
            axis.set_ylabel(f"{component} error (m)")
            axis.set_title(f"{title} · {component}", loc="left", fontweight="bold")
            axis.tick_params(axis="both", labelsize=8)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.968), ncol=3, frameon=False)
    figure.text(
        0.5,
        0.008,
        "ENU = solution ECEF minus station reference, rotated to local East/North/Up; x marks indicate status disagreement among the three files.",
        ha="center",
        fontsize=9,
    )
    figure.tight_layout(rect=(0.03, 0.035, 0.99, 0.945))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    print(f"wrote {output.resolve()}")


def main() -> None:
    args = build_parser().parse_args()
    plot_comparison(args.root, args.osb_root, args.output, args.minutes, args.start_sow)


if __name__ == "__main__":
    main()
