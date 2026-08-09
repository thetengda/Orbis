#!/usr/bin/env python3
"""Analyze GREAT-PIFGO FGO/FLT solution files over a complete data span.

The script intentionally uses only the Python standard library so it can be
reused on build and test machines without installing scientific packages.
It produces a Markdown report on stdout and can additionally write Markdown
and JSON artifacts.

Example:

  python scripts/analyze_fgo_solution.py \
    --site GODN \
    --reference 1130760.6931 -4831298.6759 3994155.1990 \
    --nominal-start-sow 259200 \
    --case NONE=sample_data/.../GODN-NONE.fgo \
    --case PARAMETER=sample_data/.../GODN-PARAMETER.fgo \
    --case CONSTRAINT=sample_data/.../GODN-CONSTRAINT.fgo \
    --baseline sample_data/.../GODN-FLT.flt \
    --markdown result/GODN-analysis.md --json result/GODN-analysis.json

Definitions used in every report:

* "threshold" means both horizontal error and absolute vertical error satisfy
  their configured limits.
* "sustained convergence" is the first complete, gap-free sample window for
  which every solution satisfies the threshold.
* "permanent convergence" is the first epoch after the final threshold
  violation. Continuity diagnostics must also pass before treating it as a
  valid uninterrupted result.
* "correct Fixed" means the file reports Fixed and the truth threshold is met;
  it is intentionally different from the first reported Fixed epoch.
* E/N/U bias and standard deviation retain their signs. Their MAE and error
  quantiles use absolute component errors. Horizontal and 3D values are
  non-negative magnitudes. Percentiles use linear interpolation between the
  nearest ordered samples.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class SolutionRow:
    sow: float
    x: float
    y: float
    z: float
    nsat: int
    status: str
    line_number: int

    @property
    def xyz(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass(frozen=True)
class ErrorRow:
    sow: float
    east: float
    north: float
    up: float
    horizontal: float
    three_d: float
    nsat: int
    status: str
    correct: bool


def parse_case(value: str) -> Tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("case must be LABEL=PATH")
    label, path = value.split("=", 1)
    label = label.strip()
    if not label or not path.strip():
        raise argparse.ArgumentTypeError("case must contain a non-empty label and path")
    return label, Path(path.strip())


def read_solution(path: Path) -> Tuple[List[SolutionRow], Dict[str, int]]:
    rows: List[SolutionRow] = []
    diagnostics = {"data_lines": 0, "malformed_lines": 0, "nonfinite_xyz": 0}
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            diagnostics["data_lines"] += 1
            fields = line.split()
            try:
                row = SolutionRow(
                    sow=float(fields[0]),
                    x=float(fields[1]),
                    y=float(fields[2]),
                    z=float(fields[3]),
                    nsat=int(fields[13]),
                    status=fields[16],
                    line_number=line_number,
                )
            except (IndexError, ValueError):
                diagnostics["malformed_lines"] += 1
                continue
            if not all(math.isfinite(value) for value in row.xyz):
                diagnostics["nonfinite_xyz"] += 1
            rows.append(row)
    rows.sort(key=lambda row: row.sow)
    return rows, diagnostics


def ecef_to_geodetic_lat_lon(x: float, y: float, z: float) -> Tuple[float, float]:
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


def ecef_delta_to_enu(
    xyz: Tuple[float, float, float],
    reference: Tuple[float, float, float],
    lat: float,
    lon: float,
) -> Tuple[float, float, float]:
    dx = xyz[0] - reference[0]
    dy = xyz[1] - reference[1]
    dz = xyz[2] - reference[2]
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    east = -sin_lon * dx + cos_lon * dy
    north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    up = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
    return east, north, up


def build_errors(
    rows: Sequence[SolutionRow],
    reference: Tuple[float, float, float],
    horizontal_threshold: float,
    vertical_threshold: float,
) -> List[ErrorRow]:
    lat, lon = ecef_to_geodetic_lat_lon(*reference)
    errors: List[ErrorRow] = []
    for row in rows:
        if not all(math.isfinite(value) for value in row.xyz):
            continue
        east, north, up = ecef_delta_to_enu(row.xyz, reference, lat, lon)
        horizontal = math.hypot(east, north)
        three_d = math.sqrt(horizontal * horizontal + up * up)
        errors.append(
            ErrorRow(
                sow=row.sow,
                east=east,
                north=north,
                up=up,
                horizontal=horizontal,
                three_d=three_d,
                nsat=row.nsat,
                status=row.status,
                correct=(horizontal <= horizontal_threshold and abs(up) <= vertical_threshold),
            )
        )
    return errors


def percentile(values: Sequence[float], fraction: float) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def rms(values: Sequence[float]) -> Optional[float]:
    if not values:
        return None
    return math.sqrt(sum(value * value for value in values) / len(values))


def sample_std(values: Sequence[float]) -> Optional[float]:
    if len(values) < 2:
        return 0.0 if values else None
    return statistics.stdev(values)


def component_statistics(
    values: Sequence[float], *, signed_component: bool
) -> Dict[str, Optional[float]]:
    magnitudes = [abs(value) for value in values]
    return {
        "count": len(values),
        "bias_m": statistics.fmean(values) if values and signed_component else None,
        "std_m": sample_std(values),
        "rmse_m": rms(values),
        "mae_m": statistics.fmean(magnitudes) if magnitudes else None,
        "p50_abs_m": percentile(magnitudes, 0.50),
        "p68_abs_m": percentile(magnitudes, 0.68),
        "p90_abs_m": percentile(magnitudes, 0.90),
        "p95_abs_m": percentile(magnitudes, 0.95),
        "p99_abs_m": percentile(magnitudes, 0.99),
        "max_abs_m": max(magnitudes) if magnitudes else None,
    }


def accuracy_metrics(errors: Sequence[ErrorRow]) -> Dict[str, object]:
    east = [row.east for row in errors]
    north = [row.north for row in errors]
    up = [row.up for row in errors]
    horizontal = [row.horizontal for row in errors]
    three_d = [row.three_d for row in errors]
    fixed_count = sum(row.status.lower() == "fixed" for row in errors)
    correct_count = sum(row.correct for row in errors)
    correct_fixed_count = sum(
        row.correct and row.status.lower() == "fixed" for row in errors
    )
    return {
        "count": len(errors),
        "east_bias_m": statistics.fmean(east) if east else None,
        "north_bias_m": statistics.fmean(north) if north else None,
        "up_bias_m": statistics.fmean(up) if up else None,
        "east_std_m": sample_std(east),
        "north_std_m": sample_std(north),
        "up_std_m": sample_std(up),
        "horizontal_rms_m": rms(horizontal),
        "vertical_rms_m": rms(up),
        "three_d_mean_m": statistics.fmean(three_d) if three_d else None,
        "three_d_mae_m": statistics.fmean(three_d) if three_d else None,
        "three_d_rms_m": rms(three_d),
        "three_d_p50_m": percentile(three_d, 0.50),
        "three_d_p68_m": percentile(three_d, 0.68),
        "three_d_p90_m": percentile(three_d, 0.90),
        "three_d_p95_m": percentile(three_d, 0.95),
        "three_d_p99_m": percentile(three_d, 0.99),
        "three_d_max_m": max(three_d) if three_d else None,
        "nsat_median": statistics.median(row.nsat for row in errors) if errors else None,
        "nsat_min": min((row.nsat for row in errors), default=None),
        "fixed_count": fixed_count,
        "fixed_fraction": fixed_count / len(errors) if errors else None,
        "correct_count": correct_count,
        "correct_fraction": correct_count / len(errors) if errors else None,
        "correct_fixed_count": correct_fixed_count,
        "correct_fixed_fraction": correct_fixed_count / len(errors) if errors else None,
        "detailed": {
            "east": component_statistics(east, signed_component=True),
            "north": component_statistics(north, signed_component=True),
            "up": component_statistics(up, signed_component=True),
            "horizontal": component_statistics(horizontal, signed_component=False),
            "three_d": component_statistics(three_d, signed_component=False),
        },
    }


def continuity_metrics(
    rows: Sequence[SolutionRow],
    interval: float,
    nominal_start: float,
    hours: float,
) -> Dict[str, object]:
    duplicate_epochs: List[float] = []
    gap_intervals: List[Dict[str, float]] = []
    for previous, current in zip(rows, rows[1:]):
        delta = current.sow - previous.sow
        if abs(delta) < 1e-8:
            duplicate_epochs.append(current.sow)
        elif abs(delta - interval) > 1e-6:
            gap_intervals.append({"after_sow": previous.sow, "before_sow": current.sow, "delta_s": delta})
    nominal_end = nominal_start + hours * 3600.0
    expected_rows = max(0, int(round(hours * 3600.0 / interval)) - 1)
    expected_epochs = [nominal_start + interval * index for index in range(1, expected_rows + 1)]
    observed_epoch_keys = {round(row.sow, 6) for row in rows}
    missing_epochs = [
        sow for sow in expected_epochs if round(sow, 6) not in observed_epoch_keys
    ]
    out_of_span_epochs = [
        row.sow
        for row in rows
        if row.sow < nominal_start + interval - 1e-6
        or row.sow > nominal_end - interval + 1e-6
    ]
    return {
        "rows": len(rows),
        "expected_rows": expected_rows,
        "coverage_fraction": len(rows) / expected_rows if expected_rows else None,
        "start_sow": rows[0].sow if rows else None,
        "end_sow": rows[-1].sow if rows else None,
        "start_delay_s": (
            rows[0].sow - (nominal_start + interval) if rows else None
        ),
        "end_shortfall_s": (
            (nominal_end - interval) - rows[-1].sow if rows else None
        ),
        "duplicate_epochs": duplicate_epochs,
        "gap_intervals": gap_intervals,
        "nominal_grid_missing_epochs": missing_epochs,
        "out_of_nominal_span_epochs": out_of_span_epochs,
        "continuous": bool(rows) and not duplicate_epochs and not gap_intervals,
        "complete_nominal_span": (
            bool(rows)
            and not duplicate_epochs
            and not gap_intervals
            and not missing_epochs
            and not out_of_span_epochs
        ),
    }


def is_consecutive(first: ErrorRow, second: ErrorRow, interval: float) -> bool:
    return abs((second.sow - first.sow) - interval) <= 1e-6


def first_matching_window(
    errors: Sequence[ErrorRow],
    window_samples: int,
    interval: float,
    predicate,
) -> Optional[float]:
    if window_samples <= 0:
        return None
    for start in range(0, len(errors) - window_samples + 1):
        window = errors[start : start + window_samples]
        if any(not is_consecutive(a, b, interval) for a, b in zip(window, window[1:])):
            continue
        if all(predicate(row) for row in window):
            return window[0].sow
    return None


def permanent_matching_epoch(errors: Sequence[ErrorRow], predicate) -> Optional[float]:
    if not errors:
        return None
    last_bad = -1
    for index, row in enumerate(errors):
        if not predicate(row):
            last_bad = index
    return errors[last_bad + 1].sow if last_bad + 1 < len(errors) else None


def convergence_metrics(
    errors: Sequence[ErrorRow], window_samples: int, interval: float
) -> Dict[str, Optional[float]]:
    first_correct = next((row.sow for row in errors if row.correct), None)
    first_fixed = next((row.sow for row in errors if row.status.lower() == "fixed"), None)
    first_correct_fixed = next(
        (row.sow for row in errors if row.correct and row.status.lower() == "fixed"), None
    )
    return {
        "first_threshold_epoch": first_correct,
        "first_sustained_convergence_epoch": first_matching_window(
            errors, window_samples, interval, lambda row: row.correct
        ),
        "permanent_convergence_epoch": permanent_matching_epoch(errors, lambda row: row.correct),
        "first_fixed_epoch": first_fixed,
        "first_correct_fixed_epoch": first_correct_fixed,
        "first_sustained_correct_fixed_epoch": first_matching_window(
            errors,
            window_samples,
            interval,
            lambda row: row.correct and row.status.lower() == "fixed",
        ),
    }


def group_events(
    events: Sequence[ErrorRow], interval: float, kind: str
) -> List[Dict[str, object]]:
    if not events:
        return []
    groups: List[List[ErrorRow]] = [[events[0]]]
    for event in events[1:]:
        if is_consecutive(groups[-1][-1], event, interval):
            groups[-1].append(event)
        else:
            groups.append([event])
    result: List[Dict[str, object]] = []
    for group in groups:
        result.append(
            {
                "kind": kind,
                "start_sow": group[0].sow,
                "end_sow": group[-1].sow,
                "samples": len(group),
                "max_horizontal_m": max(row.horizontal for row in group),
                "max_abs_vertical_m": max(abs(row.up) for row in group),
                "max_three_d_m": max(row.three_d for row in group),
                "min_nsat": min(row.nsat for row in group),
                "statuses": sorted(set(row.status for row in group)),
            }
        )
    return result


def anomaly_metrics(
    errors: Sequence[ErrorRow],
    convergence: Dict[str, Optional[float]],
    interval: float,
    jump_threshold: float,
    three_d_anomaly_threshold: float,
    min_nsat: int,
) -> Dict[str, object]:
    sustained = convergence["first_sustained_convergence_epoch"]
    after_convergence = [row for row in errors if sustained is not None and row.sow >= sustained]
    threshold_exceedances = [row for row in after_convergence if not row.correct]
    three_d_exceedances = [
        row for row in after_convergence if row.three_d > three_d_anomaly_threshold
    ]
    low_satellites = [row for row in errors if row.nsat < min_nsat]

    jumps: List[Dict[str, float]] = []
    for previous, current in zip(errors, errors[1:]):
        if not is_consecutive(previous, current, interval):
            continue
        jump = math.sqrt(
            (current.east - previous.east) ** 2
            + (current.north - previous.north) ** 2
            + (current.up - previous.up) ** 2
        )
        if jump > jump_threshold and (sustained is None or previous.sow >= sustained):
            jumps.append({"from_sow": previous.sow, "to_sow": current.sow, "jump_m": jump})

    status_transitions: List[Dict[str, object]] = []
    for previous, current in zip(errors, errors[1:]):
        if previous.status != current.status:
            status_transitions.append(
                {"sow": current.sow, "from": previous.status, "to": current.status}
            )

    return {
        "post_convergence_threshold_intervals": group_events(
            threshold_exceedances, interval, "threshold_exceedance"
        ),
        "post_convergence_3d_anomaly_intervals": group_events(
            three_d_exceedances, interval, "three_d_anomaly"
        ),
        "low_satellite_intervals": group_events(low_satellites, interval, "low_satellite"),
        "coordinate_jumps": jumps,
        "status_transitions": status_transitions,
        "baseline_delta_exceedances": [],
        "largest_error_epochs": [
            {
                "sow": row.sow,
                "horizontal_m": row.horizontal,
                "vertical_m": row.up,
                "three_d_m": row.three_d,
                "nsat": row.nsat,
                "status": row.status,
            }
            for row in sorted(errors, key=lambda item: item.three_d, reverse=True)[:10]
        ],
    }


def baseline_metrics(
    rows: Sequence[SolutionRow], baseline_by_epoch: Dict[float, SolutionRow]
) -> Tuple[Dict[str, Optional[float]], List[Dict[str, float]]]:
    deltas: List[Tuple[float, float]] = []
    for row in rows:
        baseline = baseline_by_epoch.get(row.sow)
        if baseline is None:
            continue
        if not all(math.isfinite(value) for value in row.xyz + baseline.xyz):
            continue
        deltas.append((row.sow, math.dist(row.xyz, baseline.xyz)))
    values = [value for _, value in deltas]
    return (
        {
            "matched_epochs": len(values),
            "three_d_mean_m": statistics.fmean(values) if values else None,
            "three_d_mae_m": statistics.fmean(values) if values else None,
            "three_d_rms_m": rms(values),
            "three_d_p50_m": percentile(values, 0.50),
            "three_d_p68_m": percentile(values, 0.68),
            "three_d_p90_m": percentile(values, 0.90),
            "three_d_p95_m": percentile(values, 0.95),
            "three_d_p99_m": percentile(values, 0.99),
            "three_d_max_m": max(values) if values else None,
        },
        [{"sow": sow, "delta_m": value} for sow, value in deltas],
    )


def group_scalar_events(
    events: Sequence[Dict[str, float]], interval: float, value_key: str, kind: str
) -> List[Dict[str, object]]:
    if not events:
        return []
    ordered = sorted(events, key=lambda item: item["sow"])
    groups: List[List[Dict[str, float]]] = [[ordered[0]]]
    for event in ordered[1:]:
        if abs((event["sow"] - groups[-1][-1]["sow"]) - interval) <= 1e-6:
            groups[-1].append(event)
        else:
            groups.append([event])
    return [
        {
            "kind": kind,
            "start_sow": group[0]["sow"],
            "end_sow": group[-1]["sow"],
            "samples": len(group),
            "max_value": max(item[value_key] for item in group),
        }
        for group in groups
    ]


def hourly_metrics(
    errors: Sequence[ErrorRow], nominal_start: float, hours: float
) -> List[Dict[str, object]]:
    result: List[Dict[str, object]] = []
    for hour in range(int(math.ceil(hours))):
        begin = nominal_start + hour * 3600.0
        end = begin + 3600.0
        subset = [row for row in errors if begin <= row.sow < end]
        if not subset:
            continue
        metrics = accuracy_metrics(subset)
        result.append(
            {
                "hour": hour,
                "start_sow": begin,
                "end_sow": end,
                **metrics,
            }
        )
    return result


def analyze_case(
    label: str,
    path: Path,
    rows: Sequence[SolutionRow],
    diagnostics: Dict[str, int],
    reference: Tuple[float, float, float],
    args: argparse.Namespace,
    baseline_by_epoch: Dict[float, SolutionRow],
) -> Dict[str, object]:
    errors = build_errors(rows, reference, args.horizontal_threshold, args.vertical_threshold)
    convergence = convergence_metrics(errors, args.window_samples, args.interval)
    baseline_summary, baseline_deltas = baseline_metrics(rows, baseline_by_epoch)
    anomalies = anomaly_metrics(
        errors,
        convergence,
        args.interval,
        args.jump_threshold,
        args.three_d_anomaly_threshold,
        args.min_nsat,
    )
    warm_start = args.nominal_start_sow + args.warmup_minutes * 60.0
    warm_errors = [row for row in errors if row.sow >= warm_start]
    warm_rows = [row for row in rows if row.sow >= warm_start]
    baseline_warm_summary, _ = baseline_metrics(warm_rows, baseline_by_epoch)
    sustained = convergence["first_sustained_convergence_epoch"]
    permanent = convergence["permanent_convergence_epoch"]
    sustained_errors = [
        row for row in errors if sustained is not None and row.sow >= sustained
    ]
    permanent_errors = [
        row for row in errors if permanent is not None and row.sow >= permanent
    ]
    baseline_exceedances = [
        item
        for item in baseline_deltas
        if item["delta_m"] > args.baseline_delta_threshold
        and (sustained is None or item["sow"] >= sustained)
    ]
    anomalies["baseline_delta_exceedances"] = baseline_exceedances
    anomalies["baseline_delta_intervals"] = group_scalar_events(
        baseline_exceedances, args.interval, "delta_m", "baseline_delta"
    )
    anomalies["largest_error_epochs_after_warmup"] = [
        {
            "sow": row.sow,
            "horizontal_m": row.horizontal,
            "vertical_m": row.up,
            "three_d_m": row.three_d,
            "nsat": row.nsat,
            "status": row.status,
        }
        for row in sorted(warm_errors, key=lambda item: item.three_d, reverse=True)[:10]
    ]
    anomalies["largest_baseline_delta_epochs"] = sorted(
        baseline_deltas, key=lambda item: item["delta_m"], reverse=True
    )[:10]
    return {
        "label": label,
        "path": str(path),
        "parse": diagnostics,
        "continuity": continuity_metrics(
            rows, args.interval, args.nominal_start_sow, args.hours
        ),
        "accuracy_all": accuracy_metrics(errors),
        "accuracy_after_warmup": accuracy_metrics(warm_errors),
        "accuracy_after_sustained_convergence": accuracy_metrics(sustained_errors),
        "accuracy_after_permanent_convergence": accuracy_metrics(permanent_errors),
        "warmup_start_sow": warm_start,
        "convergence": convergence,
        "baseline_comparison_all": baseline_summary,
        "baseline_comparison_after_warmup": baseline_warm_summary,
        "hourly": hourly_metrics(errors, args.nominal_start_sow, args.hours),
        "anomalies": anomalies,
    }


def format_number(value: object, digits: int = 4) -> str:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            return str(value)
        return f"{value:.{digits}f}"
    return str(value)


def format_epoch(value: Optional[float], nominal_start: float) -> str:
    if value is None:
        return "—"
    minutes = (value - nominal_start) / 60.0
    return f"{value:.0f} ({minutes:.1f} min)"


def markdown_report(report: Dict[str, object]) -> str:
    config = report["configuration"]
    nominal_start = float(config["nominal_start_sow"])
    lines = [
        f"# {report['site']} FGO/FLT 全天分析",
        "",
        (
            f"参考坐标 ECEF: `{config['reference_ecef_m']}`；期望采样间隔 "
            f"{config['interval_s']} s；正确门限为水平 ≤ {config['horizontal_threshold_m']} m、"
            f"|高程| ≤ {config['vertical_threshold_m']} m；持续窗口 "
            f"{config['window_minutes']} min（{config['window_samples']} 个连续样本）。"
        ),
        "",
        "## 总览",
        "",
        "| 模式 | 行数/期望/连续 | 全段 3D RMS | 热启动后 3D RMS / MAE / p95 / p99 / max | Fixed | 热启动正确率 | 热启动正确 Fixed 率 | 首次达标 | 首次持续收敛 | 永久收敛 | 首次 Fixed | 首次正确 Fixed | 热启动后相对 FLT RMS |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case in report["cases"]:
        all_accuracy = case["accuracy_all"]
        warm = case["accuracy_after_warmup"]
        convergence = case["convergence"]
        baseline = case["baseline_comparison_after_warmup"]
        fixed_fraction = all_accuracy["fixed_fraction"]
        lines.append(
            "| {label} | {rows}/{continuous} | {all_rms} | {warm_rms} / {mae} / {p95} / {p99} / {maxv} | "
            "{fixed} | {correct} | {correct_fixed} | {first} | {sustained} | {permanent} | {first_fixed} | "
            "{first_correct_fixed} | {baseline_rms} |".format(
                label=case["label"],
                rows=(
                    f"{case['continuity']['rows']}/"
                    f"{case['continuity']['expected_rows']}"
                ),
                continuous=(
                    "complete"
                    if case["continuity"]["complete_nominal_span"]
                    else ("internal yes" if case["continuity"]["continuous"] else "no")
                ),
                all_rms=format_number(all_accuracy["three_d_rms_m"]),
                warm_rms=format_number(warm["three_d_rms_m"]),
                mae=format_number(warm["three_d_mae_m"]),
                p95=format_number(warm["three_d_p95_m"]),
                p99=format_number(warm["three_d_p99_m"]),
                maxv=format_number(warm["three_d_max_m"]),
                fixed=(
                    f"{100.0 * fixed_fraction:.2f}% "
                    f"({all_accuracy['fixed_count']}/{all_accuracy['count']})"
                    if fixed_fraction is not None
                    else "—"
                ),
                correct=(
                    f"{100.0 * warm['correct_fraction']:.2f}% "
                    f"({warm['correct_count']}/{warm['count']})"
                    if warm["correct_fraction"] is not None
                    else "—"
                ),
                correct_fixed=(
                    f"{100.0 * warm['correct_fixed_fraction']:.2f}% "
                    f"({warm['correct_fixed_count']}/{warm['count']})"
                    if warm["correct_fixed_fraction"] is not None
                    else "—"
                ),
                first=format_epoch(convergence["first_threshold_epoch"], nominal_start),
                sustained=format_epoch(convergence["first_sustained_convergence_epoch"], nominal_start),
                permanent=format_epoch(convergence["permanent_convergence_epoch"], nominal_start),
                first_fixed=format_epoch(convergence["first_fixed_epoch"], nominal_start),
                first_correct_fixed=format_epoch(convergence["first_correct_fixed_epoch"], nominal_start),
                baseline_rms=format_number(baseline["three_d_rms_m"]),
            )
        )

    lines.extend(["", "## 分模式统计与异常", ""])
    for case in report["cases"]:
        lines.extend(
            [
                f"### {case['label']}",
                "",
                "| 时段 | E/N/U bias | E/N/U std | H RMS | U RMS | 3D RMS / p95 / max | NSat median/min |",
                "|---|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for name, accuracy in (
            ("全段", case["accuracy_all"]),
            (f"SOW ≥ {case['warmup_start_sow']:.0f}", case["accuracy_after_warmup"]),
        ):
            lines.append(
                "| {name} | {eb}/{nb}/{ub} | {es}/{ns}/{us} | {hr} | {ur} | "
                "{r3}/{p95}/{mx} | {med}/{mn} |".format(
                    name=name,
                    eb=format_number(accuracy["east_bias_m"]),
                    nb=format_number(accuracy["north_bias_m"]),
                    ub=format_number(accuracy["up_bias_m"]),
                    es=format_number(accuracy["east_std_m"]),
                    ns=format_number(accuracy["north_std_m"]),
                    us=format_number(accuracy["up_std_m"]),
                    hr=format_number(accuracy["horizontal_rms_m"]),
                    ur=format_number(accuracy["vertical_rms_m"]),
                    r3=format_number(accuracy["three_d_rms_m"]),
                    p95=format_number(accuracy["three_d_p95_m"]),
                    mx=format_number(accuracy["three_d_max_m"]),
                    med=format_number(accuracy["nsat_median"], 1),
                    mn=format_number(accuracy["nsat_min"], 0),
                )
            )

        lines.extend(
            [
                "",
                (
                    "下表误差单位均为 m；E/N/U 分位数和最大值按绝对误差统计；H/3D 本身为误差幅值。"
                    "bias 仅对有符号的 E/N/U 定义。"
                ),
                "",
                "| 时段 | 分量 | n | bias | std | RMSE | MAE | P50 | P68 | P90 | P95 | P99 | max |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for period_name, accuracy in (
            ("全段", case["accuracy_all"]),
            (f"SOW ≥ {case['warmup_start_sow']:.0f}", case["accuracy_after_warmup"]),
            (
                "首次持续收敛后",
                case["accuracy_after_sustained_convergence"],
            ),
            (
                "永久收敛后",
                case["accuracy_after_permanent_convergence"],
            ),
        ):
            for component_name, stats in accuracy["detailed"].items():
                lines.append(
                    "| {period} | {component} | {count} | {bias} | {std} | {rmse} | "
                    "{mae} | {p50} | {p68} | {p90} | {p95} | {p99} | {maxv} |".format(
                        period=period_name,
                        component=component_name,
                        count=stats["count"],
                        bias=format_number(stats["bias_m"]),
                        std=format_number(stats["std_m"]),
                        rmse=format_number(stats["rmse_m"]),
                        mae=format_number(stats["mae_m"]),
                        p50=format_number(stats["p50_abs_m"]),
                        p68=format_number(stats["p68_abs_m"]),
                        p90=format_number(stats["p90_abs_m"]),
                        p95=format_number(stats["p95_abs_m"]),
                        p99=format_number(stats["p99_abs_m"]),
                        maxv=format_number(stats["max_abs_m"]),
                    )
                )

        baseline = case["baseline_comparison_after_warmup"]
        lines.extend(
            [
                "",
                "| 热启动后同历元相对 FLT | 匹配历元 | MAE | RMS | P50 | P68 | P90 | P95 | P99 | max |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
                (
                    "| 3D 差值 | {matched} | {mae} | {rmsv} | {p50} | {p68} | {p90} | "
                    "{p95} | {p99} | {maxv} |"
                ).format(
                    matched=baseline["matched_epochs"],
                    mae=format_number(baseline["three_d_mae_m"]),
                    rmsv=format_number(baseline["three_d_rms_m"]),
                    p50=format_number(baseline["three_d_p50_m"]),
                    p68=format_number(baseline["three_d_p68_m"]),
                    p90=format_number(baseline["three_d_p90_m"]),
                    p95=format_number(baseline["three_d_p95_m"]),
                    p99=format_number(baseline["three_d_p99_m"]),
                    maxv=format_number(baseline["three_d_max_m"]),
                ),
            ]
        )

        anomalies = case["anomalies"]
        continuity = case["continuity"]
        parse = case["parse"]
        convergence = case["convergence"]
        lines.extend(
            [
                "",
                (
                    "收敛时刻：首次达标 "
                    f"{format_epoch(convergence['first_threshold_epoch'], nominal_start)}；首次持续收敛 "
                    f"{format_epoch(convergence['first_sustained_convergence_epoch'], nominal_start)}；"
                    f"永久收敛 {format_epoch(convergence['permanent_convergence_epoch'], nominal_start)}；"
                    f"首次 Fixed {format_epoch(convergence['first_fixed_epoch'], nominal_start)}；"
                    f"首次正确 Fixed {format_epoch(convergence['first_correct_fixed_epoch'], nominal_start)}；"
                    "首次持续正确 Fixed "
                    f"{format_epoch(convergence['first_sustained_correct_fixed_epoch'], nominal_start)}。"
                ),
                "",
                (
                    f"完整性：malformed={parse['malformed_lines']}，nonfinite={parse['nonfinite_xyz']}，"
                    f"duplicate={len(continuity['duplicate_epochs'])}，gap={len(continuity['gap_intervals'])}，"
                    f"名义网格缺失={len(continuity['nominal_grid_missing_epochs'])}，"
                    f"起点延迟={format_number(continuity['start_delay_s'], 1)} s，"
                    f"末端缺失={format_number(continuity['end_shortfall_s'], 1)} s，"
                    f"覆盖率={format_number(100.0 * continuity['coverage_fraction'], 2)}%。"
                ),
                "",
                (
                    f"异常摘要：收敛后越界区间 {len(anomalies['post_convergence_threshold_intervals'])}，"
                    f"3D>{config['three_d_anomaly_threshold_m']} m 区间 "
                    f"{len(anomalies['post_convergence_3d_anomaly_intervals'])}，"
                    f"坐标跳变 {len(anomalies['coordinate_jumps'])}，低于 {config['min_nsat']} 颗卫星区间 "
                    f"{len(anomalies['low_satellite_intervals'])}，状态切换 "
                    f"{len(anomalies['status_transitions'])}，相对 FLT 超过 "
                    f"{config['baseline_delta_threshold_m']} m 的历元 "
                    f"{len(anomalies['baseline_delta_exceedances'])}。"
                ),
                "",
                "| 收敛后异常类型 | 起点 SOW | 终点 SOW | 样本 | 最大 H | 最大 |U| | 最大 3D/差值 |",
                "|---|---:|---:|---:|---:|---:|---:|",
            ]
        )
        combined_intervals = list(anomalies["post_convergence_threshold_intervals"])
        combined_intervals.extend(anomalies["post_convergence_3d_anomaly_intervals"])
        combined_intervals.extend(anomalies["low_satellite_intervals"])
        for item in combined_intervals[:30]:
            lines.append(
                f"| {item['kind']} | {item['start_sow']:.0f} | {item['end_sow']:.0f} | "
                f"{item['samples']} | {item['max_horizontal_m']:.4f} | "
                f"{item['max_abs_vertical_m']:.4f} | {item['max_three_d_m']:.4f} |"
            )
        for item in anomalies["baseline_delta_intervals"][:30]:
            lines.append(
                f"| baseline_delta | {item['start_sow']:.0f} | {item['end_sow']:.0f} | "
                f"{item['samples']} | — | — | {item['max_value']:.4f} |"
            )
        if not combined_intervals and not anomalies["baseline_delta_intervals"]:
            lines.append("| 无 | — | — | — | — | — | — |")

        lines.extend(
            [
                "",
                "| 收敛后坐标跳变 | 到达 SOW | 跳变量 |",
                "|---:|---:|---:|",
            ]
        )
        for item in anomalies["coordinate_jumps"][:30]:
            lines.append(
                f"| {item['from_sow']:.0f} | {item['to_sow']:.0f} | {item['jump_m']:.4f} |"
            )
        if not anomalies["coordinate_jumps"]:
            lines.append("| — | — | 无 |")

        lines.extend(
            [
                "",
                "| 全段最大误差 SOW | H | U | 3D | NSat | 状态 |",
                "|---:|---:|---:|---:|---:|---|",
            ]
        )
        for item in anomalies["largest_error_epochs"]:
            lines.append(
                f"| {item['sow']:.0f} | {item['horizontal_m']:.4f} | "
                f"{item['vertical_m']:.4f} | {item['three_d_m']:.4f} | "
                f"{item['nsat']} | {item['status']} |"
            )

        lines.extend(
            [
                "",
                "| 热启动后最大误差 SOW | H | U | 3D | NSat | 状态 |",
                "|---:|---:|---:|---:|---:|---|",
            ]
        )
        for item in anomalies["largest_error_epochs_after_warmup"]:
            lines.append(
                f"| {item['sow']:.0f} | {item['horizontal_m']:.4f} | "
                f"{item['vertical_m']:.4f} | {item['three_d_m']:.4f} | "
                f"{item['nsat']} | {item['status']} |"
            )

        lines.extend(
            [
                "",
                "| 相对 FLT 最大 3D 差值 SOW | 差值 |",
                "|---:|---:|",
            ]
        )
        for item in anomalies["largest_baseline_delta_epochs"]:
            lines.append(f"| {item['sow']:.0f} | {item['delta_m']:.4f} |")
        if not anomalies["largest_baseline_delta_epochs"]:
            lines.append("| — | 无基线匹配 |")

        lines.extend(
            [
                "",
                "| 小时 | n | Fixed | H RMS | U RMS | 3D RMS | p95 | max | NSat median/min |",
                "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for hour in case["hourly"]:
            fixed_fraction = hour["fixed_fraction"]
            lines.append(
                f"| {hour['hour']:02d} | {hour['count']} | "
                f"{(100.0 * fixed_fraction if fixed_fraction is not None else 0.0):.1f}% | "
                f"{format_number(hour['horizontal_rms_m'])} | {format_number(hour['vertical_rms_m'])} | "
                f"{format_number(hour['three_d_rms_m'])} | {format_number(hour['three_d_p95_m'])} | "
                f"{format_number(hour['three_d_max_m'])} | "
                f"{format_number(hour['nsat_median'], 1)}/{format_number(hour['nsat_min'], 0)} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", required=True, help="site label used in the report")
    parser.add_argument(
        "--reference", nargs=3, type=float, required=True, metavar=("X", "Y", "Z"),
        help="reference ECEF coordinate in metres",
    )
    parser.add_argument(
        "--case", action="append", type=parse_case, required=True, metavar="LABEL=PATH",
        help="solution case; may be provided more than once",
    )
    parser.add_argument("--baseline", type=Path, help="FLT solution used for same-epoch comparison")
    parser.add_argument("--nominal-start-sow", type=float, required=True)
    parser.add_argument("--interval", type=float, default=30.0)
    parser.add_argument(
        "--hours",
        type=float,
        default=24.0,
        help="nominal analysis span in hours; fractional values are supported",
    )
    parser.add_argument("--warmup-minutes", type=float, default=30.0)
    parser.add_argument("--window-minutes", type=float, default=5.0)
    parser.add_argument("--horizontal-threshold", type=float, default=0.10)
    parser.add_argument("--vertical-threshold", type=float, default=0.20)
    parser.add_argument("--jump-threshold", type=float, default=0.10)
    parser.add_argument("--three-d-anomaly-threshold", type=float, default=0.10)
    parser.add_argument("--min-nsat", type=int, default=10)
    parser.add_argument("--baseline-delta-threshold", type=float, default=0.10)
    parser.add_argument("--markdown", type=Path, help="optional Markdown output path")
    parser.add_argument("--json", type=Path, help="optional JSON output path")
    parser.add_argument("--quiet", action="store_true", help="do not print Markdown to stdout")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.interval <= 0.0 or args.window_minutes <= 0.0:
        parser.error("interval and window-minutes must be positive")
    if args.horizontal_threshold <= 0.0 or args.vertical_threshold <= 0.0:
        parser.error("convergence thresholds must be positive")
    if args.jump_threshold <= 0.0 or args.three_d_anomaly_threshold <= 0.0:
        parser.error("anomaly thresholds must be positive")
    args.window_samples = max(1, int(math.ceil(args.window_minutes * 60.0 / args.interval)))

    baseline_rows: List[SolutionRow] = []
    if args.baseline:
        if not args.baseline.is_file():
            parser.error(f"baseline does not exist: {args.baseline}")
        baseline_rows, _ = read_solution(args.baseline)
    baseline_by_epoch = {row.sow: row for row in baseline_rows}

    cases: List[Dict[str, object]] = []
    seen_labels = set()
    for label, path in args.case:
        if label in seen_labels:
            parser.error(f"duplicate case label: {label}")
        seen_labels.add(label)
        if not path.is_file():
            parser.error(f"case does not exist: {path}")
        rows, diagnostics = read_solution(path)
        cases.append(
            analyze_case(
                label,
                path,
                rows,
                diagnostics,
                tuple(args.reference),
                args,
                baseline_by_epoch,
            )
        )

    report: Dict[str, object] = {
        "schema_version": 1,
        "site": args.site,
        "configuration": {
            "reference_ecef_m": list(args.reference),
            "nominal_start_sow": args.nominal_start_sow,
            "interval_s": args.interval,
            "hours": args.hours,
            "warmup_minutes": args.warmup_minutes,
            "window_minutes": args.window_minutes,
            "window_samples": args.window_samples,
            "horizontal_threshold_m": args.horizontal_threshold,
            "vertical_threshold_m": args.vertical_threshold,
            "jump_threshold_m": args.jump_threshold,
            "three_d_anomaly_threshold_m": args.three_d_anomaly_threshold,
            "min_nsat": args.min_nsat,
            "baseline_delta_threshold_m": args.baseline_delta_threshold,
            "baseline": str(args.baseline) if args.baseline else None,
        },
        "cases": cases,
    }
    markdown = markdown_report(report)
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown, encoding="utf-8")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(
            json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
    if not args.quiet:
        sys.stdout.write(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
