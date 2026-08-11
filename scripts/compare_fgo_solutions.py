#!/usr/bin/env python3
"""Compare two GREAT-PIFGO solution files (.fgo or .flt) for differences.

Both inputs are parsed by the same reader: their data rows share one layout
(columns 0-3 are SOW/X/Y/Z, column 13 is NSat, column 16 is AmbStatus; header
comments start with '#'). Statistics reuse ``analyze_fgo_solution`` so the
definitions stay identical to the per-file analyzer.

Reports:
  * per-file overview, continuity, and optional nominal-grid coverage;
  * common-epoch coordinate differences: 3D distance (always), E/N/U components
    and horizontal/vertical (when a reference is supplied), largest differences,
    and intervals above ``--delta-threshold``;
  * Fixed/Float status consistency cross-table;
  * each file's own accuracy against the reference truth (when supplied).

Example:

  python scripts/compare_fgo_solutions.py --site GODN \
    --a PARAMETER=result/GODN-CODEX_RAWPAR2_UPD_PARAMETER.fgo \
    --b FLT=result/GODN-CODEX_RAWPAR2_UPD_PARAMETER.flt \
    --reference 1130760.6931 -4831298.6759 3994155.1990 \
    --nominal-start-sow 259200 --interval 30 --hours 24 \
    --markdown result/GODN-PARAMETER-vs-FLT.md --json result/GODN-PARAMETER-vs-FLT.json
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from analyze_fgo_solution import (
    SolutionRow,
    accuracy_metrics,
    build_errors,
    component_statistics,
    continuity_metrics,
    ecef_to_geodetic_lat_lon,
    format_epoch,
    format_fraction,
    format_number,
    normalize_status,
    parse_case,
    percentile,
    read_solution,
    rms,
)


def ecef_vector_to_enu(
    dx: float, dy: float, dz: float, lat: float, lon: float
) -> Tuple[float, float, float]:
    """Rotate an ECEF difference vector into the local E/N/U frame."""
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    east = -sin_lon * dx + cos_lon * dy
    north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    up = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
    return east, north, up


def magnitude_statistics(values: Sequence[float]) -> Dict[str, Optional[float]]:
    """Statistics for a non-negative magnitude series (e.g. 3D distance)."""
    return {
        "count": len(values),
        "mean_m": statistics.fmean(values) if values else None,
        "mae_m": statistics.fmean(values) if values else None,
        "rms_m": rms(values),
        "p50_m": percentile(values, 0.50),
        "p68_m": percentile(values, 0.68),
        "p90_m": percentile(values, 0.90),
        "p95_m": percentile(values, 0.95),
        "p99_m": percentile(values, 0.99),
        "max_m": max(values) if values else None,
    }


def group_intervals(
    epochs_values: Sequence[Tuple[float, float]], interval: float
) -> List[Dict[str, float]]:
    """Group (sow, value) pairs into consecutive-sample runs."""
    ordered = sorted(epochs_values, key=lambda item: item[0])
    groups: List[List[Tuple[float, float]]] = []
    for sow, value in ordered:
        if groups and abs((sow - groups[-1][-1][0]) - interval) <= 1e-6:
            groups[-1].append((sow, value))
        else:
            groups.append([(sow, value)])
    return [
        {
            "start_sow": group[0][0],
            "end_sow": group[-1][0],
            "samples": len(group),
            "max_m": max(value for _, value in group),
        }
        for group in groups
    ]


def is_fixed(status: str) -> bool:
    """True when an ambiguity status string denotes a fixed solution."""
    return normalize_status(status) == "FIXED"


def accuracy_pair(
    rows: Sequence[SolutionRow],
    reference: Tuple[float, float, float],
    args: argparse.Namespace,
) -> Tuple[Dict[str, object], Dict[str, object]]:
    """Per-file accuracy (all span / after warmup) against the reference."""
    errors = build_errors(rows, reference, args.horizontal_threshold, args.vertical_threshold)
    all_metrics = accuracy_metrics(errors)
    nominal_start = args.nominal_start_sow if args.nominal_start_sow is not None else (
        rows[0].sow if rows else None
    )
    warm_start = (
        nominal_start + args.warmup_minutes * 60.0 if nominal_start is not None else None
    )
    warm_metrics = accuracy_metrics(
        [row for row in errors if warm_start is not None and row.sow >= warm_start]
    )
    return all_metrics, warm_metrics


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for the difference-evaluation script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", required=True, help="site label used in the report")
    parser.add_argument(
        "--a", type=parse_case, required=True, metavar="LABEL=PATH",
        help="first solution file (.fgo or .flt)",
    )
    parser.add_argument(
        "--b", type=parse_case, required=True, metavar="LABEL=PATH",
        help="second solution file (.fgo or .flt)",
    )
    parser.add_argument(
        "--reference", nargs=3, type=float, metavar=("X", "Y", "Z"),
        help="reference ECEF coordinate; enables E/N/U decomposition and per-file accuracy",
    )
    parser.add_argument("--nominal-start-sow", type=float)
    parser.add_argument("--interval", type=float, default=30.0)
    parser.add_argument(
        "--hours", type=float, default=24.0,
        help="nominal analysis span in hours (for grid coverage)",
    )
    parser.add_argument("--warmup-minutes", type=float, default=30.0)
    parser.add_argument("--delta-threshold", type=float, default=0.10)
    parser.add_argument("--horizontal-threshold", type=float, default=0.10)
    parser.add_argument("--vertical-threshold", type=float, default=0.20)
    parser.add_argument("--markdown", type=Path, help="optional Markdown output path")
    parser.add_argument("--json", type=Path, help="optional JSON output path")
    parser.add_argument("--quiet", action="store_true", help="do not print Markdown to stdout")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point: parse args, compare the two files, emit report."""
    parser = build_parser()
    args = parser.parse_args(argv)
    finite_values = (
        args.interval,
        args.hours,
        args.warmup_minutes,
        args.delta_threshold,
        args.horizontal_threshold,
        args.vertical_threshold,
        *(args.reference or ()),
    )
    if args.nominal_start_sow is not None:
        finite_values += (args.nominal_start_sow,)
    if not all(math.isfinite(value) for value in finite_values):
        parser.error("reference and comparison parameters must be finite")
    if args.reference and math.sqrt(sum(value * value for value in args.reference)) <= 1e-9:
        parser.error("reference ECEF coordinate must be non-zero")
    if args.interval <= 0.0 or args.horizontal_threshold <= 0.0 or args.vertical_threshold <= 0.0:
        parser.error("interval and thresholds must be positive")
    if args.hours <= 0.0:
        parser.error("hours must be positive")
    if args.delta_threshold <= 0.0:
        parser.error("delta-threshold must be positive")
    if args.warmup_minutes < 0.0:
        parser.error("warmup-minutes must be non-negative")
    if args.a[0] == args.b[0]:
        parser.error("--a and --b must use different labels")
    for label, path in (args.a, args.b):
        if not path.is_file():
            parser.error(f"case does not exist: {label}={path}")

    label_a, path_a = args.a
    label_b, path_b = args.b
    rows_a, diag_a = read_solution(path_a)
    rows_b, diag_b = read_solution(path_b)
    if not rows_a or not rows_b:
        parser.error("both files must contain at least one data row")

    a_by_sow = {row.sow: row for row in rows_a}
    b_by_sow = {row.sow: row for row in rows_b}
    common = sorted(set(a_by_sow) & set(b_by_sow))
    only_a = sorted(set(a_by_sow) - set(b_by_sow))
    only_b = sorted(set(b_by_sow) - set(a_by_sow))
    # Match the analyzer's baseline definition: skip epochs where either side
    # has a non-finite coordinate (baseline_metrics filters the same way).
    aligned = [
        (sow, a_by_sow[sow], b_by_sow[sow])
        for sow in common
        if all(math.isfinite(v) for v in a_by_sow[sow].xyz + b_by_sow[sow].xyz)
    ]

    lat = lon = None
    if args.reference:
        lat, lon = ecef_to_geodetic_lat_lon(*args.reference)

    three_d = [math.dist(a.xyz, b.xyz) for _, a, b in aligned]
    components: Optional[Dict[str, object]] = None
    enu_diffs: List[Tuple[float, float, float]] = []
    if lat is not None:
        for _, a, b in aligned:
            enu_diffs.append(
                ecef_vector_to_enu(a.x - b.x, a.y - b.y, a.z - b.z, lat, lon)
            )
        east = [item[0] for item in enu_diffs]
        north = [item[1] for item in enu_diffs]
        up = [item[2] for item in enu_diffs]
        horizontal = [math.hypot(e, n) for e, n, _ in enu_diffs]
        components = {
            "east": component_statistics(east, signed_component=True),
            "north": component_statistics(north, signed_component=True),
            "up": component_statistics(up, signed_component=True),
            "horizontal": component_statistics(horizontal, signed_component=False),
            "three_d": magnitude_statistics(three_d),
        }

    differences: Dict[str, object] = {"three_d": magnitude_statistics(three_d)}
    if components is not None:
        differences["components"] = components

    differences["largest"] = [
        {
            "sow": aligned[i][0],
            "three_d_m": three_d[i],
            "east_m": enu_diffs[i][0] if enu_diffs else None,
            "north_m": enu_diffs[i][1] if enu_diffs else None,
            "up_m": enu_diffs[i][2] if enu_diffs else None,
        }
        for _, i in sorted(
            ((three_d[i], i) for i in range(len(aligned))),
            key=lambda item: item[0],
            reverse=True,
        )[:10]
    ]
    over_threshold = [
        (sow, three_d[i])
        for i, (sow, _, _) in enumerate(aligned)
        if three_d[i] > args.delta_threshold
    ]
    differences["over_threshold_intervals"] = group_intervals(over_threshold, args.interval)
    differences["over_threshold_epochs"] = len(over_threshold)

    cross = {"both_fixed": 0, "both_float": 0, "a_fixed_b_float": 0, "a_float_b_fixed": 0}
    mismatches: List[Dict[str, str]] = []
    for sow, a, b in aligned:
        a_fixed = is_fixed(a.status)
        b_fixed = is_fixed(b.status)
        if a_fixed and b_fixed:
            cross["both_fixed"] += 1
        elif not a_fixed and not b_fixed:
            cross["both_float"] += 1
        elif a_fixed and not b_fixed:
            cross["a_fixed_b_float"] += 1
        else:
            cross["a_float_b_fixed"] += 1
        if a_fixed != b_fixed:
            mismatches.append({"sow": sow, "a_status": a.status, "b_status": b.status})
    status = {
        "cross_table": cross,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches[:100],
    }

    accuracy = None
    if args.reference:
        reference = tuple(args.reference)
        accuracy = {
            "a": {"label": label_a, "all": None, "after_warmup": None},
            "b": {"label": label_b, "all": None, "after_warmup": None},
        }
        for key, rows in (("a", rows_a), ("b", rows_b)):
            all_metrics, warm_metrics = accuracy_pair(rows, reference, args)
            accuracy[key]["all"] = {
                "three_d_rms_m": all_metrics["three_d_rms_m"],
                "three_d_mae_m": all_metrics["three_d_mae_m"],
                "three_d_p95_m": all_metrics["three_d_p95_m"],
                "three_d_p99_m": all_metrics["three_d_p99_m"],
                "three_d_max_m": all_metrics["three_d_max_m"],
                "fixed_fraction": all_metrics["fixed_fraction"],
                "correct_fraction": all_metrics["correct_fraction"],
            }
            accuracy[key]["after_warmup"] = {
                "three_d_rms_m": warm_metrics["three_d_rms_m"],
                "three_d_mae_m": warm_metrics["three_d_mae_m"],
                "three_d_p95_m": warm_metrics["three_d_p95_m"],
                "three_d_p99_m": warm_metrics["three_d_p99_m"],
                "three_d_max_m": warm_metrics["three_d_max_m"],
                "fixed_fraction": warm_metrics["fixed_fraction"],
                "correct_fraction": warm_metrics["correct_fraction"],
            }

    nominal_start = args.nominal_start_sow
    grid = {
        "common_epochs": len(common),
        "only_a": {"count": len(only_a), "epochs": only_a[:100]},
        "only_b": {"count": len(only_b), "epochs": only_b[:100]},
    }
    for key, rows in (("a", rows_a), ("b", rows_b)):
        grid[key] = {
            "rows": len(rows),
            "first_sow": rows[0].sow,
            "last_sow": rows[-1].sow,
            "parse": diag_a if key == "a" else diag_b,
        }
        if nominal_start is not None:
            grid[key]["continuity"] = continuity_metrics(
                rows, args.interval, nominal_start, args.hours
            )

    report: Dict[str, object] = {
        "schema_version": 1,
        "site": args.site,
        "configuration": {
            "label_a": label_a,
            "path_a": str(path_a),
            "label_b": label_b,
            "path_b": str(path_b),
            "reference_ecef_m": list(args.reference) if args.reference else None,
            "nominal_start_sow": nominal_start,
            "interval_s": args.interval,
            "hours": args.hours,
            "warmup_minutes": args.warmup_minutes,
            "delta_threshold_m": args.delta_threshold,
            "horizontal_threshold_m": args.horizontal_threshold,
            "vertical_threshold_m": args.vertical_threshold,
        },
        "grid": grid,
        "differences": differences,
        "status": status,
        "accuracy": accuracy,
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


def show_epoch(value: float, nominal_start: Optional[float]) -> str:
    """Format an epoch as SOW (+offset from nominal start, when available)."""
    if nominal_start is not None:
        return format_epoch(value, nominal_start)
    return f"{value:.0f}"


def markdown_report(report: Dict[str, object]) -> str:
    """Render the comparison report as Chinese Markdown (analyzer style)."""
    config = report["configuration"]
    nominal_start = config["nominal_start_sow"]
    label_a, label_b = config["label_a"], config["label_b"]
    interval = float(config["interval_s"])
    intro_parts = [f"采样间隔 {interval} s"]
    if config["reference_ecef_m"]:
        intro_parts.insert(0, f"参考坐标 ECEF: `{config['reference_ecef_m']}`")
    intro_parts.append(f"差异阈值 {format_number(config['delta_threshold_m'])} m")
    lines = [
        f"# {report['site']} 结果差异：{label_a} vs {label_b}",
        "",
        "；".join(intro_parts) + "。",
        "",
        "## 输入概览",
        "",
        "| 文件 | 行数 | SOW 范围 | malformed | nonfinite SOW | nonfinite XYZ |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key in ("a", "b"):
        grid = report["grid"][key]
        lines.append(
            f"| {config['label_' + key]} | {grid['rows']} | "
            f"{grid['first_sow']:.0f}--{grid['last_sow']:.0f} | "
            f"{grid['parse']['malformed_lines']} | {grid['parse']['nonfinite_sow']} | "
            f"{grid['parse']['nonfinite_xyz']} |"
        )

    lines.extend(["", "## 时间网格对比", ""])
    grid = report["grid"]
    lines.append(
        f"共同历元 **{grid['common_epochs']}**；"
        f"仅 {label_a} 有 {grid['only_a']['count']} 个"
        + (
            f"（前 {len(grid['only_a']['epochs'])} 个 SOW: "
            + ", ".join(show_epoch(s, nominal_start) for s in grid['only_a']['epochs'])
            + "）" if grid["only_a"]["epochs"] else ""
        )
        + f"；仅 {label_b} 有 {grid['only_b']['count']} 个"
        + (
            f"（前 {len(grid['only_b']['epochs'])} 个 SOW: "
            + ", ".join(show_epoch(s, nominal_start) for s in grid['only_b']['epochs'])
            + "）" if grid["only_b"]["epochs"] else ""
        )
        + "。"
    )
    for key in ("a", "b"):
        continuity = grid[key].get("continuity")
        if not continuity:
            continue
        lines.append(
            f"- {config['label_' + key]}：覆盖率 "
            f"{format_fraction(continuity['coverage_fraction'])}，"
            f"重复 {len(continuity['duplicate_epochs'])}，间隔 "
            f"{len(continuity['gap_intervals'])}，名义网格缺失 "
            f"{len(continuity['nominal_grid_missing_epochs'])}，"
            f"完整={continuity['complete_nominal_span']}。"
        )

    lines.extend(["", "## 坐标差异（共同历元）", ""])
    differences = report["differences"]
    three_d = differences["three_d"]
    lines.extend(
        [
            "| 指标 | 3D 距离 (m) |",
            "|---|---:|",
            f"| count | {three_d['count']} |",
            f"| mean | {format_number(three_d['mean_m'])} |",
            f"| RMS | {format_number(three_d['rms_m'])} |",
            f"| P50 | {format_number(three_d['p50_m'])} |",
            f"| P68 | {format_number(three_d['p68_m'])} |",
            f"| P90 | {format_number(three_d['p90_m'])} |",
            f"| P95 | {format_number(three_d['p95_m'])} |",
            f"| P99 | {format_number(three_d['p99_m'])} |",
            f"| max | {format_number(three_d['max_m'])} |",
        ]
    )

    components = differences.get("components")
    if components:
        lines.extend(
            [
                "",
                "E/N/U 分量差（A−B，单位 m；bias/std 保留符号，分位数按绝对分量误差）：",
                "",
                "| 分量 | bias | std | RMSE | MAE | P95 | P99 | max |",
                "|---|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for name, key in (("E", "east"), ("N", "north"), ("U", "up"), ("H", "horizontal")):
            stats = components[key]
            lines.append(
                f"| {name} | {format_number(stats['bias_m'])} | {format_number(stats['std_m'])} | "
                f"{format_number(stats['rmse_m'])} | {format_number(stats['mae_m'])} | "
                f"{format_number(stats['p95_abs_m'])} | {format_number(stats['p99_abs_m'])} | "
                f"{format_number(stats['max_abs_m'])} |"
            )

    lines.extend(
        [
            "",
            f"超过 {format_number(config['delta_threshold_m'])} m 的历元："
            f"**{differences['over_threshold_epochs']}** 个"
            + (
                f"，区间 {len(differences['over_threshold_intervals'])} 段"
                if differences["over_threshold_intervals"]
                else ""
            )
            + "。",
            "",
        ]
    )
    if differences["over_threshold_intervals"]:
        lines.extend(
            [
                "| 超差区间 | 起点 SOW | 终点 SOW | 样本 | 最大差值 (m) |",
                "|---|---:|---:|---:|---:|",
            ]
        )
        for item in differences["over_threshold_intervals"][:30]:
            lines.append(
                f"| — | {show_epoch(item['start_sow'], nominal_start)} | "
                f"{show_epoch(item['end_sow'], nominal_start)} | "
                f"{item['samples']} | {format_number(item['max_m'])} |"
            )
        lines.append("")

    lines.extend(
        [
            "| 最大 3D 差 SOW | 3D | E | N | U |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for item in differences["largest"]:
        east = format_number(item["east_m"]) if item["east_m"] is not None else "—"
        north = format_number(item["north_m"]) if item["north_m"] is not None else "—"
        up = format_number(item["up_m"]) if item["up_m"] is not None else "—"
        lines.append(
            f"| {show_epoch(item['sow'], nominal_start)} | "
            f"{format_number(item['three_d_m'])} | {east} | {north} | {up} |"
        )

    lines.extend(["", "## 状态一致性（共同历元）", ""])
    status = report["status"]
    cross = status["cross_table"]
    lines.extend(
        [
            "| 组合 | 历元数 |",
            "|---|---:|",
            f"| {label_a} Fixed / {label_b} Fixed | {cross['both_fixed']} |",
            f"| {label_a} Fixed / {label_b} Float | {cross['a_fixed_b_float']} |",
            f"| {label_a} Float / {label_b} Fixed | {cross['a_float_b_fixed']} |",
            f"| {label_a} Float / {label_b} Float | {cross['both_float']} |",
        ]
    )
    if status["mismatches"]:
        lines.append(
            "",
        )
        lines.append(
            f"状态不一致 {status['mismatch_count']} 个历元（前 "
            f"{len(status['mismatches'])} 个）："
        )
        lines.append("")
        lines.extend(
            f"- SOW {show_epoch(item['sow'], nominal_start)}: {label_a}={item['a_status']} / "
            f"{label_b}={item['b_status']}"
            for item in status["mismatches"][:20]
        )

    accuracy = report.get("accuracy")
    if accuracy:
        lines.extend(["", "## 各自相对参考真值精度", ""])
        lines.extend(
            [
                "| 文件 | 时段 | 3D RMS | MAE | P95 | P99 | max | Fixed | 正确率 |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for key in ("a", "b"):
            entry = accuracy[key]
            for period, field in (("全段", "all"), ("热启动后", "after_warmup")):
                stats = entry[field]
                fixed = stats["fixed_fraction"]
                fixed_text = (
                    f"{100.0 * fixed:.2f}%" if fixed is not None else "—"
                )
                correct = stats["correct_fraction"]
                correct_text = (
                    f"{100.0 * correct:.2f}%" if correct is not None else "—"
                )
                lines.append(
                    f"| {config['label_' + key]} | {period} | "
                    f"{format_number(stats['three_d_rms_m'])} | "
                    f"{format_number(stats['three_d_mae_m'])} | "
                    f"{format_number(stats['three_d_p95_m'])} | "
                    f"{format_number(stats['three_d_p99_m'])} | "
                    f"{format_number(stats['three_d_max_m'])} | "
                    f"{fixed_text} | {correct_text} |"
                )

    lines.extend(
        [
            "",
            "## Artifacts",
            "",
            f"- 机器可读结果: `{config['path_a']}` vs `{config['path_b']}`",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
