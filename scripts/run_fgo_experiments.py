#!/usr/bin/env python3
"""Launch, monitor, validate, and summarize GREAT-PIFGO experiments.

The runner accepts either one or more ``--case NAME=CONFIG`` arguments or a
JSON manifest.  It deliberately treats a zero process exit code as necessary
but insufficient: every expected FGO file must be fresh, complete, continuous,
finite, and free from configured fatal log patterns.

Examples
--------

Run one existing configuration (the data directory is inferred from its bias
product):

  python scripts/run_fgo_experiments.py run \
    --case UPD_FF_NONE=build/parallel30/UPD_NONE.xml \
    --reference GODN=1130760.6931,-4831298.6759,3994155.1990

Run a reusable matrix:

  python scripts/run_fgo_experiments.py run \
    --manifest scripts/fgo_experiment_matrix.example.json

Inspect an interrupted or concurrently running test from another terminal:

  python scripts/run_fgo_experiments.py status --run-dir build/fgo_runs/<run-id>
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import socket
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parent.parent
DEFAULT_EXECUTABLE = REPO_ROOT / "build" / "Bin" / "Release" / "GREAT_PVT.exe"
DEFAULT_ANALYZER = SCRIPT_PATH.with_name("analyze_fgo_solution.py")
GPS_EPOCH = dt.datetime(1980, 1, 6)

# These patterns describe conditions for which a numerically produced file is
# not trustworthy.  A manifest can replace them with ``failure_patterns``.
DEFAULT_FAILURE_PATTERNS = (
    r"segmentation fault",
    r"access violation",
    r"terminate called",
    r"ceres[^\r\n]*(?:failure|failed|error)",
    r"covariance[^\r\n]*(?:failure|failed|rank deficient)",
    r"raw[^\r\n]*prepar[^\r\n]*failed",
    r"\bnon[- ]?finite\b",
    r"\bnan\b",
    r"\binf\b",
    r"capacity[^\r\n]*(?:exceed|overflow)",
    r"rank deficient",
    r"pseudo(?:inverse|-inverse)",
)


class ConfigurationError(ValueError):
    """Raised for an invalid experiment definition."""


@dataclass(frozen=True)
class OutputSpec:
    station: str
    fgo: Path
    ppp: Optional[Path]


@dataclass
class Experiment:
    name: str
    config: Path
    workdir: Path
    executable: Path
    stations: List[str]
    begin: dt.datetime
    end: dt.datetime
    interval_s: float
    mode: str
    frequency: Optional[int]
    solver_threads: int
    product: str
    outputs: List[OutputSpec]
    diagnostic_logs: List[Path]
    arguments: List[str] = field(default_factory=list)
    tags: Dict[str, str] = field(default_factory=dict)

    @property
    def expected_rows(self) -> int:
        # Configuration begin/end are inclusive processing epochs.
        return max(
            0,
            int(math.floor((self.end - self.begin).total_seconds() / self.interval_s + 1e-9)) + 1,
        )

    @property
    def nominal_hours(self) -> float:
        return ((self.end - self.begin).total_seconds() + self.interval_s) / 3600.0


@dataclass
class ProcessEntry:
    experiment: Experiment
    process: subprocess.Popen
    stdout_stream: object
    stderr_stream: object
    stdout_path: Path
    stderr_path: Path
    started_wall: float
    started_iso: str
    output_before: Dict[str, Optional[Tuple[int, int]]]
    timed_out: bool = False
    peak_rss_bytes: int = 0
    last_rss_bytes: int = 0


def resolve_path(value: str | Path, base: Path = REPO_ROOT) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def parse_datetime(value: str, field_name: str) -> dt.datetime:
    value = value.strip().replace("[GPS]", "")
    try:
        return dt.datetime.fromisoformat(value)
    except ValueError as exc:
        raise ConfigurationError(f"invalid {field_name} time: {value}") from exc


def xml_text(root: ET.Element, path: str, default: Optional[str] = None) -> Optional[str]:
    node = root.find(path)
    if node is None or node.text is None or not node.text.strip():
        return default
    return node.text.strip()


def substitute_station(path_text: str, station: str) -> str:
    return path_text.replace("$(rec)", station).replace("${rec}", station)


def classify_product(root: ET.Element) -> str:
    bias = (xml_text(root, "./inputs/bias", "") or "").upper()
    upd_mode = (xml_text(root, "./ambiguity/upd_mode", "") or "").upper()
    if "OSB" in bias or "WUM0MGX" in bias or upd_mode == "OSB":
        return "OSB"
    return "UPD"


def infer_workdir(root: ET.Element) -> Path:
    product = classify_product(root)
    directory = "PPPFLT_2023305_OSB" if product == "OSB" else "PPPFLT_2023305"
    return REPO_ROOT / "sample_data" / directory


def parse_output_specs(
    root: ET.Element,
    workdir: Path,
    stations: Sequence[str],
    explicit: Optional[Mapping[str, object]] = None,
) -> List[OutputSpec]:
    fgo_template = xml_text(root, "./outputs/fgo")
    ppp_template = xml_text(root, "./outputs/ppp")
    outputs: List[OutputSpec] = []
    for station in stations:
        override = explicit.get(station) if explicit else None
        if isinstance(override, str):
            fgo_path = resolve_path(override, workdir)
            ppp_path = None
        elif isinstance(override, Mapping):
            if "fgo" not in override:
                raise ConfigurationError(f"explicit output for {station} has no fgo path")
            fgo_path = resolve_path(str(override["fgo"]), workdir)
            ppp_path = (
                resolve_path(str(override["ppp"]), workdir)
                if override.get("ppp")
                else None
            )
        else:
            if not fgo_template:
                raise ConfigurationError("configuration has no <outputs><fgo> path")
            fgo_path = resolve_path(substitute_station(fgo_template, station), workdir)
            ppp_path = (
                resolve_path(substitute_station(ppp_template, station), workdir)
                if ppp_template
                else None
            )
        outputs.append(OutputSpec(station=station, fgo=fgo_path, ppp=ppp_path))
    return outputs


def load_experiment(
    raw: Mapping[str, object],
    defaults: Mapping[str, object],
    manifest_base: Path,
) -> Experiment:
    if not raw.get("name") or not raw.get("config"):
        raise ConfigurationError("each experiment requires name and config")
    name = str(raw["name"]).strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
        raise ConfigurationError(
            f"invalid experiment name {name!r}; use letters, digits, dot, dash, or underscore"
        )
    config = resolve_path(str(raw["config"]), manifest_base)
    if not config.is_file():
        raise ConfigurationError(f"configuration does not exist: {config}")
    try:
        root = ET.parse(config).getroot()
    except ET.ParseError as exc:
        raise ConfigurationError(f"invalid XML {config}: {exc}") from exc

    configured_workdir = raw.get("workdir", defaults.get("workdir"))
    workdir = (
        resolve_path(str(configured_workdir), manifest_base)
        if configured_workdir
        else infer_workdir(root).resolve()
    )
    executable_value = raw.get("executable", defaults.get("executable", DEFAULT_EXECUTABLE))
    executable = resolve_path(str(executable_value), manifest_base)
    if not workdir.is_dir():
        raise ConfigurationError(f"working directory does not exist: {workdir}")
    if not executable.is_file():
        raise ConfigurationError(f"executable does not exist: {executable}")

    begin = parse_datetime(xml_text(root, "./gen/beg", "") or "", "begin")
    end = parse_datetime(xml_text(root, "./gen/end", "") or "", "end")
    interval_s = float(xml_text(root, "./gen/int", "30") or "30")
    if interval_s <= 0 or end <= begin:
        raise ConfigurationError(f"invalid time span or interval in {config}")
    station_value = raw.get("stations") or (xml_text(root, "./gen/rec", "") or "").split()
    if isinstance(station_value, str):
        station_value = station_value.split()
    stations = [str(item).strip().upper() for item in station_value if str(item).strip()]
    if not stations:
        raise ConfigurationError(f"no stations configured in {config}")

    outputs = parse_output_specs(root, workdir, stations, raw.get("outputs"))
    log_node = root.find("./outputs/log")
    diagnostic_logs: List[Path] = []
    if log_node is not None and log_node.get("name"):
        log_name = log_node.get("name", "")
        if "$(rec)" in log_name or "${rec}" in log_name:
            diagnostic_logs.extend(
                resolve_path(substitute_station(log_name, station), workdir)
                for station in stations
            )
        else:
            diagnostic_logs.append(resolve_path(log_name, workdir))
    mode = (xml_text(root, "./fgo/ambiguity_feedback_mode", "NONE") or "NONE").upper()
    frequency_text = xml_text(root, "./process/frequency")
    tags = {str(key): str(value) for key, value in dict(raw.get("tags", {})).items()}
    return Experiment(
        name=name,
        config=config,
        workdir=workdir,
        executable=executable,
        stations=stations,
        begin=begin,
        end=end,
        interval_s=interval_s,
        mode=mode,
        frequency=int(frequency_text) if frequency_text else None,
        solver_threads=int(xml_text(root, "./fgo/gnss_num_threads", "1") or "1"),
        product=classify_product(root),
        outputs=outputs,
        diagnostic_logs=diagnostic_logs,
        arguments=[str(item) for item in raw.get("arguments", [])],
        tags=tags,
    )


def parse_reference(value: str) -> Tuple[str, Tuple[float, float, float]]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("reference must be SITE=X,Y,Z")
    site, coordinates = value.split("=", 1)
    try:
        xyz = tuple(float(item.strip()) for item in coordinates.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("reference coordinates must be numbers") from exc
    if len(xyz) != 3 or not all(math.isfinite(item) for item in xyz):
        raise argparse.ArgumentTypeError("reference must contain three finite coordinates")
    return site.strip().upper(), xyz  # type: ignore[return-value]


def parse_named_path(value: str, kind: str) -> Tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"{kind} must be NAME=PATH")
    name, path = value.split("=", 1)
    if not name.strip() or not path.strip():
        raise argparse.ArgumentTypeError(f"{kind} must contain a name and path")
    return name.strip(), path.strip()


def read_manifest(args: argparse.Namespace) -> Tuple[dict, List[Experiment], Dict[str, Tuple[float, float, float]]]:
    if args.manifest:
        manifest_path = resolve_path(args.manifest)
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ConfigurationError(f"cannot read manifest {manifest_path}: {exc}") from exc
        if manifest.get("schema_version", 1) != 1:
            raise ConfigurationError("unsupported manifest schema_version")
        manifest_base = manifest_path.parent
        raw_experiments = manifest.get("experiments", [])
        defaults = manifest.get("defaults", {})
    else:
        if not args.case:
            raise ConfigurationError("run requires --manifest or at least one --case")
        manifest_path = None
        manifest_base = REPO_ROOT
        workdirs = dict(args.work_dir or [])
        raw_experiments = []
        for name, config in args.case:
            raw: Dict[str, object] = {"name": name, "config": config}
            if name in workdirs:
                raw["workdir"] = workdirs[name]
            elif args.default_work_dir:
                raw["workdir"] = args.default_work_dir
            raw_experiments.append(raw)
        defaults = {"executable": str(args.executable or DEFAULT_EXECUTABLE)}
        manifest = {"schema_version": 1, "experiments": raw_experiments, "defaults": defaults}

    if not isinstance(raw_experiments, list) or not raw_experiments:
        raise ConfigurationError("manifest has no experiments")
    if not isinstance(defaults, Mapping):
        raise ConfigurationError("manifest defaults must be an object")
    experiments = [load_experiment(item, defaults, manifest_base) for item in raw_experiments]
    names = [experiment.name for experiment in experiments]
    if len(names) != len(set(names)):
        raise ConfigurationError("experiment names must be unique")
    owners: Dict[Path, str] = {}
    for experiment in experiments:
        critical_argument = next(
            (
                value for value in experiment.arguments
                if re.search(r"config:(?:gen:(?:beg|end|int|rec)|outputs)", value, re.IGNORECASE)
            ),
            None,
        )
        if critical_argument:
            raise ConfigurationError(
                f"{experiment.name}: command-line override changes monitored XML metadata: "
                f"{critical_argument}; use a dedicated XML instead"
            )
        paths = [output.fgo for output in experiment.outputs]
        paths.extend(output.ppp for output in experiment.outputs if output.ppp)
        paths.extend(experiment.diagnostic_logs)
        for path in paths:
            if path in owners:
                raise ConfigurationError(
                    f"experiments {owners[path]} and {experiment.name} share output {path}"
                )
            owners[path] = experiment.name

    references: Dict[str, Tuple[float, float, float]] = {}
    for site, coordinates in dict(manifest.get("references", {})).items():
        if not isinstance(coordinates, list) or len(coordinates) != 3:
            raise ConfigurationError(f"reference for {site} must contain X, Y, Z")
        xyz = tuple(float(value) for value in coordinates)
        if not all(math.isfinite(value) for value in xyz):
            raise ConfigurationError(f"reference for {site} is not finite")
        references[str(site).upper()] = xyz  # type: ignore[assignment]
    for site, xyz in args.reference or []:
        references[site] = xyz
    return manifest, experiments, references


def file_signature(path: Path) -> Optional[Tuple[int, int]]:
    try:
        stat = path.stat()
    except OSError:
        return None
    return stat.st_mtime_ns, stat.st_size


def solution_progress(path: Path) -> Dict[str, object]:
    rows = 0
    malformed = 0
    nonfinite = 0
    epochs: List[float] = []
    if not path.is_file():
        return {"exists": False, "rows": 0, "malformed": 0, "nonfinite": 0, "last_sow": None}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            for line in stream:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.split()
                try:
                    sow = float(fields[0])
                    xyz = (float(fields[1]), float(fields[2]), float(fields[3]))
                except (IndexError, ValueError):
                    malformed += 1
                    continue
                rows += 1
                epochs.append(sow)
                if not all(math.isfinite(value) for value in xyz):
                    nonfinite += 1
    except OSError:
        return {"exists": False, "rows": 0, "malformed": 0, "nonfinite": 0, "last_sow": None}
    return {
        "exists": True,
        "rows": rows,
        "malformed": malformed,
        "nonfinite": nonfinite,
        "first_sow": epochs[0] if epochs else None,
        "last_sow": epochs[-1] if epochs else None,
        "epochs": epochs,
    }


def process_resource_snapshot(pid: int) -> Dict[str, int]:
    """Return portable best-effort RSS counters without a third-party package."""
    if pid <= 0:
        return {}
    if os.name == "nt":
        try:
            import ctypes
            from ctypes import wintypes

            class Counters(ctypes.Structure):
                _fields_ = [
                    ("cb", wintypes.DWORD),
                    ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t),
                ]

            kernel32 = ctypes.windll.kernel32
            psapi = ctypes.windll.psapi
            handle = kernel32.OpenProcess(0x1000 | 0x0400, False, pid)
            if not handle:
                return {}
            counters = Counters()
            counters.cb = ctypes.sizeof(Counters)
            try:
                if not psapi.GetProcessMemoryInfo(
                    handle, ctypes.byref(counters), counters.cb
                ):
                    return {}
                return {
                    "rss_bytes": int(counters.WorkingSetSize),
                    "peak_rss_bytes": int(counters.PeakWorkingSetSize),
                }
            finally:
                kernel32.CloseHandle(handle)
        except (AttributeError, OSError, ValueError):
            return {}
    status_path = Path(f"/proc/{pid}/status")
    try:
        values = {}
        for line in status_path.read_text(encoding="utf-8").splitlines():
            if line.startswith("VmRSS:") or line.startswith("VmHWM:"):
                fields = line.split()
                if len(fields) >= 2:
                    values[line.split(":", 1)[0]] = int(fields[1]) * 1024
        return {
            "rss_bytes": values.get("VmRSS", 0),
            "peak_rss_bytes": values.get("VmHWM", values.get("VmRSS", 0)),
        }
    except (OSError, ValueError):
        return {}


def sample_process_resources(entry: ProcessEntry) -> None:
    sample = process_resource_snapshot(entry.process.pid)
    if not sample:
        return
    entry.last_rss_bytes = sample.get("rss_bytes", 0)
    entry.peak_rss_bytes = max(
        entry.peak_rss_bytes,
        sample.get("peak_rss_bytes", 0),
        sample.get("rss_bytes", 0),
    )


def entry_output_progress(entry: ProcessEntry, output: OutputSpec) -> Dict[str, object]:
    """Hide a pre-existing result until this process has changed the file."""
    if file_signature(output.fgo) == entry.output_before.get(str(output.fgo)):
        return {"exists": False, "rows": 0, "malformed": 0, "nonfinite": 0, "last_sow": None}
    return solution_progress(output.fgo)


def validate_output(
    experiment: Experiment,
    output: OutputSpec,
    before: Optional[Tuple[int, int]],
) -> Tuple[Dict[str, object], List[str]]:
    progress = solution_progress(output.fgo)
    result = {key: value for key, value in progress.items() if key != "epochs"}
    result["path"] = str(output.fgo)
    result["expected_rows"] = experiment.expected_rows
    result["fresh"] = file_signature(output.fgo) != before
    errors: List[str] = []
    if not result["exists"]:
        errors.append(f"{output.station}: missing FGO output")
        return result, errors
    if not result["fresh"]:
        errors.append(f"{output.station}: output was not updated by this run")
    if result["rows"] != experiment.expected_rows:
        errors.append(
            f"{output.station}: expected {experiment.expected_rows} rows, got {result['rows']}"
        )
    if result["malformed"]:
        errors.append(f"{output.station}: {result['malformed']} malformed row(s)")
    if result["nonfinite"]:
        errors.append(f"{output.station}: {result['nonfinite']} non-finite row(s)")

    epochs = progress.get("epochs", [])
    if epochs:
        expected_first = gps_sow(experiment.begin)
        expected_last = gps_sow(experiment.end)
        if abs(epochs[0] - expected_first) > 1e-5:
            errors.append(f"{output.station}: first SOW {epochs[0]} != {expected_first}")
        if abs(epochs[-1] - expected_last) > 1e-5:
            errors.append(f"{output.station}: last SOW {epochs[-1]} != {expected_last}")
        gaps = [
            (left, right)
            for left, right in zip(epochs, epochs[1:])
            if abs((right - left) - experiment.interval_s) > 1e-5
        ]
        result["gap_count"] = len(gaps)
        if gaps:
            errors.append(f"{output.station}: {len(gaps)} duplicate/gap interval(s)")
    return result, errors


def gps_sow(value: dt.datetime) -> float:
    return (value - GPS_EPOCH).total_seconds() % 604800.0


def scan_patterns(paths: Iterable[Path], patterns: Sequence[re.Pattern]) -> List[Dict[str, object]]:
    matches: List[Dict[str, object]] = []
    seen_paths = set()
    for path in paths:
        if path in seen_paths or not path.is_file():
            continue
        seen_paths.add(path)
        try:
            with path.open("r", encoding="utf-8", errors="replace") as stream:
                for line_number, line in enumerate(stream, 1):
                    for pattern in patterns:
                        if pattern.search(line):
                            matches.append(
                                {
                                    "path": str(path),
                                    "line": line_number,
                                    "pattern": pattern.pattern,
                                    "text": line.strip()[:500],
                                }
                            )
                            break
                    if len(matches) >= 100:
                        return matches
        except OSError:
            continue
    return matches


def parse_runtime_diagnostics(paths: Iterable[Path]) -> Dict[str, object]:
    """Extract optional GREAT/FGO timing summaries from fresh process logs."""
    spent: List[float] = []
    profiles: List[Dict[str, object]] = []
    event_counts = {
        "feedback_accepted": 0,
        "feedback_rejected": 0,
        "outlier": 0,
        "cycle_slip": 0,
        "raw_prepare_failed": 0,
        "covariance_fallback": 0,
        "rank_deficient": 0,
        "pseudo_inverse": 0,
    }
    feedback_events: List[Dict[str, object]] = []
    phase_pattern = re.compile(
        r"^\s*(?P<phase>[A-Za-z][A-Za-z0-9_. ]*?)\s*\|\s*"
        r"(?P<count>\d+)\s*\|\s*(?P<avg>[-+0-9.eE]+)\s*\|\s*"
        r"(?P<min>[-+0-9.eE]+)\s*\|\s*(?P<max>[-+0-9.eE]+)\s*\|\s*"
        r"(?P<total>[-+0-9.eE]+)\s*$"
    )
    profile_header = re.compile(r"\[FGO\] per-window phase timing over (\d+) processed window")
    spent_pattern = re.compile(r"Spent\s*([0-9]+(?:\.[0-9]+)?)\s*seconds", re.IGNORECASE)
    raw_average = re.compile(
        r"RAW graph averages:\s*([0-9.eE+-]+)\s+build/solve attempt\(s\) per window,\s*"
        r"([0-9.eE+-]+)\s+parameter scalar\(s\),\s*"
        r"([0-9.eE+-]+)\s+residual scalar\(s\),\s*"
        r"([0-9.eE+-]+)\s+Ceres iteration\(s\) per solve"
    )
    event_patterns = {
        "feedback_accepted": re.compile(r"feedback\s+accepted", re.IGNORECASE),
        "feedback_rejected": re.compile(r"feedback\s+rejected|rejected\s+fixed\s+candidate", re.IGNORECASE),
        "outlier": re.compile(r"outlier", re.IGNORECASE),
        "cycle_slip": re.compile(r"cycle\s*slip|slip\s+detect", re.IGNORECASE),
        "raw_prepare_failed": re.compile(r"raw[^\r\n]*prepar[^\r\n]*failed", re.IGNORECASE),
        "covariance_fallback": re.compile(r"covariance[^\r\n]*(?:fallback|regulariz|ridge)", re.IGNORECASE),
        "rank_deficient": re.compile(r"rank\s*(?:deficient|=\s*N-\d+)", re.IGNORECASE),
        "pseudo_inverse": re.compile(r"pseudo(?:inverse|-inverse)", re.IGNORECASE),
    }
    feedback_pattern = re.compile(
        r"feedback\s+(?P<status>accepted|rejected)[^\r\n]*?"
        r"(?:selected_candidates=(?P<selected>\d+)/(?:\s*)?(?P<total>\d+))?[^\r\n]*?"
        r"(?:NIS=(?P<nis>[-+0-9.eE]+)/(?P<nis_limit>[-+0-9.eE]+))?[^\r\n]*?"
        r"(?:delta_cost=(?P<delta>[-+0-9.eE]+)\s+limit=(?P<delta_limit>[-+0-9.eE]+))?",
        re.IGNORECASE,
    )
    seen_paths = set()
    for path in paths:
        if path in seen_paths or not path.is_file():
            continue
        seen_paths.add(path)
        current: Optional[Dict[str, object]] = None
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for line in lines:
            for name, pattern in event_patterns.items():
                if pattern.search(line):
                    event_counts[name] += 1
            feedback_match = feedback_pattern.search(line)
            if feedback_match:
                event: Dict[str, object] = {"status": feedback_match.group("status").upper()}
                for field, group in (
                    ("selected_candidates", "selected"),
                    ("candidate_total", "total"),
                    ("nis", "nis"),
                    ("nis_limit", "nis_limit"),
                    ("delta_cost", "delta"),
                    ("delta_cost_limit", "delta_limit"),
                ):
                    value = feedback_match.group(group)
                    if value is not None:
                        event[field] = int(value) if field in ("selected_candidates", "candidate_total") else float(value)
                feedback_events.append(event)
            spent_match = spent_pattern.search(line)
            if spent_match:
                spent.append(float(spent_match.group(1)))
            header_match = profile_header.search(line)
            if header_match:
                if current is not None:
                    profiles.append(current)
                current = {
                    "windows": int(header_match.group(1)),
                    "phases": {},
                }
                continue
            if current is None:
                continue
            phase_match = phase_pattern.match(line)
            if phase_match:
                phase = phase_match.group("phase").strip()
                current["phases"][phase] = {
                    "count": int(phase_match.group("count")),
                    "avg_ms": float(phase_match.group("avg")),
                    "min_ms": float(phase_match.group("min")),
                    "max_ms": float(phase_match.group("max")),
                    "total_ms": float(phase_match.group("total")),
                }
                continue
            raw_match = raw_average.search(line)
            if raw_match:
                current["raw_graph_averages"] = {
                    "attempts_per_window": float(raw_match.group(1)),
                    "parameters": float(raw_match.group(2)),
                    "residuals": float(raw_match.group(3)),
                    "ceres_iterations_per_solve": float(raw_match.group(4)),
                }
        if current is not None:
            profiles.append(current)
    result: Dict[str, object] = {
        "spent_seconds": spent[-1] if spent else None,
        "spent_samples": spent,
        "fgo_profile_count": len(profiles),
        "fgo_profiles": profiles,
        "event_counts": event_counts,
        "feedback_events": feedback_events,
    }
    if len(profiles) == 1:
        result["fgo_profile"] = profiles[0]
    return result


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_provenance() -> Dict[str, object]:
    def run_git(*arguments: str) -> Optional[str]:
        try:
            completed = subprocess.run(
                ["git", *arguments], cwd=REPO_ROOT, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=10,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        return completed.stdout.strip() if completed.returncode == 0 else None

    status = run_git("status", "--porcelain")
    return {
        "commit": run_git("rev-parse", "HEAD"),
        "branch": run_git("branch", "--show-current"),
        "dirty": bool(status) if status is not None else None,
    }


def experiment_description(experiment: Experiment) -> Dict[str, object]:
    return {
        "name": experiment.name,
        "config": str(experiment.config),
        "workdir": str(experiment.workdir),
        "executable": str(experiment.executable),
        "stations": experiment.stations,
        "begin": experiment.begin.isoformat(sep=" "),
        "end": experiment.end.isoformat(sep=" "),
        "interval_s": experiment.interval_s,
        "expected_rows_per_station": experiment.expected_rows,
        "product": experiment.product,
        "mode": experiment.mode,
        "frequency": experiment.frequency,
        "solver_threads": experiment.solver_threads,
        "arguments": experiment.arguments,
        "tags": experiment.tags,
        "outputs": [
            {"station": item.station, "fgo": str(item.fgo), "ppp": str(item.ppp) if item.ppp else None}
            for item in experiment.outputs
        ],
        "diagnostic_logs": [str(path) for path in experiment.diagnostic_logs],
    }


def build_initial_state(run_id: str, run_dir: Path, experiments: Sequence[Experiment]) -> Dict[str, object]:
    executables = sorted({experiment.executable for experiment in experiments})
    return {
        "schema_version": 1,
        "run_id": run_id,
        "run_dir": str(run_dir),
        "runner_pid": os.getpid(),
        "status": "PENDING",
        "started_at": dt.datetime.now().astimezone().isoformat(),
        "updated_at": dt.datetime.now().astimezone().isoformat(),
        "host": {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "processor": platform.processor(),
            "logical_cpus": os.cpu_count(),
            "python": sys.version.split()[0],
        },
        "git": git_provenance(),
        "executables": [
            {
                "path": str(path),
                "size": path.stat().st_size,
                "mtime": dt.datetime.fromtimestamp(path.stat().st_mtime).astimezone().isoformat(),
                "sha256": sha256(path),
            }
            for path in executables
        ],
        "experiments": {
            experiment.name: {
                "status": "PENDING",
                "definition": experiment_description(experiment),
                "progress": {},
            }
            for experiment in experiments
        },
    }


def update_progress(state: dict, active: Mapping[int, ProcessEntry]) -> None:
    now = time.monotonic()
    for entry in active.values():
        sample_process_resources(entry)
        progress = {}
        for output in entry.experiment.outputs:
            snapshot = entry_output_progress(entry, output)
            progress[output.station] = {
                key: value for key, value in snapshot.items() if key != "epochs"
            }
            progress[output.station]["expected_rows"] = entry.experiment.expected_rows
        item = state["experiments"][entry.experiment.name]
        item["progress"] = progress
        item["elapsed_seconds"] = round(now - entry.started_wall, 3)
    state["updated_at"] = dt.datetime.now().astimezone().isoformat()


def status_line(entry: ProcessEntry) -> str:
    parts = []
    for output in entry.experiment.outputs:
        progress = entry_output_progress(entry, output)
        parts.append(f"{output.station}={progress['rows']}/{entry.experiment.expected_rows}")
    elapsed = time.monotonic() - entry.started_wall
    return f"RUN   {entry.experiment.name:<28} {elapsed:7.1f}s  " + " ".join(parts)


def start_experiment(experiment: Experiment, run_dir: Path, state: dict) -> ProcessEntry:
    log_dir = run_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = log_dir / f"{experiment.name}.stdout.log"
    stderr_path = log_dir / f"{experiment.name}.stderr.log"
    stdout_stream = stdout_path.open("w", encoding="utf-8")
    stderr_stream = stderr_path.open("w", encoding="utf-8")
    config_argument = os.path.relpath(experiment.config, experiment.workdir)
    command = [str(experiment.executable), "-x", config_argument, *experiment.arguments]
    flags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    # Snapshot before launch so even a very fast truncate/write is recognized.
    output_before: Dict[str, Optional[Tuple[int, int]]] = {}
    for output in experiment.outputs:
        output_before[str(output.fgo)] = file_signature(output.fgo)
        if output.ppp:
            output_before[str(output.ppp)] = file_signature(output.ppp)
    for path in experiment.diagnostic_logs:
        output_before[str(path)] = file_signature(path)
    try:
        process = subprocess.Popen(
            command,
            cwd=experiment.workdir,
            stdout=stdout_stream,
            stderr=stderr_stream,
            creationflags=flags,
        )
    except Exception:
        stdout_stream.close()
        stderr_stream.close()
        raise
    started_iso = dt.datetime.now().astimezone().isoformat()
    entry = ProcessEntry(
        experiment=experiment,
        process=process,
        stdout_stream=stdout_stream,
        stderr_stream=stderr_stream,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        started_wall=time.monotonic(),
        started_iso=started_iso,
        output_before=output_before,
    )
    item = state["experiments"][experiment.name]
    item.update(
        {
            "status": "RUNNING",
            "pid": process.pid,
            "started_at": started_iso,
            "command": command,
            "stdout": str(stdout_path),
            "stderr": str(stderr_path),
        }
    )
    print(f"START {experiment.name} pid={process.pid}", flush=True)
    return entry


def finish_experiment(
    entry: ProcessEntry,
    state: dict,
    compiled_patterns: Sequence[re.Pattern],
) -> Dict[str, object]:
    sample_process_resources(entry)
    entry.process.wait()
    entry.stdout_stream.close()
    entry.stderr_stream.close()
    elapsed = time.monotonic() - entry.started_wall
    output_results: Dict[str, object] = {}
    errors: List[str] = []
    for output in entry.experiment.outputs:
        validation, output_errors = validate_output(
            entry.experiment, output, entry.output_before[str(output.fgo)]
        )
        output_results[output.station] = validation
        errors.extend(output_errors)

    scan_paths = [entry.stdout_path, entry.stderr_path]
    scan_paths.extend(
        output.ppp
        for output in entry.experiment.outputs
        if output.ppp and file_signature(output.ppp) != entry.output_before.get(str(output.ppp))
    )
    scan_paths.extend(
        path
        for path in entry.experiment.diagnostic_logs
        if file_signature(path) != entry.output_before.get(str(path))
    )
    fatal_matches = scan_patterns(scan_paths, compiled_patterns)
    if entry.process.returncode != 0:
        errors.append(f"process exit code {entry.process.returncode}")
    if entry.timed_out:
        errors.append("process exceeded timeout")
    if fatal_matches:
        errors.append(f"fatal log patterns matched {len(fatal_matches)} time(s)")

    total_rows = sum(int(item["rows"]) for item in output_results.values())
    data_seconds = total_rows * entry.experiment.interval_s
    performance = {
        "wall_seconds": round(elapsed, 3),
        "station_rows": total_rows,
        "station_epochs_per_second": total_rows / elapsed if elapsed > 0 else None,
        "processed_data_seconds": data_seconds,
        "realtime_factor": data_seconds / elapsed if elapsed > 0 else None,
        "solver_threads": entry.experiment.solver_threads,
        "peak_rss_mb": (
            entry.peak_rss_bytes / (1024.0 * 1024.0)
            if entry.peak_rss_bytes else None
        ),
        "last_rss_mb": (
            entry.last_rss_bytes / (1024.0 * 1024.0)
            if entry.last_rss_bytes else None
        ),
        "runtime_diagnostics": parse_runtime_diagnostics(
            [
                entry.stdout_path,
                entry.stderr_path,
                *(
                    path
                    for path in entry.experiment.diagnostic_logs
                    if file_signature(path) != entry.output_before.get(str(path))
                ),
            ]
        ),
    }
    status = "PASSED" if not errors else "FAILED"
    result = {
        "status": status,
        "exit_code": entry.process.returncode,
        "timed_out": entry.timed_out,
        "started_at": entry.started_iso,
        "finished_at": dt.datetime.now().astimezone().isoformat(),
        "outputs": output_results,
        "fatal_log_matches": fatal_matches,
        "errors": errors,
        "performance": performance,
    }
    state["experiments"][entry.experiment.name].update(result)
    print(
        f"{status:<6} {entry.experiment.name:<28} exit={entry.process.returncode} "
        f"wall={elapsed:.2f}s rows={total_rows}",
        flush=True,
    )
    for error in errors:
        print(f"       - {error}", flush=True)
    return result


def performance_gate_config(manifest: Mapping[str, object]) -> Dict[str, object]:
    configured = manifest.get("performance_gates")
    if isinstance(configured, Mapping):
        return dict(configured)
    gates = manifest.get("gates")
    if isinstance(gates, Mapping) and isinstance(gates.get("performance"), Mapping):
        return dict(gates["performance"])
    return {}


def apply_performance_gates(state: dict, manifest: Mapping[str, object]) -> bool:
    """Apply optional runtime gates after all child processes have finished."""
    gates = performance_gate_config(manifest)
    if not gates:
        return False
    any_failed = False
    phase_avg_limits = gates.get("max_profile_avg_ms", {})
    phase_max_limits = gates.get("max_profile_max_ms", {})
    if not isinstance(phase_avg_limits, Mapping):
        phase_avg_limits = {}
    if not isinstance(phase_max_limits, Mapping):
        phase_max_limits = {}
    for name, item in state.get("experiments", {}).items():
        performance = item.get("performance", {})
        diagnostics = performance.get("runtime_diagnostics", {})
        failures: List[str] = []

        def check_max(key: str, label: str) -> None:
            limit = gates.get(key)
            if limit is None:
                return
            if key == "max_wall_seconds":
                value = performance.get("wall_seconds")
            elif key == "max_spent_seconds":
                value = diagnostics.get("spent_seconds")
            else:
                value = performance.get(key)
            if value is None or float(value) > float(limit):
                failures.append(f"{label}={value} exceeds {limit}")

        check_max("max_wall_seconds", "wall seconds")
        check_max("max_spent_seconds", "Spent seconds")
        limit_rss = gates.get("max_peak_rss_mb")
        if limit_rss is not None:
            value_rss = performance.get("peak_rss_mb")
            if value_rss is None or float(value_rss) > float(limit_rss):
                failures.append(f"peak RSS MB={value_rss} exceeds {limit_rss}")
        for key, label in (
            ("min_realtime_factor", "realtime factor"),
            ("min_station_epochs_per_second", "station epochs/s"),
        ):
            limit = gates.get(key)
            if limit is None:
                continue
            value_key = "realtime_factor" if key == "min_realtime_factor" else "station_epochs_per_second"
            value = performance.get(value_key)
            if value is None or float(value) < float(limit):
                failures.append(f"{label}={value} below {limit}")

        profiles = diagnostics.get("fgo_profiles", [])
        event_counts = diagnostics.get("event_counts", {})
        for gate_key, event_key, label in (
            ("max_feedback_rejected", "feedback_rejected", "feedback rejected"),
            ("max_covariance_fallback", "covariance_fallback", "covariance fallback"),
            ("max_rank_deficient", "rank_deficient", "rank deficient"),
            ("max_pseudo_inverse", "pseudo_inverse", "pseudo-inverse"),
            ("max_raw_prepare_failed", "raw_prepare_failed", "RAW preparation failures"),
        ):
            limit = gates.get(gate_key)
            if limit is not None:
                value = event_counts.get(event_key, 0)
                if value > int(limit):
                    failures.append(f"{label}={value} exceeds {limit}")
        if (phase_avg_limits or phase_max_limits or
                gates.get("max_profile_attempts_per_window") is not None or
                gates.get("max_profile_iterations") is not None):
            if not profiles:
                failures.append("FGO profile is missing")
            for profile in profiles:
                phases = profile.get("phases", {})
                for phase, limit in phase_avg_limits.items():
                    value = phases.get(phase, {}).get("avg_ms")
                    if value is None or float(value) > float(limit):
                        failures.append(f"{phase} avg ms={value} exceeds {limit}")
                for phase, limit in phase_max_limits.items():
                    value = phases.get(phase, {}).get("max_ms")
                    if value is None or float(value) > float(limit):
                        failures.append(f"{phase} max ms={value} exceeds {limit}")
                averages = profile.get("raw_graph_averages", {})
                for key, label in (
                    ("max_profile_attempts_per_window", "RAW attempts/window"),
                    ("max_profile_iterations", "Ceres iterations/solve"),
                ):
                    limit = gates.get(key)
                    if limit is None:
                        continue
                    value_key = (
                        "attempts_per_window"
                        if key == "max_profile_attempts_per_window"
                        else "ceres_iterations_per_solve"
                    )
                    value = averages.get(value_key)
                    if value is None or float(value) > float(limit):
                        failures.append(f"{label}={value} exceeds {limit}")

        gate_result = {
            "enabled": True,
            "passed": not failures,
            "limits": gates,
            "failures": failures,
        }
        item["performance_gate"] = gate_result
        if failures:
            any_failed = True
            item["status"] = "FAILED"
            item.setdefault("errors", []).extend(f"performance gate: {failure}" for failure in failures)
    return any_failed


def terminate_process(entry: ProcessEntry) -> None:
    if entry.process.poll() is not None:
        return
    entry.process.terminate()
    try:
        entry.process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        entry.process.kill()
        entry.process.wait(timeout=10)


def run_analysis(
    experiments: Sequence[Experiment],
    state: dict,
    references: Mapping[str, Tuple[float, float, float]],
    run_dir: Path,
    analyzer: Path,
    manifest: Mapping[str, object],
) -> List[Dict[str, object]]:
    if not references:
        return []
    if not analyzer.is_file():
        return [{"status": "FAILED", "error": f"analyzer not found: {analyzer}"}]

    # The existing analyzer accepts cases sharing one nominal time grid.  Group
    # by station and grid so mixed-duration matrices remain valid.
    groups: Dict[Tuple[str, dt.datetime, dt.datetime, float], List[Tuple[Experiment, OutputSpec]]] = {}
    for experiment in experiments:
        if state["experiments"][experiment.name]["status"] != "PASSED":
            continue
        for output in experiment.outputs:
            if output.station not in references:
                continue
            key = (output.station, experiment.begin, experiment.end, experiment.interval_s)
            groups.setdefault(key, []).append((experiment, output))

    analysis_options = dict(manifest.get("analysis", {}))
    quality_options = analysis_options.get("quality_gates", {})
    if isinstance(quality_options, Mapping):
        analysis_options.update(quality_options)
    baseline_map = {str(key).upper(): value for key, value in dict(analysis_options.get("baselines", {})).items()}
    analysis_dir = run_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    reports: List[Dict[str, object]] = []
    for index, (key, cases) in enumerate(sorted(groups.items(), key=lambda item: str(item[0]))):
        station, begin, end, interval_s = key
        label = f"{station}_{begin:%Y%m%dT%H%M%S}_{index + 1}"
        markdown_path = analysis_dir / f"{label}.md"
        json_path = analysis_dir / f"{label}.json"
        command = [
            sys.executable,
            str(analyzer),
            "--site", station,
            "--reference", *(str(value) for value in references[station]),
            "--nominal-start-sow", str(gps_sow(begin)),
            "--interval", str(interval_s),
            "--hours", str(((end - begin).total_seconds() + interval_s) / 3600.0),
        ]
        for experiment, output in cases:
            command.extend(["--case", f"{experiment.name}={output.fgo}"])
        if station in baseline_map:
            command.extend(["--baseline", str(resolve_path(str(baseline_map[station]), REPO_ROOT))])
        for option, cli_name in (
            ("warmup_minutes", "--warmup-minutes"),
            ("window_minutes", "--window-minutes"),
            ("horizontal_threshold", "--horizontal-threshold"),
            ("vertical_threshold", "--vertical-threshold"),
            ("jump_threshold", "--jump-threshold"),
            ("three_d_anomaly_threshold", "--three-d-anomaly-threshold"),
            ("min_nsat", "--min-nsat"),
            ("baseline_delta_threshold", "--baseline-delta-threshold"),
            ("max_sustained_minutes", "--max-sustained-minutes"),
            ("max_permanent_minutes", "--max-permanent-minutes"),
            ("min_fixed_fraction", "--min-fixed-fraction"),
            ("min_correct_fixed_fraction", "--min-correct-fixed-fraction"),
            ("max_post_convergence_3d", "--max-post-convergence-3d"),
            ("max_coordinate_jumps", "--max-coordinate-jumps"),
            ("max_status_transitions", "--max-status-transitions"),
            ("max_baseline_delta", "--max-baseline-delta"),
            ("max_baseline_delta_exceedances", "--max-baseline-delta-exceedances"),
        ):
            if option in analysis_options:
                command.extend([cli_name, str(analysis_options[option])])
        if analysis_options.get("require_complete"):
            command.append("--require-complete")
        if analysis_options.get("fail_on_gate"):
            command.append("--fail-on-gate")
        command.extend(["--markdown", str(markdown_path), "--json", str(json_path), "--quiet"])
        completed = subprocess.run(command, cwd=REPO_ROOT, text=True, capture_output=True, check=False)
        report: Dict[str, object] = {
            "station": station,
            "cases": [experiment.name for experiment, _ in cases],
            "status": "PASSED" if completed.returncode == 0 else "FAILED",
            "markdown": str(markdown_path),
            "json": str(json_path),
            "exit_code": completed.returncode,
        }
        if json_path.is_file():
            try:
                payload = json.loads(json_path.read_text(encoding="utf-8"))
                gate_failures = []
                for analyzed_case in payload.get("cases", []):
                    gate = analyzed_case.get("quality_gate", {})
                    for failure in gate.get("failures", []):
                        gate_failures.append(
                            f"{analyzed_case.get('label', 'case')}: {failure}"
                        )
                if gate_failures:
                    report["quality_gate_failures"] = gate_failures
            except (OSError, json.JSONDecodeError):
                pass
        if completed.returncode != 0:
            report["error"] = (completed.stderr or completed.stdout).strip()[-2000:]
        reports.append(report)
    return reports


def minutes_from_start(epoch: Optional[float], nominal_start: float) -> Optional[float]:
    if epoch is None:
        return None
    delta = epoch - nominal_start
    if delta < -302400:
        delta += 604800
    elif delta > 302400:
        delta -= 604800
    return delta / 60.0


def analysis_rows(reports: Sequence[Mapping[str, object]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for report in reports:
        try:
            payload = json.loads(Path(str(report["json"])).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        nominal_start = float(payload["configuration"]["nominal_start_sow"])
        for case in payload.get("cases", []):
            accuracy = case.get("accuracy_after_permanent_convergence") or case.get("accuracy_all", {})
            convergence = case.get("convergence", {})
            rows.append(
                {
                    "site": payload.get("site"),
                    "case": case.get("label"),
                    "complete": case.get("continuity", {}).get("complete_nominal_span"),
                    "fixed_fraction": case.get("accuracy_all", {}).get("fixed_fraction"),
                    "first_sustained_min": minutes_from_start(
                        convergence.get("first_sustained_convergence_epoch"), nominal_start
                    ),
                    "permanent_min": minutes_from_start(
                        convergence.get("permanent_convergence_epoch"), nominal_start
                    ),
                    "first_correct_fixed_min": minutes_from_start(
                        convergence.get("first_correct_fixed_epoch"), nominal_start
                    ),
                    "rms_3d_m": accuracy.get("three_d_rms_m"),
                    "mae_3d_m": accuracy.get("three_d_mae_m"),
                    "p95_3d_m": accuracy.get("three_d_p95_m"),
                    "p99_3d_m": accuracy.get("three_d_p99_m"),
                    "max_3d_m": accuracy.get("three_d_max_m"),
                }
            )
    return rows


def fmt(value: object, digits: int = 3) -> str:
    if value is None:
        return "N/A"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int, float)):
        return f"{value:.{digits}f}"
    return str(value)


def pid_is_alive(pid: object) -> bool:
    try:
        numeric_pid = int(pid)
        if numeric_pid <= 0:
            return False
        os.kill(numeric_pid, 0)
    except (OSError, TypeError, ValueError):
        return False
    return True


def markdown_report(state: Mapping[str, object], analysis_reports: Sequence[Mapping[str, object]]) -> str:
    experiments = state["experiments"]
    passed = sum(item["status"] == "PASSED" for item in experiments.values())
    failed = sum(item["status"] == "FAILED" for item in experiments.values())
    lines = [
        f"# FGO experiment report: {state['run_id']}",
        "",
        f"- Overall status: **{state['status']}**",
        f"- Experiments: {passed} passed, {failed} failed, {len(experiments)} total",
        f"- Started: {state['started_at']}",
        f"- Finished: {state.get('finished_at', 'N/A')}",
        f"- Git: `{state.get('git', {}).get('commit')}` (dirty={state.get('git', {}).get('dirty')})",
        "",
        "## Execution and performance",
        "",
        "| Experiment | Product | Freq | Mode | Stations | Status | Wall s | Spent s | Peak RSS MB | Rows | Rows/s | Realtime x | FGO profiles |",
        "|---|---|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, item in experiments.items():
        definition = item["definition"]
        performance = item.get("performance", {})
        runtime = performance.get("runtime_diagnostics", {})
        lines.append(
            f"| {name} | {definition['product']} | {fmt(definition['frequency'], 0)} | "
            f"{definition['mode']} | {len(definition['stations'])} | {item['status']} | "
            f"{fmt(performance.get('wall_seconds'))} | {fmt(runtime.get('spent_seconds'))} | "
            f"{fmt(performance.get('peak_rss_mb'))} | {fmt(performance.get('station_rows'), 0)} | "
            f"{fmt(performance.get('station_epochs_per_second'))} | "
            f"{fmt(performance.get('realtime_factor'))} | {fmt(runtime.get('fgo_profile_count'), 0)} |"
        )
    matrix = state.get("performance", {})
    lines.extend(
        [
            "",
            f"Matrix wall time: **{fmt(matrix.get('matrix_wall_seconds'))} s**; "
            f"aggregate throughput: **{fmt(matrix.get('aggregate_station_epochs_per_second'))} station epochs/s**.",
            "",
        ]
    )

    accuracy = analysis_rows(analysis_reports)
    if accuracy:
        lines.extend(
            [
                "## Accuracy and convergence",
                "",
                "Statistics use the existing analyzer. Accuracy is taken after permanent convergence when available.",
                "",
                "| Site | Experiment | Complete | Fixed | Sustained min | Permanent min | First correct fixed min | 3D RMS | MAE | P95 | P99 | Max (m) |",
                "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for row in accuracy:
            fixed = row["fixed_fraction"] * 100 if row["fixed_fraction"] is not None else None
            fixed_text = f"{fmt(fixed, 1)}%" if fixed is not None else "N/A"
            lines.append(
                f"| {row['site']} | {row['case']} | {fmt(row['complete'])} | {fixed_text} | "
                f"{fmt(row['first_sustained_min'], 1)} | {fmt(row['permanent_min'], 1)} | "
                f"{fmt(row['first_correct_fixed_min'], 1)} | {fmt(row['rms_3d_m'], 4)} | "
                f"{fmt(row['mae_3d_m'], 4)} | {fmt(row['p95_3d_m'], 4)} | "
                f"{fmt(row['p99_3d_m'], 4)} | {fmt(row['max_3d_m'], 4)} |"
            )
        lines.append("")

    failures = [(name, item) for name, item in experiments.items() if item["status"] == "FAILED"]
    if failures:
        lines.extend(["## Failures", ""])
        for name, item in failures:
            lines.append(f"### {name}")
            lines.append("")
            for error in item.get("errors", []):
                lines.append(f"- {error}")
            for match in item.get("fatal_log_matches", [])[:10]:
                lines.append(f"- `{match['path']}:{match['line']}`: {match['text']}")
            lines.append("")

    analysis_failures = [
        report for report in analysis_reports
        if report.get("status") == "FAILED"
    ]
    if analysis_failures:
        lines.extend(["## Analysis gate failures", ""])
        for report in analysis_failures:
            lines.append(f"### {report.get('station', 'analysis')}")
            lines.append("")
            for failure in report.get("quality_gate_failures", []):
                lines.append(f"- {failure}")
            if report.get("error"):
                lines.append(f"- {report['error']}")
            lines.append("")

    lines.extend(["## Artifacts", ""])
    lines.append(f"- Machine-readable state and results: `{Path(str(state['run_dir'])) / 'state.json'}`")
    for report in analysis_reports:
        if report.get("markdown"):
            lines.append(
                f"- {report['station']} accuracy report: `{report['markdown']}` "
                f"({report.get('status', 'UNKNOWN')})"
            )
    return "\n".join(lines) + "\n"


def print_plan(experiments: Sequence[Experiment], max_parallel: int) -> None:
    print(f"Validated {len(experiments)} experiment(s); max_parallel={max_parallel}")
    for experiment in experiments:
        outputs = ", ".join(f"{item.station}:{item.fgo}" for item in experiment.outputs)
        print(
            f"- {experiment.name}: {experiment.product} F{experiment.frequency or '?'} "
            f"{experiment.mode}, {experiment.begin}..{experiment.end}, "
            f"{experiment.expected_rows} rows/station, workdir={experiment.workdir}"
        )
        print(f"  outputs: {outputs}")


def run_command(args: argparse.Namespace) -> int:
    manifest, experiments, references = read_manifest(args)
    max_parallel = int(args.max_parallel or manifest.get("max_parallel", 1))
    poll_seconds = float(args.poll_seconds or manifest.get("poll_seconds", 1.0))
    status_seconds = float(args.status_seconds or manifest.get("status_seconds", 10.0))
    timeout_seconds = float(args.timeout_seconds or manifest.get("timeout_seconds", 0.0))
    if max_parallel < 1 or poll_seconds <= 0 or status_seconds <= 0 or timeout_seconds < 0:
        raise ConfigurationError("parallelism and intervals must be positive; timeout may be zero")
    print_plan(experiments, max_parallel)
    busiest_threads = sum(
        sorted((item.solver_threads for item in experiments), reverse=True)[:max_parallel]
    )
    if os.cpu_count() and busiest_threads > os.cpu_count():
        print(
            f"WARNING active Ceres threads may reach {busiest_threads}, above "
            f"{os.cpu_count()} logical CPUs",
            flush=True,
        )
    if args.dry_run:
        return 0

    run_id = args.run_id or dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    if Path(run_id).name != run_id or run_id in (".", ".."):
        raise ConfigurationError("run-id must be a single directory name")
    run_root_value = args.run_root or manifest.get("run_root", "build/fgo_runs")
    run_dir = resolve_path(run_root_value) / run_id
    if run_dir.exists() and any(run_dir.iterdir()):
        raise ConfigurationError(f"run directory is not empty: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "manifest.resolved.json").write_text(
        json.dumps(
            {
                "source": manifest,
                "experiments": [experiment_description(item) for item in experiments],
                "references": references,
            },
            ensure_ascii=False,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    state = build_initial_state(run_id, run_dir, experiments)
    state_path = run_dir / "state.json"
    state["status"] = "RUNNING"
    atomic_json(state_path, state)

    pattern_texts = manifest.get("failure_patterns", DEFAULT_FAILURE_PATTERNS)
    try:
        compiled_patterns = [re.compile(str(value), re.IGNORECASE) for value in pattern_texts]
    except re.error as exc:
        raise ConfigurationError(f"invalid failure pattern: {exc}") from exc

    pending = list(experiments)
    active: Dict[int, ProcessEntry] = {}
    matrix_started = time.monotonic()
    next_status = matrix_started
    interrupted = False
    try:
        while pending or active:
            while pending and len(active) < max_parallel:
                experiment = pending.pop(0)
                try:
                    entry = start_experiment(experiment, run_dir, state)
                except OSError as exc:
                    item = state["experiments"][experiment.name]
                    item.update({"status": "FAILED", "errors": [f"start failed: {exc}"]})
                    print(f"FAILED {experiment.name}: start failed: {exc}", flush=True)
                    continue
                active[entry.process.pid] = entry
                atomic_json(state_path, state)

            now = time.monotonic()
            for pid, entry in list(active.items()):
                if timeout_seconds and now - entry.started_wall > timeout_seconds and entry.process.poll() is None:
                    entry.timed_out = True
                    terminate_process(entry)
                if entry.process.poll() is None:
                    continue
                finish_experiment(entry, state, compiled_patterns)
                del active[pid]
                atomic_json(state_path, state)

            update_progress(state, active)
            if now >= next_status and active:
                for entry in active.values():
                    print(status_line(entry), flush=True)
                next_status = now + status_seconds
                atomic_json(state_path, state)
            if pending or active:
                time.sleep(poll_seconds)
    except KeyboardInterrupt:
        interrupted = True
        print("Interrupt received; terminating active experiments...", flush=True)
        for entry in active.values():
            terminate_process(entry)
            entry.stdout_stream.close()
            entry.stderr_stream.close()
            state["experiments"][entry.experiment.name].update(
                {"status": "INTERRUPTED", "errors": ["runner interrupted"]}
            )
        for experiment in pending:
            state["experiments"][experiment.name]["status"] = "NOT_RUN"
    except Exception as exc:
        for entry in active.values():
            terminate_process(entry)
            entry.stdout_stream.close()
            entry.stderr_stream.close()
            state["experiments"][entry.experiment.name].update(
                {"status": "ERROR", "errors": [f"runner error: {exc}"]}
            )
        state["status"] = "ERROR"
        state["updated_at"] = dt.datetime.now().astimezone().isoformat()
        atomic_json(state_path, state)
        raise

    matrix_wall = time.monotonic() - matrix_started
    completed_items = list(state["experiments"].values())
    total_rows = sum(
        item.get("performance", {}).get("station_rows", 0) for item in completed_items
    )
    state["performance"] = {
        "matrix_wall_seconds": round(matrix_wall, 3),
        "aggregate_station_rows": total_rows,
        "aggregate_station_epochs_per_second": total_rows / matrix_wall if matrix_wall > 0 else None,
        "max_parallel": max_parallel,
        "poll_seconds": poll_seconds,
    }

    performance_gate_failed = apply_performance_gates(state, manifest)

    analyzer_value = args.analyzer or manifest.get("analyzer", DEFAULT_ANALYZER)
    analysis_reports = [] if interrupted else run_analysis(
        experiments,
        state,
        references,
        run_dir,
        resolve_path(str(analyzer_value)),
        manifest,
    )
    analysis_failed = any(report.get("status") == "FAILED" for report in analysis_reports)
    execution_failed = performance_gate_failed or any(
        item["status"] not in ("PASSED",) for item in state["experiments"].values()
    )
    state["analysis"] = analysis_reports
    state["status"] = "INTERRUPTED" if interrupted else ("FAILED" if execution_failed or analysis_failed else "PASSED")
    state["finished_at"] = dt.datetime.now().astimezone().isoformat()
    state["updated_at"] = state["finished_at"]
    atomic_json(state_path, state)
    report_path = run_dir / "report.md"
    report_path.write_text(markdown_report(state, analysis_reports), encoding="utf-8")
    print(f"REPORT {report_path}")
    print(f"OVERALL {state['status']} wall={matrix_wall:.2f}s")
    return 0 if state["status"] == "PASSED" else 1


def status_command(args: argparse.Namespace) -> int:
    run_dir = resolve_path(args.run_dir)
    state_path = run_dir / "state.json"
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"cannot read state {state_path}: {exc}") from exc
    saved_status = state.get("status")
    runner_alive = pid_is_alive(state.get("runner_pid")) if saved_status == "RUNNING" else None
    display_status = "STALE" if saved_status == "RUNNING" and not runner_alive else saved_status
    print(
        f"Run {state.get('run_id')} status={display_status} updated={state.get('updated_at')}"
        + (f" runner_alive={runner_alive}" if runner_alive is not None else "")
    )
    if display_status == "STALE":
        print("The saved runner is no longer alive; RUNNING is the last checkpoint, not an active test.")
    for name, item in state.get("experiments", {}).items():
        progress = item.get("progress", {})
        detail = " ".join(
            f"{site}={value.get('rows', 0)}/{value.get('expected_rows', '?')}"
            for site, value in progress.items()
        )
        wall = item.get("performance", {}).get("wall_seconds", item.get("elapsed_seconds"))
        print(f"{item.get('status', '?'):<11} {name:<28} wall={fmt(wall)}s {detail}")
        for error in item.get("errors", []):
            print(f"  - {error}")
    report_path = run_dir / "report.md"
    if report_path.is_file():
        print(f"Report: {report_path}")
    return 1 if display_status in ("FAILED", "INTERRUPTED", "ERROR", "STALE") else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="run a single experiment or matrix")
    run_parser.add_argument("--manifest", type=Path, help="JSON experiment matrix")
    run_parser.add_argument(
        "--case", action="append", type=lambda value: parse_named_path(value, "case"),
        metavar="NAME=CONFIG", help="single/direct experiment; repeat for a small matrix",
    )
    run_parser.add_argument(
        "--work-dir", action="append", type=lambda value: parse_named_path(value, "work-dir"),
        metavar="NAME=DIR", help="working directory for a named --case",
    )
    run_parser.add_argument("--default-work-dir", type=Path)
    run_parser.add_argument("--executable", type=Path, help="executable for direct --case runs")
    run_parser.add_argument("--analyzer", type=Path, help="existing solution analyzer")
    run_parser.add_argument(
        "--reference", action="append", type=parse_reference, metavar="SITE=X,Y,Z",
        help="enable accuracy analysis for a station; repeat as needed",
    )
    run_parser.add_argument("--max-parallel", type=int)
    run_parser.add_argument("--poll-seconds", type=float)
    run_parser.add_argument("--status-seconds", type=float)
    run_parser.add_argument("--timeout-seconds", type=float)
    run_parser.add_argument("--run-root", type=Path)
    run_parser.add_argument("--run-id", help="unique output directory name")
    run_parser.add_argument("--dry-run", action="store_true", help="validate and print without launching")
    run_parser.set_defaults(handler=run_command)

    status_parser = subparsers.add_parser("status", help="show a saved or running matrix state")
    status_parser.add_argument("--run-dir", type=Path, required=True)
    status_parser.set_defaults(handler=status_command)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except ConfigurationError as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
