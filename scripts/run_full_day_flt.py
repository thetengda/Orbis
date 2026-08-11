#!/usr/bin/env python3
"""Run the isolated full-day FLT matrix and record progress/results."""

from __future__ import annotations

import datetime as dt
import json
import os
import re
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "build" / "full_day_tests_20260811" / "flt_runs"
CONFIG_ROOT = ROOT / "build" / "full_day_tests_20260811" / "configs"
EXE = ROOT / "build" / "Bin" / "Release" / "GREAT_PVT.exe"
CASES = {
    "FLT_UPD_DF": ROOT / "sample_data" / "PPPFLT_2023305",
    "FLT_UPD_FF": ROOT / "sample_data" / "PPPFLT_2023305",
    "FLT_OSB_DF": ROOT / "sample_data" / "PPPFLT_2023305_OSB",
    "FLT_OSB_FF": ROOT / "sample_data" / "PPPFLT_2023305_OSB",
}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat()


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def rows(path: Path) -> dict:
    result = {"exists": path.is_file(), "rows": 0, "malformed": 0, "nonfinite": 0, "first_sow": None, "last_sow": None}
    if not path.is_file():
        return result
    first = last = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if not fields or not re.fullmatch(r"[-+]?\d+(?:\.\d+)?", fields[0]):
            continue
        if len(fields) < 4:
            result["malformed"] += 1
            continue
        try:
            sow = float(fields[0])
            xyz = [float(item) for item in fields[1:4]]
        except ValueError:
            result["malformed"] += 1
            continue
        if not all(map(lambda item: item == item and abs(item) != float("inf"), [sow, *xyz])):
            result["nonfinite"] += 1
            continue
        result["rows"] += 1
        first = sow if first is None else first
        last = sow
    result["first_sow"] = first
    result["last_sow"] = last
    return result


def main() -> None:
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    log_root = RUN_ROOT / "logs"
    log_root.mkdir(parents=True, exist_ok=True)
    processes = {}
    state = {
        "status": "RUNNING",
        "started_at": now(),
        "executable": str(EXE),
        "cases": {},
    }
    for name, workdir in CASES.items():
        product, frequency = name.removeprefix("FLT_").split("_")
        config = CONFIG_ROOT / f"FLT_{product}_{frequency}.xml"
        output = workdir / "result" / f"CODEX_FULLDAY_FLT_{product}_{frequency}_{{station}}.flt"
        stdout = (log_root / f"{name}.stdout.log").open("w", encoding="utf-8")
        stderr = (log_root / f"{name}.stderr.log").open("w", encoding="utf-8")
        proc = subprocess.Popen([str(EXE), "-x", str(config)], cwd=workdir, stdout=stdout, stderr=stderr)
        processes[name] = (proc, stdout, stderr, workdir, product, frequency)
        state["cases"][name] = {
            "status": "RUNNING",
            "pid": proc.pid,
            "config": str(config),
            "workdir": str(workdir),
            "started_at": now(),
            "outputs": {},
        }
    atomic_json(RUN_ROOT / "state.json", state)

    try:
        while processes:
            for name, (proc, stdout, stderr, workdir, product, frequency) in list(processes.items()):
                outputs = {
                    station: rows(workdir / "result" / f"CODEX_FULLDAY_FLT_{product}_{frequency}_{station}.flt")
                    for station in ("HARB", "GODN")
                }
                state["cases"][name]["outputs"] = outputs
                state["cases"][name]["elapsed_seconds"] = round(time.monotonic() - state.setdefault("_started_mono", time.monotonic()), 3)
                code = proc.poll()
                if code is None:
                    continue
                stdout.close()
                stderr.close()
                state["cases"][name]["status"] = "COMPLETED" if code == 0 else "FAILED"
                state["cases"][name]["exit_code"] = code
                state["cases"][name]["finished_at"] = now()
                del processes[name]
            atomic_json(RUN_ROOT / "state.json", state)
            if processes:
                time.sleep(5)
    except KeyboardInterrupt:
        for proc, stdout, stderr, *_ in processes.values():
            proc.terminate()
        state["status"] = "INTERRUPTED"
        atomic_json(RUN_ROOT / "state.json", state)
        raise

    state.pop("_started_mono", None)
    state["status"] = "COMPLETED" if all(item["status"] == "COMPLETED" for item in state["cases"].values()) else "FAILED"
    state["finished_at"] = now()
    atomic_json(RUN_ROOT / "state.json", state)
    (RUN_ROOT / "report.md").write_text(
        "# FLT full-day regression\n\n"
        + "| Case | Exit | HARB rows | GODN rows |\n|---|---:|---:|---:|\n"
        + "\n".join(
            f"| {name} | {item.get('exit_code')} | {item.get('outputs', {}).get('HARB', {}).get('rows', 0)} | {item.get('outputs', {}).get('GODN', {}).get('rows', 0)} |"
            for name, item in state["cases"].items()
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
