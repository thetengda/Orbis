#!/usr/bin/env python3
"""Run the OSB half of the full-day FGO matrix independently of the runner."""

from __future__ import annotations

import datetime as dt
import json
import math
import re
import subprocess
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "build" / "full_day_tests_20260811" / "osb_fgo_runs"
CONFIG_ROOT = ROOT / "build" / "full_day_tests_20260811" / "configs"
EXE = ROOT / "build" / "Bin" / "Release" / "GREAT_PVT.exe"
WORKDIR = ROOT / "sample_data" / "PPPFLT_2023305_OSB"
CASES = ("OSB_DF_PARAMETER", "OSB_DF_CONSTRAINT", "OSB_FF_PARAMETER", "OSB_FF_CONSTRAINT")


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).astimezone().isoformat()


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def count_rows(path: Path) -> dict:
    out = {"exists": path.is_file(), "rows": 0, "malformed": 0, "nonfinite": 0, "first_sow": None, "last_sow": None}
    if not path.is_file():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if not fields or not re.fullmatch(r"[-+]?\d+(?:\.\d+)?", fields[0]):
            continue
        if len(fields) < 4:
            out["malformed"] += 1
            continue
        try:
            values = [float(item) for item in fields[:4]]
        except ValueError:
            out["malformed"] += 1
            continue
        if not all(math.isfinite(item) for item in values):
            out["nonfinite"] += 1
            continue
        out["rows"] += 1
        out["first_sow"] = values[0] if out["first_sow"] is None else out["first_sow"]
        out["last_sow"] = values[0]
    return out


def main() -> None:
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    log_root = RUN_ROOT / "logs"
    log_root.mkdir(parents=True, exist_ok=True)
    state = {"status": "RUNNING", "started_at": now(), "executable": str(EXE), "cases": {}}
    processes = {}
    for name in CASES:
        config = CONFIG_ROOT / f"{name}.xml"
        stdout = (log_root / f"{name}.stdout.log").open("w", encoding="utf-8")
        stderr = (log_root / f"{name}.stderr.log").open("w", encoding="utf-8")
        proc = subprocess.Popen([str(EXE), "-x", str(config)], cwd=WORKDIR, stdout=stdout, stderr=stderr)
        product, frequency, mode = name.split("_")
        processes[name] = (proc, stdout, stderr, frequency)
        state["cases"][name] = {"status": "RUNNING", "pid": proc.pid, "config": str(config), "started_at": now(), "outputs": {}}
    atomic_json(RUN_ROOT / "state.json", state)
    started_mono = time.monotonic()
    while processes:
        for name, (proc, stdout, stderr, frequency) in list(processes.items()):
            outputs = {
                station: count_rows(WORKDIR / "result" / f"CODEX_FULLDAY_OSB_{frequency}_{name.split('_')[-1]}_{station}.fgo")
                for station in ("HARB", "GODN")
            }
            state["cases"][name]["outputs"] = outputs
            state["cases"][name]["elapsed_seconds"] = round(time.monotonic() - started_mono, 3)
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
    state["status"] = "COMPLETED" if all(item["status"] == "COMPLETED" for item in state["cases"].values()) else "FAILED"
    state["finished_at"] = now()
    atomic_json(RUN_ROOT / "state.json", state)


if __name__ == "__main__":
    main()
