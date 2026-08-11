# FGO experiment runner

`run_fgo_experiments.py` provides one entry point for a single experiment or a
parallel comparison matrix. It uses `analyze_fgo_solution.py` for accuracy and
convergence statistics instead of implementing a second statistical definition.

## Success criteria

An experiment passes only when all conditions below hold:

1. the process exits with code 0 before its timeout;
2. every configured station's `.fgo` file is changed by the current process;
3. the result has the exact XML-derived number of epochs, expected first/last
   SOW, a continuous sampling grid, valid rows, and finite coordinates;
4. stdout, stderr, and the newly written PPP log contain no configured fatal
   pattern;
5. when reference coordinates are supplied, the existing analyzer completes
   and writes its Markdown and JSON reports.

The final report includes process and matrix wall time, station epochs per
second, real-time factor, fixed rate, first sustained/permanent convergence,
first correct fixed epoch, and 3D RMS/MAE/P95/P99/maximum error.

The JSON analysis also includes status counts, baseline epoch matching, hourly
statistics, anomaly intervals, coordinate jumps, and optional quality gates.
The runner extracts `Spent` and DEBUG FGO phase timing when those lines are
present, including RAW/Ceres phase min/mean/max/total times, graph dimensions,
and average solver iterations.

## Single experiment

Paths passed on the command line are relative to the repository root. For the
current UPD/OSB sample data, the working directory is inferred from the XML
bias product. It can be set explicitly with `--work-dir NAME=DIR`.

```powershell
.\.venv\Scripts\python.exe scripts\run_fgo_experiments.py run `
  --case UPD_FF_NONE=build\parallel30\UPD_NONE.xml `
  --reference GODN=1130760.6931,-4831298.6759,3994155.1990 `
  --timeout-seconds 300
```

Multiple `--case` options form a small matrix. Use `--max-parallel` to bound
the number of independent GREAT processes; each process retains the Ceres
thread count from its XML.

## Reusable matrix

Copy `fgo_experiment_matrix.example.json`, then add one experiment per XML.
Configuration, executable, and working-directory paths in a manifest are
relative to the manifest file. `run_root` is relative to the repository root.

```powershell
.\.venv\Scripts\python.exe scripts\run_fgo_experiments.py run `
  --manifest scripts\fgo_experiment_matrix.example.json `
  --run-id df_ff_30m_01
```

The important manifest fields are:

```json
{
  "schema_version": 1,
  "run_root": "build/fgo_runs",
  "max_parallel": 4,
  "timeout_seconds": 900,
  "defaults": {"executable": "../build/Bin/Release/GREAT_PVT.exe"},
  "references": {"SITE": [1.0, 2.0, 3.0]},
  "analysis": {
    "window_minutes": 5,
    "horizontal_threshold": 0.1,
    "vertical_threshold": 0.2
  },
  "experiments": [
    {
      "name": "UPD_DF_NONE",
      "config": "path/to/unique-output-config.xml",
      "workdir": "../sample_data/PPPFLT_2023305",
      "tags": {"frequency_set": "DF"}
    }
  ]
}
```

Quality gates are opt-in so existing reports remain diagnostic-only. They can
be placed directly under `analysis` or under `analysis.quality_gates`:

```json
{
  "analysis": {
    "require_complete": true,
    "fail_on_gate": true,
    "max_sustained_minutes": 5.0,
    "max_permanent_minutes": 10.0,
    "min_correct_fixed_fraction": 0.95,
    "max_post_convergence_3d": 0.10,
    "max_coordinate_jumps": 0,
    "max_status_transitions": 2
  },
  "performance_gates": {
    "max_wall_seconds": 120.0,
    "min_realtime_factor": 10.0,
    "max_peak_rss_mb": 1500.0,
    "max_profile_avg_ms": {"raw.solve": 120.0},
    "max_profile_iterations": 25.0,
    "max_feedback_rejected": 0,
    "max_pseudo_inverse": 0
  }
}
```

Without `fail_on_gate`, quality-gate failures are written to the report but do
not change the process exit code. Performance gates are enforced whenever the
`performance_gates` block is present.

Each concurrently running configuration must use distinct output paths. The
runner intentionally does not delete or rename existing result files; it checks
that the executable overwrote them and rejects stale results.

Use `--dry-run` to validate the complete matrix without starting GREAT.

## Monitoring and artifacts

The run command prints station progress as `rows/expected`. A second terminal
can inspect the atomic state file while the matrix is running:

```powershell
.\.venv\Scripts\python.exe scripts\run_fgo_experiments.py status `
  --run-dir build\fgo_runs\df_ff_30m_01
```

Each run directory contains:

- `state.json`: current state and final machine-readable result;
- `manifest.resolved.json`: resolved configuration and provenance;
- `logs/`: separate stdout and stderr for every process;
- `analysis/`: detailed Markdown and JSON accuracy reports;
- `report.md`: matrix pass/fail, performance, convergence, and precision summary.

`Ctrl+C` terminates active child processes and records the matrix as
`INTERRUPTED`. A fixed `--run-id` is never silently reused, preventing results
from two executions from being mixed.
