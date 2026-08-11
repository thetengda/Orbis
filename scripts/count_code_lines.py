#!/usr/bin/env python3
"""Count source files and lines, grouped by the directory containing each file.

By default the script scans C/C++, Python, and CMake files.  ``nonblank`` is
reported as the practical code-line count; it includes comment lines but does
not include blank lines.  The physical line count is reported as ``lines``.

Examples::

    python scripts/count_code_lines.py
    python scripts/count_code_lines.py src --extensions .cpp,.h,.py
    python scripts/count_code_lines.py --json > code-line-count.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_EXTENSIONS = {
    ".c",
    ".cc",
    ".cpp",
    ".cxx",
    ".h",
    ".hh",
    ".hpp",
    ".hxx",
    ".py",
}
DEFAULT_FILENAMES = {"CMakeLists.txt"}
DEFAULT_EXCLUDED_DIRS = {
    ".git",
    ".venv",
    "__pycache__",
    "build",
    "doc",
    "sample_data",
}


@dataclass
class LineStats:
    """Statistics for one directory or for the complete scan."""

    files: int = 0
    lines: int = 0
    blank: int = 0
    nonblank: int = 0

    def add_file(self, lines: Iterable[str]) -> None:
        line_list = list(lines)
        self.files += 1
        self.lines += len(line_list)
        self.blank += sum(not line.strip() for line in line_list)
        self.nonblank += sum(bool(line.strip()) for line in line_list)

    def add(self, other: "LineStats") -> None:
        self.files += other.files
        self.lines += other.lines
        self.blank += other.blank
        self.nonblank += other.nonblank


def parse_csv(value: str, *, leading_dot: bool = False) -> set[str]:
    """Parse a comma-separated command-line option into a normalized set."""

    items = {item.strip() for item in value.split(",") if item.strip()}
    if leading_dot:
        items = {item if item.startswith(".") else f".{item}" for item in items}
        items = {item.lower() for item in items}
    return items


def iter_source_files(
    root: Path,
    extensions: set[str],
    filenames: set[str],
    excluded_dirs: set[str],
) -> Iterable[Path]:
    """Yield matching files while pruning excluded directories early."""

    # os.walk is avoided here so that this script only needs pathlib; the
    # directory tree is still pruned before descending into excluded dirs.
    pending = [root]
    while pending:
        current = pending.pop()
        try:
            children = list(current.iterdir())
        except OSError as exc:
            print(f"warning: cannot read {current}: {exc}", file=sys.stderr)
            continue
        for child in children:
            if child.is_dir():
                if child.name not in excluded_dirs and not child.is_symlink():
                    pending.append(child)
                continue
            if not child.is_file():
                continue
            if child.name in filenames or child.suffix.lower() in extensions:
                yield child


def count_lines(
    root: Path,
    *,
    extensions: set[str],
    filenames: set[str],
    excluded_dirs: set[str],
) -> tuple[dict[str, LineStats], LineStats]:
    """Return per-directory and total line statistics for ``root``."""

    by_directory: dict[str, LineStats] = defaultdict(LineStats)
    total = LineStats()
    for path in iter_source_files(root, extensions, filenames, excluded_dirs):
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError as exc:
            print(f"warning: cannot read {path}: {exc}", file=sys.stderr)
            continue
        relative_directory = path.parent.relative_to(root).as_posix() or "."
        by_directory[relative_directory].add_file(lines)
        total.add_file(lines)
    return dict(by_directory), total


def format_table(by_directory: dict[str, LineStats], total: LineStats) -> str:
    """Format the result as a readable fixed-width table."""

    rows = [(directory, stats) for directory, stats in sorted(by_directory.items())]
    rows.append(("TOTAL", total))
    headers = ("Directory", "Files", "Lines", "Blank", "Nonblank")
    values = [
        (directory, str(stats.files), str(stats.lines), str(stats.blank), str(stats.nonblank))
        for directory, stats in rows
    ]
    widths = [
        max(len(headers[index]), *(len(row[index]) for row in values))
        for index in range(len(headers))
    ]
    output = [
        "  ".join(header.ljust(width) for header, width in zip(headers, widths)),
        "  ".join("-" * width for width in widths),
    ]
    output.extend(
        "  ".join(value.rjust(width) if index else value.ljust(width) for index, (value, width) in enumerate(zip(row, widths)))
        for row in values
    )
    return "\n".join(output)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "root",
        nargs="?",
        type=Path,
        default=Path("."),
        help="directory to scan (default: current directory)",
    )
    parser.add_argument(
        "--extensions",
        default=",".join(sorted(DEFAULT_EXTENSIONS)),
        help="comma-separated extensions to include (default: C/C++ and Python)",
    )
    parser.add_argument(
        "--filenames",
        default=",".join(sorted(DEFAULT_FILENAMES)),
        help="comma-separated exact filenames to include (default: CMakeLists.txt)",
    )
    parser.add_argument(
        "--exclude-dir",
        action="append",
        default=[],
        metavar="NAME",
        help=(
            "directory name to skip; may be repeated "
            f"(defaults also include: {', '.join(sorted(DEFAULT_EXCLUDED_DIRS))})"
        ),
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="write machine-readable JSON instead of the table",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        raise SystemExit(f"error: scan root is not a directory: {root}")

    extensions = parse_csv(args.extensions, leading_dot=True)
    filenames = parse_csv(args.filenames)
    excluded_dirs = DEFAULT_EXCLUDED_DIRS | set(args.exclude_dir)
    by_directory, total = count_lines(
        root,
        extensions=extensions,
        filenames=filenames,
        excluded_dirs=excluded_dirs,
    )

    if args.json:
        result = {
            "root": str(root),
            "extensions": sorted(extensions),
            "filenames": sorted(filenames),
            "excluded_dirs": sorted(excluded_dirs),
            "directories": {name: asdict(stats) for name, stats in sorted(by_directory.items())},
            "total": asdict(total),
        }
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print(f"Root: {root}")
        print("Nonblank is used as the practical code-line count; comments are included.")
        print(format_table(by_directory, total))


if __name__ == "__main__":
    main()
