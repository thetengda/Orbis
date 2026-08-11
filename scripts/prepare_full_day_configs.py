#!/usr/bin/env python3
"""Prepare isolated full-day FGO regression configurations.

The checked-in sample XML files are also used by historical runs, so this
helper derives test-only copies under ``build/full_day_tests_20260811``.  It
does not modify sample data or the checked-in XML files.
"""

from __future__ import annotations

import copy
import json
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = ROOT / "build" / "full_day_tests_20260811"
CONFIG_ROOT = RUN_ROOT / "configs"

BASES = {
    ("UPD", "DF"): ROOT / "sample_data" / "PPPFLT_2023305" / "xml" / "GREAT_PPPFGO_kinematic_DF_Fixed.xml",
    ("UPD", "FF"): ROOT / "sample_data" / "PPPFLT_2023305" / "xml" / "GREAT_PPPFLT_kinematic_FF_Fixed.xml",
    ("OSB", "DF"): ROOT / "sample_data" / "PPPFLT_2023305_OSB" / "xml" / "GREAT_PPPFGO_kinematic_DF_Fixed.xml",
    ("OSB", "FF"): ROOT / "sample_data" / "PPPFLT_2023305_OSB" / "xml" / "GREAT_PPPFLT_kinematic_FF_Fixed.xml",
}

FLT_BASES = {
    ("UPD", "DF"): ROOT / "sample_data" / "PPPFLT_2023305" / "xml" / "GREAT_PPPFLT_kinematic_DF_Fixed.xml",
    ("UPD", "FF"): ROOT / "sample_data" / "PPPFLT_2023305" / "xml" / "GREAT_PPPFLT_kinematic_FF_Fixed.xml",
    ("OSB", "DF"): ROOT / "sample_data" / "PPPFLT_2023305_OSB" / "xml" / "GREAT_PPPFLT_kinematic_DF_Fixed.xml",
    ("OSB", "FF"): ROOT / "sample_data" / "PPPFLT_2023305_OSB" / "xml" / "GREAT_PPPFLT_kinematic_FF_Fixed.xml",
}


def child(parent: ET.Element, name: str, value: str) -> ET.Element:
    node = parent.find(name)
    if node is None:
        node = ET.SubElement(parent, name)
    node.text = value
    return node


def ensure_fgo(root: ET.Element) -> ET.Element:
    node = root.find("fgo")
    if node is None:
        node = ET.Element("fgo")
        gen = root.find("gen")
        insert_at = list(root).index(gen) + 1 if gen is not None else 0
        root.insert(insert_at, node)
    child(node, "gnss_enable", "1")
    child(node, "amb_propagation", "1")
    child(node, "gnss_window_size", "2")
    child(node, "gnss_num_threads", "4")
    return node


def set_output(root: ET.Element, name: str) -> None:
    outputs = root.find("outputs")
    if outputs is None:
        outputs = ET.SubElement(root, "outputs", {"append": "false", "verb": "0"})
    log = outputs.find("log")
    if log is None:
        log = ET.SubElement(outputs, "log")
    log.set("type", "BASIC")
    log.set("name", f"FULLDAY_{name}.log")
    log.set("level", "INFO")
    child(outputs, "ppp", f"./result/CODEX_FULLDAY_{name}_$(rec).ppp.log")
    child(outputs, "flt", f"./result/CODEX_FULLDAY_{name}_$(rec).flt")
    child(outputs, "fgo", f"./result/CODEX_FULLDAY_{name}_$(rec).fgo")
    child(outputs, "fgo_ar", f"./result/CODEX_FULLDAY_{name}_$(rec).fgo.ar")
    child(outputs, "ratio", f"./result/CODEX_FULLDAY_{name}_$(rec).ratio")


def make_config(product: str, frequency_set: str, mode: str,
                ambiguity_fix_factor_enable: bool) -> Path:
    base = BASES[(product, frequency_set)]
    root = ET.parse(base).getroot()
    gen = root.find("gen")
    if gen is None:
        raise RuntimeError(f"{base} has no <gen>")
    child(gen, "beg", "2023-11-01 00:00:00")
    child(gen, "end", "2023-11-01 23:59:30")
    child(gen, "int", "30")
    child(gen, "rec", "HARB GODN")
    child(gen, "est", "FGO")

    fgo = ensure_fgo(root)
    child(fgo, "ambiguity_fix_factor_enable",
          "true" if ambiguity_fix_factor_enable else "false")

    process = root.find("process")
    if process is None:
        process = ET.SubElement(root, "process")
    child(process, "frequency", "2" if frequency_set == "DF" else "5")

    ambiguity = root.find("ambiguity")
    if ambiguity is None:
        ambiguity = ET.SubElement(root, "ambiguity")
    child(ambiguity, "fix_mode", "SEARCH")
    child(ambiguity, "upd_mode", product)

    name = f"{product}_{frequency_set}_{mode}"
    set_output(root, name)
    target = CONFIG_ROOT / f"{name}.xml"
    target.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root, space="    ")
    ET.ElementTree(root).write(target, encoding="utf-8", xml_declaration=True)
    return target


def make_flt_config(product: str, frequency_set: str) -> Path:
    base = FLT_BASES[(product, frequency_set)]
    root = ET.parse(base).getroot()
    gen = root.find("gen")
    if gen is None:
        raise RuntimeError(f"{base} has no <gen>")
    child(gen, "beg", "2023-11-01 00:00:00")
    child(gen, "end", "2023-11-01 23:59:30")
    child(gen, "int", "30")
    child(gen, "rec", "HARB GODN")
    child(gen, "est", "FLT")

    process = root.find("process")
    if process is None:
        process = ET.SubElement(root, "process")
    child(process, "frequency", "2" if frequency_set == "DF" else "5")

    ambiguity = root.find("ambiguity")
    if ambiguity is None:
        ambiguity = ET.SubElement(root, "ambiguity")
    child(ambiguity, "fix_mode", "SEARCH")
    child(ambiguity, "upd_mode", product)

    outputs = root.find("outputs")
    if outputs is None:
        outputs = ET.SubElement(root, "outputs", {"append": "false", "verb": "0"})
    log = outputs.find("log")
    if log is None:
        log = ET.SubElement(outputs, "log")
    log.set("type", "BASIC")
    log.set("name", f"FULLDAY_FLT_{product}_{frequency_set}.log")
    log.set("level", "INFO")
    child(outputs, "ppp", f"./result/CODEX_FULLDAY_FLT_{product}_{frequency_set}_$(rec).ppp.log")
    child(outputs, "flt", f"./result/CODEX_FULLDAY_FLT_{product}_{frequency_set}_$(rec).flt")
    child(outputs, "ratio", f"./result/CODEX_FULLDAY_FLT_{product}_{frequency_set}_$(rec).ratio")
    # A stale FGO node in a source template is harmless for FLT, but removing
    # it prevents accidental FGO output creation if the parser is extended.
    for node_name in ("fgo", "fgo_ar"):
        node = outputs.find(node_name)
        if node is not None:
            outputs.remove(node)

    target = CONFIG_ROOT / f"FLT_{product}_{frequency_set}.xml"
    target.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(root, space="    ")
    ET.ElementTree(root).write(target, encoding="utf-8", xml_declaration=True)
    return target


def main() -> None:
    experiments = []
    for product in ("UPD", "OSB"):
        workdir = "../../sample_data/PPPFLT_2023305_OSB" if product == "OSB" else "../../sample_data/PPPFLT_2023305"
        for frequency_set in ("DF", "FF"):
            for mode, factor_enabled in (("PARAMETER", False), ("CONSTRAINT", True)):
                config = make_config(
                    product, frequency_set, mode, factor_enabled)
                experiments.append(
                    {
                        "name": f"{product}_{frequency_set}_{mode}",
                        "config": f"configs/{config.name}",
                        "workdir": workdir,
                        "tags": {"product": product, "frequency_set": frequency_set, "feedback": mode},
                    }
                )

    manifest = {
        "schema_version": 1,
        "run_root": "build/fgo_runs",
        "max_parallel": 4,
        "timeout_seconds": 21600,
        "poll_seconds": 2,
        "status_seconds": 60,
        "defaults": {"executable": "../../build/Bin/Release/GREAT_PVT.exe"},
        "references": {
            "GODN": [1130760.6931, -4831298.6759, 3994155.1990],
            "HARB": [5084657.6078, 2670325.4787, -2768480.8416],
        },
        "analysis": {
            "window_minutes": 30,
            "horizontal_threshold": 0.10,
            "vertical_threshold": 0.20,
            "require_complete": True,
        },
        "experiments": experiments,
    }
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    (RUN_ROOT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    flt_experiments = []
    for product in ("UPD", "OSB"):
        workdir = "../../sample_data/PPPFLT_2023305_OSB" if product == "OSB" else "../../sample_data/PPPFLT_2023305"
        for frequency_set in ("DF", "FF"):
            config = make_flt_config(product, frequency_set)
            flt_experiments.append(
                {
                    "name": f"FLT_{product}_{frequency_set}",
                    "config": f"configs/{config.name}",
                    "workdir": workdir,
                    "product": product,
                    "frequency_set": frequency_set,
                }
            )
    flt_manifest = {
        "schema_version": 1,
        "executable": "../../build/Bin/Release/GREAT_PVT.exe",
        "run_root": "build/full_day_tests_20260811/flt_runs",
        "experiments": flt_experiments,
    }
    (RUN_ROOT / "flt_manifest.json").write_text(json.dumps(flt_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(RUN_ROOT / "manifest.json")


if __name__ == "__main__":
    main()
