#!/usr/bin/env python3
"""
Batch GEPIA correlation PDF runner for HEM/chromatin/control pipeline.

This uses the GEPIA Python package API pattern confirmed locally:
    c = gepia.correlation()
    c.setParams({...})
    c.query()  # returns a PDF path

It saves/renames each PDF systematically and writes an index CSV.

Example:
    python run_gepia_pdf_batch.py --phase 1 --dataset BRCA_Tumor BLCA_Tumor --sleep 8
    python run_gepia_pdf_batch.py --all --dataset BRCA_Tumor BLCA_Tumor --sleep 10
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import string
import sys
import time
from pathlib import Path
from typing import Iterable, List, Tuple, Dict

import gepia

HEM_GENES = ["ALAS1", "ALAD", "HMBS", "UROS", "UROD", "CPOX", "PPOX", "FECH"]
CHROMATIN_REGULATORS = ["H2AFZ", "ATAD2", "BRD4", "BRD2"]
PROLIFERATION_CONTROLS = ["MKI67", "PCNA"]
GENERAL_CHROMATIN_CONTROLS = ["HDAC1"]
SPECIFICITY_CONTROLS = ["FASN", "ACACA", "TFRC", "FTH1", "ACTB", "GAPDH"]

CONSECUTIVE_HEM_PAIRS = [
    ("ALAS1", "ALAD"),
    ("ALAD", "HMBS"),
    ("HMBS", "UROS"),
    ("UROS", "UROD"),
    ("UROD", "CPOX"),
    ("CPOX", "PPOX"),
    ("PPOX", "FECH"),
]
NONCONSECUTIVE_HEM_PAIRS = [
    ("HMBS", "FECH"),
    ("ALAS1", "FECH"),
    ("HMBS", "PPOX"),
]
POSITIVE_CONTROLS = [
    ("E2F1", "CCNE1"),
    ("CDK2", "CCNE1"),
]

# Common TCGA tumor dataset labels in GEPIA package style. Edit if needed.
DEFAULT_CANCER_DATASETS = [
    "ACC_Tumor", "BLCA_Tumor", "BRCA_Tumor", "CESC_Tumor", "CHOL_Tumor", "COAD_Tumor",
    "DLBC_Tumor", "ESCA_Tumor", "GBM_Tumor", "HNSC_Tumor", "KICH_Tumor", "KIRC_Tumor",
    "KIRP_Tumor", "LAML_Tumor", "LGG_Tumor", "LIHC_Tumor", "LUAD_Tumor", "LUSC_Tumor",
    "MESO_Tumor", "OV_Tumor", "PAAD_Tumor", "PCPG_Tumor", "PRAD_Tumor", "READ_Tumor",
    "SARC_Tumor", "SKCM_Tumor", "STAD_Tumor", "TGCT_Tumor", "THCA_Tumor", "THYM_Tumor",
    "UCEC_Tumor", "UCS_Tumor", "UVM_Tumor",
]


def safe_name(s: str) -> str:
    keep = string.ascii_letters + string.digits + "_-"
    return "".join(ch if ch in keep else "_" for ch in s)


def build_tasks(phases: Iterable[str], per_cancer: bool = False, cancers: List[str] | None = None) -> List[Dict[str, str]]:
    tasks: List[Dict[str, str]] = []
    phases = set(str(p) for p in phases)

    if "1" in phases:
        for g1, g2 in CONSECUTIVE_HEM_PAIRS + NONCONSECUTIVE_HEM_PAIRS:
            tasks.append({"phase": "1A", "analysis": "normal_or_selected_within_heme", "gene1": g1, "gene2": g2})
        for reg in CHROMATIN_REGULATORS:
            for hem in HEM_GENES:
                tasks.append({"phase": "1B", "analysis": "normal_or_selected_chromatin_vs_heme", "gene1": reg, "gene2": hem})
        for ctrl in PROLIFERATION_CONTROLS + GENERAL_CHROMATIN_CONTROLS:
            for hem in HEM_GENES:
                tasks.append({"phase": "1C", "analysis": "confound_vs_heme", "gene1": ctrl, "gene2": hem})
        for ctrl in SPECIFICITY_CONTROLS:
            tasks.append({"phase": "1C", "analysis": "h2afz_specificity", "gene1": "H2AFZ", "gene2": ctrl})

    if "2" in phases:
        for g1, g2 in CONSECUTIVE_HEM_PAIRS + NONCONSECUTIVE_HEM_PAIRS:
            tasks.append({"phase": "2A", "analysis": "pancancer_within_heme", "gene1": g1, "gene2": g2})
        for reg in CHROMATIN_REGULATORS:
            for hem in HEM_GENES:
                tasks.append({"phase": "2B", "analysis": "pancancer_chromatin_vs_heme", "gene1": reg, "gene2": hem})
        for ctrl in PROLIFERATION_CONTROLS + GENERAL_CHROMATIN_CONTROLS:
            for hem in HEM_GENES:
                tasks.append({"phase": "2C", "analysis": "pancancer_confound_vs_heme", "gene1": ctrl, "gene2": hem})
        for ctrl in SPECIFICITY_CONTROLS:
            tasks.append({"phase": "2C", "analysis": "pancancer_h2afz_specificity", "gene1": "H2AFZ", "gene2": ctrl})
        for g1, g2 in POSITIVE_CONTROLS:
            tasks.append({"phase": "2D", "analysis": "positive_control", "gene1": g1, "gene2": g2})

    if "3" in phases:
        selected = cancers or DEFAULT_CANCER_DATASETS
        for ds in selected:
            tasks.append({"phase": "3A", "analysis": "per_cancer_h2afz_hmbs", "gene1": "H2AFZ", "gene2": "HMBS", "dataset_override": ds})
            tasks.append({"phase": "3B", "analysis": "per_cancer_hmbs_uros", "gene1": "HMBS", "gene2": "UROS", "dataset_override": ds})

    if "5" in phases:
        selected = cancers or ["GBM_Tumor", "BLCA_Tumor", "ESCA_Tumor", "HNSC_Tumor"]
        for ds in selected:
            for reg in CHROMATIN_REGULATORS:
                for hem in HEM_GENES:
                    tasks.append({"phase": "5A", "analysis": "focused_chromatin_vs_heme", "gene1": reg, "gene2": hem, "dataset_override": ds})
            for g1, g2 in CONSECUTIVE_HEM_PAIRS:
                tasks.append({"phase": "5B", "analysis": "focused_within_heme", "gene1": g1, "gene2": g2, "dataset_override": ds})
            for ctrl in ["MKI67", "HDAC1"]:
                tasks.append({"phase": "5C", "analysis": "focused_confound_hmbs", "gene1": ctrl, "gene2": "HMBS", "dataset_override": ds})

    return tasks


def run_one(gene1: str, gene2: str, dataset: List[str], method: str, out_dir: Path, prefix: str) -> Tuple[str, str]:
    c = gepia.correlation()
    c.setOutDir(str(out_dir))
    c.setParams({
        "dataset": dataset,
        "methodoption": method,
        "signature1": [gene1],
        "signature1_norm": "",
        "signature2": [gene2],
        "signature2_norm": "",
    })
    result = c.query()
    if result is None:
        return "failed", "GEPIA returned None. Check dataset labels/parameters."

    src = Path(str(result))
    if not src.is_absolute():
        # GEPIA may return ./file.pdf relative to cwd or out_dir.
        cwd_candidate = Path.cwd() / src
        out_candidate = out_dir / src.name
        if cwd_candidate.exists():
            src = cwd_candidate
        elif out_candidate.exists():
            src = out_candidate

    if not src.exists():
        return "unknown", f"GEPIA returned {result}, but file was not found by runner."

    dest = out_dir / f"{safe_name(prefix)}__{gene1}_vs_{gene2}.pdf"
    if src.resolve() != dest.resolve():
        shutil.move(str(src), str(dest))
    return "ok", str(dest)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", action="append", choices=["1", "2", "3", "5"], help="Phase to run. Repeatable. Example: --phase 2 --phase 3")
    ap.add_argument("--all", action="store_true", help="Run phases 1, 2, 3, and 5")
    ap.add_argument("--dataset", nargs="+", default=["BRCA_Tumor", "BLCA_Tumor"], help="GEPIA dataset labels for pooled analyses")
    ap.add_argument("--cancers", nargs="*", default=None, help="Cancer datasets for per-cancer/focused phases, e.g. GBM_Tumor BLCA_Tumor")
    ap.add_argument("--method", default="pearson", choices=["pearson", "spearman"], help="Correlation method")
    ap.add_argument("--out", default="gepia_pdf_results", help="Output directory")
    ap.add_argument("--sleep", type=float, default=8.0, help="Seconds between GEPIA requests")
    ap.add_argument("--limit", type=int, default=None, help="Limit number of tasks for testing")
    ap.add_argument("--resume", action="store_true", help="Skip tasks whose output PDF already exists")
    args = ap.parse_args()

    phases = ["1", "2", "3", "5"] if args.all else (args.phase or ["2"])
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    tasks = build_tasks(phases, cancers=args.cancers)
    if args.limit:
        tasks = tasks[: args.limit]

    index_path = out_dir / "gepia_pdf_index.csv"
    write_header = not index_path.exists()

    print(f"Running {len(tasks)} GEPIA PDF queries")
    print(f"Default pooled dataset: {args.dataset}")
    print(f"Output: {out_dir.resolve()}")
    print("WARNING: GEPIA rate-limits. Use --sleep 8 or higher for large runs.\n")

    with open(index_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["phase", "analysis", "gene1", "gene2", "dataset", "method", "status", "path_or_message"])
        if write_header:
            writer.writeheader()

        for i, task in enumerate(tasks, 1):
            gene1 = task["gene1"]
            gene2 = task["gene2"]
            ds = [task["dataset_override"]] if "dataset_override" in task else args.dataset
            prefix = f"{task['phase']}__{task['analysis']}__{'_'.join(ds)}"
            expected = out_dir / f"{safe_name(prefix)}__{gene1}_vs_{gene2}.pdf"
            if args.resume and expected.exists():
                status, msg = "skipped", str(expected)
            else:
                print(f"[{i}/{len(tasks)}] {task['phase']} {task['analysis']}: {gene1} vs {gene2} on {ds}")
                try:
                    status, msg = run_one(gene1, gene2, ds, args.method, out_dir, prefix)
                except Exception as e:
                    status, msg = "error", repr(e)
                print(f"    {status}: {msg}")

            writer.writerow({
                "phase": task["phase"],
                "analysis": task["analysis"],
                "gene1": gene1,
                "gene2": gene2,
                "dataset": ";".join(ds),
                "method": args.method,
                "status": status,
                "path_or_message": msg,
            })
            fh.flush()

            if i < len(tasks):
                time.sleep(args.sleep + random.uniform(0, 2.0))

    print(f"\nDone. Index: {index_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
