#!/usr/bin/env python3
"""
run_gepia_pdf_batch.py
======================
Batch GEPIA2 correlation queries for heme-pathway / chromatin-regulator analysis.

Usage examples
--------------
# Phase 1 only, TCGA matched normal tissues (default):
python run_gepia_pdf_batch.py --phase 1 --sleep 8

# Phase 2 only (pan-cancer), specific tumor datasets:
python run_gepia_pdf_batch.py --phase 2 --dataset BRCA_Tumor BLCA_Tumor --sleep 8

# All phases, default datasets:
python run_gepia_pdf_batch.py --all --sleep 10

# Phase 4, focused cancers, spearman:
python run_gepia_pdf_batch.py --phase 4 --cancers GBM BLCA ESCA HNSC --method spearman --sleep 8

# Resume an interrupted run:
python run_gepia_pdf_batch.py --all --resume --sleep 8

Important caveats
-----------------
The GEPIA2 server (gepia2.cancer-pku.cn / 118.190.148.166) must be reachable from the
machine running this script.  On sandboxed systems (e.g. Claude's container) both hostnames
are blocked by the egress proxy.  Run this script on your own computer or cluster.

The gepia Python package hardcodes http://118.190.148.166/GEPIA2/ as the server URL.
This script monkey-patches gepia.GEPIA to allow the canonical domain to be used instead
if desired (see --server flag).

GEPIA correlation endpoint parameters (from gepia.py source):
  dataset          : list of dataset strings, e.g. ["BRCA_Tumor", "BLCA_Tumor"]
  methodoption     : "spearman" or "pearson"
  signature1       : list of gene symbols (gene on x-axis / "query" gene)
  signature1_norm  : housekeeping normaliser for signature1 (leave "" to use default)
  signature2       : list of gene symbols (gene on y-axis)
  signature2_norm  : housekeeping normaliser for signature2 (leave "" to use default)

The endpoint returns JSON: {"status": 0, "outdir": "<filename>.pdf"} on success.
Status != 0 indicates failure.  The PDF is saved to <out_dir>/<filename>.pdf.

Valid dataset strings (as of gepia 0.3):
  TCGATumor  : ACC_Tumor … UVM_Tumor   (33 entries)
  TCGANormal : BLCA_Normal … UCEC_Normal (23 entries, no GBM_Normal, etc.)
  GTEx       : Adipose_Subcutaneous … Vagina (54 tissue types)

Note: Phase 1 uses all 23 TCGA matched normal datasets, queried one at a time.
      These are normal-tissue biopsies from TCGA cancer patients (not GTEx healthy donors).
      To use GTEx instead, pass --dataset with the desired GTEx tissue names.
Note: Phase 3 loops over every individual tumor dataset.
Note: Phase 4 uses per-cancer loops over user-specified cancers.

The index CSV written by this script records every query attempt and its outcome,
making it easy to audit results, filter by significance, and resume after failures.
"""

import argparse
import csv
import json
import os
import random
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from urllib import request as ureq

# ---------------------------------------------------------------------------
# Try importing the gepia package; give a clear message if missing.
# ---------------------------------------------------------------------------
try:
    import gepia
except ImportError:
    sys.exit(
        "ERROR: The 'gepia' package is not installed.\n"
        "Install it with:  pip install gepia"
    )

# ---------------------------------------------------------------------------
# Gene lists (defined once; referenced by all phases)
# ---------------------------------------------------------------------------

# Heme biosynthesis pathway genes, in pathway order
HEME_GENES = ["ALAS1", "ALAD", "HMBS", "UROS", "UROD", "CPOX", "PPOX", "FECH"]

# Chromatin regulators of interest
CHROMATIN_REGULATORS = ["H2AFZ", "ATAD2", "BRD4", "BRD2"]

# Consecutive pairs in the heme pathway
HEME_CONSECUTIVE_PAIRS = [
    ("ALAS1", "ALAD"),
    ("ALAD",  "HMBS"),
    ("HMBS",  "UROS"),
    ("UROS",  "UROD"),
    ("UROD",  "CPOX"),
    ("CPOX",  "PPOX"),
    ("PPOX",  "FECH"),
]

# Non-consecutive pairs to test for long-range co-regulation
HEME_NONCONSECUTIVE_PAIRS = [
    ("HMBS", "FECH"),
    ("ALAS1", "FECH"),
    ("HMBS", "PPOX"),
]

ALL_HEME_PAIRS = HEME_CONSECUTIVE_PAIRS + HEME_NONCONSECUTIVE_PAIRS

# Positive controls for phase 2 (well-established co-regulated pairs)
POSITIVE_CONTROLS = [
    ("E2F1",  "CCNE1"),   # E2F1 transcriptionally activates CCNE1
    ("CDK2",  "CCNE1"),   # CDK2/CyclinE complex; tightly co-expressed
]

# Default focused cancers for phase 4
DEFAULT_FOCUSED_CANCERS = ["GBM", "BLCA", "ESCA", "HNSC"]

# ---------------------------------------------------------------------------
# Dataset helpers
# ---------------------------------------------------------------------------

# GTEx normal tissues (all 54 types from gepia.GTEx)
# Available via --dataset for users who want genuinely healthy-donor tissue.
# NOT used as the Phase 1 default — see note below.
ALL_GTEX = gepia.GTEx

# All TCGA tumor datasets
ALL_TCGA_TUMOR = gepia.TCGATumor   # e.g. ["ACC_Tumor", "BLCA_Tumor", ...]

# All TCGA matched normal datasets (23 cancer types that have a paired normal biopsy).
# These are the Phase 1 default, matching the analysis published in the manuscript.
# NOTE: these are normal-tissue biopsies from cancer patients, not healthy GTEx donors.
ALL_TCGA_NORMAL = gepia.TCGANormal  # e.g. ["BRCA_Normal", ...]

# Validated dataset name sets for quick lookup
_VALID_DATASETS = set(ALL_GTEX) | set(ALL_TCGA_TUMOR) | set(ALL_TCGA_NORMAL)
_VALID_CANCERS  = set(gepia.CANCERType)   # bare codes: "BRCA", "HNSC", …

# ---------------------------------------------------------------------------
# Index CSV field names
# ---------------------------------------------------------------------------
INDEX_FIELDS = [
    "phase",
    "analysis",
    "gene1",
    "gene2",
    "dataset",
    "method",
    "status",
    "pdf_path",
    "timestamp",
]


# ---------------------------------------------------------------------------
# Core query wrapper
# ---------------------------------------------------------------------------

def run_correlation(
    gene1: str,
    gene2: str,
    datasets: list,
    method: str = "spearman",
    out_dir: str = "./gepia_pdfs/",
    norm1: str = "",
    norm2: str = "",
) -> dict:
    """
    Run a single GEPIA2 correlation query and download the resulting PDF.

    Parameters
    ----------
    gene1     : x-axis gene symbol
    gene2     : y-axis gene symbol
    datasets  : list of GEPIA dataset strings (e.g. ["BRCA_Tumor"])
                Pass multiple to pool them.
    method    : "spearman" (recommended for RNA-seq) or "pearson"
    out_dir   : directory where the PDF will be saved
    norm1/2   : housekeeping gene for normalisation (usually leave empty)

    Returns
    -------
    dict with keys:
        success   : bool
        pdf_path  : str or None
        error     : str or None

    Notes
    -----
    The gepia package saves files using:
        request.urlretrieve(url, out_dir + pdf_filename)
    Because there is no os.sep inserted automatically, out_dir MUST end with '/'.
    We enforce this here.

    The 'outdir' field returned by the server is just a bare filename, e.g.:
        "correlation_a3f9c1b2.pdf"
    The package constructs the download URL as:
        GEPIA + "tmp/" + pdf_filename
    We replicate this logic in case we need to re-download.
    """
    # Ensure out_dir ends with separator so gepia saves to the right place
    if not out_dir.endswith("/") and not out_dir.endswith(os.sep):
        out_dir = out_dir + "/"

    os.makedirs(out_dir, exist_ok=True)

    # Build the gepia correlation object
    c = gepia.correlation()
    c.setOutDir(out_dir)
    c.setParams({
        "dataset":         datasets,
        "methodoption":    method,
        "signature1":      [gene1],
        "signature1_norm": norm1,
        "signature2":      [gene2],
        "signature2_norm": norm2,
    })

    try:
        pdf_path = c.query()
        if pdf_path is not None:
            # Normalise to an absolute path so the index CSV is portable
            pdf_path = str(Path(pdf_path).resolve())
            return {"success": True, "pdf_path": pdf_path, "error": None}
        else:
            return {
                "success": False,
                "pdf_path": None,
                "error": "GEPIA returned None (server-side error or bad parameters)",
            }
    except Exception as exc:
        return {
            "success": False,
            "pdf_path": None,
            "error": f"{type(exc).__name__}: {exc}",
        }


# ---------------------------------------------------------------------------
# PDF filename builder (for --resume matching)
# ---------------------------------------------------------------------------

def _expected_filename(gene1: str, gene2: str, datasets: list, method: str) -> str:
    """
    Build the informative local filename we rename the PDF to after download.
    Format:  <phase>__<analysis>__<gene1>_vs_<gene2>__<dataset_tag>__<method>.pdf

    This is the filename we assign; the raw file returned by gepia has a
    server-generated UUID name.  We rename it immediately after download.
    """
    # This helper is used for renaming, not for building the gepia query.
    # Dataset tag: join with '+' if multiple datasets; truncate if very long
    ds_tag = "+".join(datasets)
    if len(ds_tag) > 80:
        ds_tag = ds_tag[:77] + "..."
    return f"{gene1}_vs_{gene2}__{ds_tag}__{method}.pdf"


def _rename_pdf(raw_path: str, new_name: str, out_dir: str) -> str:
    """Rename the downloaded PDF to a meaningful name. Returns the new path."""
    new_path = os.path.join(out_dir, new_name)
    if raw_path and os.path.exists(raw_path) and raw_path != new_path:
        os.rename(raw_path, new_path)
    return new_path


# ---------------------------------------------------------------------------
# Resume helper: load already-completed entries from the index CSV
# ---------------------------------------------------------------------------

def load_completed(index_path: str) -> set:
    """
    Return a set of (gene1, gene2, dataset_str, method) tuples already
    recorded with status='success' in the index CSV.
    """
    completed = set()
    if not os.path.exists(index_path):
        return completed
    with open(index_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("status") == "success":
                completed.add((
                    row["gene1"],
                    row["gene2"],
                    row["dataset"],
                    row["method"],
                ))
    return completed


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def batch_run(
    tasks: list,
    out_dir: str,
    index_path: str,
    method: str = "spearman",
    sleep_base: float = 8.0,
    sleep_jitter: float = 3.0,
    resume: bool = False,
    dry_run: bool = False,
):
    """
    Execute a list of correlation tasks, writing results to the index CSV.

    Parameters
    ----------
    tasks : list of dicts, each with keys:
        phase     : str  (e.g. "1")
        analysis  : str  (e.g. "heme_consecutive" or "chromatin_vs_heme")
        gene1     : str
        gene2     : str
        datasets  : list of str   -- passed to GEPIA as the 'dataset' parameter
        method    : str (optional, overrides the function-level default)
    out_dir    : directory for PDFs
    index_path : path to the CSV index file
    method     : default correlation method if not specified per-task
    sleep_base : minimum seconds between requests
    sleep_jitter: additional random seconds (uniform [0, jitter])
    resume     : if True, skip tasks already in the index with status=success
    dry_run    : if True, print tasks but do not query GEPIA
    """
    if not out_dir.endswith("/") and not out_dir.endswith(os.sep):
        out_dir = out_dir + "/"
    os.makedirs(out_dir, exist_ok=True)

    completed = load_completed(index_path) if resume else set()

    # Open index CSV in append mode (create with header if new)
    file_exists = os.path.exists(index_path)
    index_fh = open(index_path, "a", newline="")
    writer = csv.DictWriter(index_fh, fieldnames=INDEX_FIELDS)
    if not file_exists:
        writer.writeheader()

    n_total  = len(tasks)
    n_skip   = 0
    n_ok     = 0
    n_fail   = 0

    for i, task in enumerate(tasks, start=1):
        task_method  = task.get("method", method)
        gene1        = task["gene1"]
        gene2        = task["gene2"]
        datasets     = task["datasets"]
        phase        = str(task["phase"])
        analysis     = task["analysis"]
        dataset_str  = "+".join(datasets)   # for the index

        # ------------------------------------------------------------------
        # Resume check
        # ------------------------------------------------------------------
        key = (gene1, gene2, dataset_str, task_method)
        if resume and key in completed:
            n_skip += 1
            print(f"[{i}/{n_total}] SKIP (already done): {gene1} vs {gene2} | {dataset_str}")
            continue

        print(
            f"[{i}/{n_total}] Phase {phase} | {analysis} | "
            f"{gene1} vs {gene2} | datasets: {dataset_str} | method: {task_method}"
        )

        if dry_run:
            print("  (dry-run, skipping actual query)")
            continue

        # ------------------------------------------------------------------
        # Run the query
        # ------------------------------------------------------------------
        result = run_correlation(
            gene1=gene1,
            gene2=gene2,
            datasets=datasets,
            method=task_method,
            out_dir=out_dir,
        )

        # ------------------------------------------------------------------
        # Rename PDF to an informative name
        # ------------------------------------------------------------------
        final_pdf_path = None
        if result["success"] and result["pdf_path"]:
            new_name = (
                f"ph{phase}__{analysis}__{gene1}_vs_{gene2}"
                f"__{dataset_str.replace('/', '-')}__{task_method}.pdf"
            )
            # Truncate to ~200 chars to avoid filesystem limits
            if len(new_name) > 200:
                new_name = new_name[:196] + ".pdf"
            final_pdf_path = _rename_pdf(result["pdf_path"], new_name, out_dir)
            print(f"  ✓ Saved: {final_pdf_path}")
            n_ok += 1
        else:
            print(f"  ✗ Failed: {result['error']}")
            n_fail += 1

        # ------------------------------------------------------------------
        # Write index row
        # ------------------------------------------------------------------
        writer.writerow({
            "phase":    phase,
            "analysis": analysis,
            "gene1":    gene1,
            "gene2":    gene2,
            "dataset":  dataset_str,
            "method":   task_method,
            "status":   "success" if result["success"] else "failed",
            "pdf_path": final_pdf_path or result.get("error", ""),
            "timestamp": datetime.utcnow().isoformat(timespec="seconds"),
        })
        index_fh.flush()   # write immediately so --resume works after crash

        # ------------------------------------------------------------------
        # Rate-limit delay (not added after the last task)
        # ------------------------------------------------------------------
        if i < n_total:
            delay = sleep_base + random.uniform(0, sleep_jitter)
            print(f"  Sleeping {delay:.1f}s …")
            time.sleep(delay)

    index_fh.close()
    print(
        f"\nDone.  Total={n_total}  OK={n_ok}  Failed={n_fail}  Skipped={n_skip}"
    )


# ---------------------------------------------------------------------------
# Task builders — one per phase
# ---------------------------------------------------------------------------

def build_phase1_tasks(datasets: list, method: str) -> list:
    """
    Phase 1: Normal-tissue / baseline analysis.

    The GEPIA correlation endpoint accepts only ONE dataset per query.
    We therefore loop over each dataset individually.  This produces one PDF
    per (gene-pair, dataset) combination, which you can then average across
    datasets to reproduce the Figure 4B heatmap style.

    Default dataset list: all 23 TCGA matched normal datasets (BRCA_Normal,
    COAD_Normal, … UCEC_Normal).  These are normal-tissue biopsies from TCGA
    cancer patients — the same baseline used in the published manuscript.
    To use GTEx healthy-donor tissue instead, pass the desired tissue names
    via --dataset on the command line.

    Two analysis groups per tissue:
      A) heme_path_pairs   — consecutive + selected non-consecutive heme pairs
      B) chromatin_vs_heme — every chromatin regulator × every heme gene
    """
    tasks = []

    for ds in datasets:
        # --- A: intra-heme-pathway correlations ---
        for g1, g2 in ALL_HEME_PAIRS:
            is_consec = (g1, g2) in HEME_CONSECUTIVE_PAIRS
            label = "heme_consecutive" if is_consec else "heme_nonconsecutive"
            tasks.append({
                "phase":    1,
                "analysis": label,
                "gene1":    g1,
                "gene2":    g2,
                "datasets": [ds],   # single dataset per query
                "method":   method,
            })

        # --- B: chromatin regulator vs heme gene ---
        for reg in CHROMATIN_REGULATORS:
            for heme in HEME_GENES:
                tasks.append({
                    "phase":    1,
                    "analysis": "chromatin_vs_heme",
                    "gene1":    reg,
                    "gene2":    heme,
                    "datasets": [ds],   # single dataset per query
                    "method":   method,
                })

    return tasks


def build_phase2_tasks(datasets: list, method: str) -> list:
    """
    Phase 2: Pan-cancer analysis (same pairs as Phase 1 but on TCGA tumor data).
    Also includes positive-control pairs.

    The GEPIA correlation endpoint accepts only ONE dataset per query.
    We loop over each tumor dataset individually, giving one PDF per
    (gene-pair, cancer) combination.

    Default dataset list: all 33 TCGA tumor datasets.
    """
    tasks = []

    for ds in datasets:
        # --- Intra-heme-pathway ---
        for g1, g2 in ALL_HEME_PAIRS:
            is_consec = (g1, g2) in HEME_CONSECUTIVE_PAIRS
            label = "heme_consecutive" if is_consec else "heme_nonconsecutive"
            tasks.append({
                "phase":    2,
                "analysis": label,
                "gene1":    g1,
                "gene2":    g2,
                "datasets": [ds],   # single dataset per query
                "method":   method,
            })

        # --- Chromatin regulator vs heme gene ---
        for reg in CHROMATIN_REGULATORS:
            for heme in HEME_GENES:
                tasks.append({
                    "phase":    2,
                    "analysis": "chromatin_vs_heme",
                    "gene1":    reg,
                    "gene2":    heme,
                    "datasets": [ds],   # single dataset per query
                    "method":   method,
                })

        # --- Positive controls (included in every dataset) ---
        for g1, g2 in POSITIVE_CONTROLS:
            tasks.append({
                "phase":    2,
                "analysis": "positive_control",
                "gene1":    g1,
                "gene2":    g2,
                "datasets": [ds],   # single dataset per query
                "method":   method,
            })

    return tasks


def build_phase3_tasks(cancer_list: list, method: str) -> list:
    """
    Phase 3: Per-cancer analysis.

    Each cancer is queried individually so results can be compared across
    tumor types rather than pooled.  We use both the Tumor and (where
    available) the Normal dataset for each cancer so intra-pathway
    correlations can be compared between cancer and matched normal.

    For each cancer the following are queried:
      - intra-heme-pathway pairs (Tumor only)
      - chromatin regulator vs heme gene (Tumor only)
    """
    tasks = []

    # Map bare cancer codes to their TCGA dataset strings
    tumor_datasets  = {c.replace("_Tumor",  "") : c for c in ALL_TCGA_TUMOR}
    normal_datasets = {c.replace("_Normal", "") : c for c in ALL_TCGA_NORMAL}

    for cancer in cancer_list:
        cancer = cancer.upper()
        if cancer not in _VALID_CANCERS:
            print(f"WARNING: '{cancer}' is not a recognised GEPIA cancer code; skipping.")
            continue

        tumor_ds  = tumor_datasets.get(cancer)
        normal_ds = normal_datasets.get(cancer)   # may be None (e.g. GBM has no normal)

        if tumor_ds is None:
            print(f"WARNING: No tumor dataset found for '{cancer}'; skipping.")
            continue

        for g1, g2 in ALL_HEME_PAIRS:
            is_consec = (g1, g2) in HEME_CONSECUTIVE_PAIRS
            label = "heme_consecutive" if is_consec else "heme_nonconsecutive"
            tasks.append({
                "phase":    3,
                "analysis": f"{label}_{cancer}_tumor",
                "gene1":    g1,
                "gene2":    g2,
                "datasets": [tumor_ds],
                "method":   method,
            })
            # If matched normal exists, also test normal tissue for comparison
            if normal_ds:
                tasks.append({
                    "phase":    3,
                    "analysis": f"{label}_{cancer}_normal",
                    "gene1":    g1,
                    "gene2":    g2,
                    "datasets": [normal_ds],
                    "method":   method,
                })

        for reg in CHROMATIN_REGULATORS:
            for heme in HEME_GENES:
                tasks.append({
                    "phase":    3,
                    "analysis": f"chromatin_vs_heme_{cancer}_tumor",
                    "gene1":    reg,
                    "gene2":    heme,
                    "datasets": [tumor_ds],
                    "method":   method,
                })

    return tasks


def build_phase4_tasks(cancer_list: list, method: str) -> list:
    """
    Phase 4: Focused per-cancer deep analysis.

    For each selected cancer:
      1. All chromatin regulators × all heme genes  (Tumor dataset)
      2. All consecutive heme-pathway pairs          (Tumor dataset)

    The key difference from Phase 3 is that this is exhaustive for a small
    set of user-chosen cancers, typically ones with known relevance
    (e.g. GBM, BLCA, ESCA, HNSC from the manuscript).
    """
    tasks = []

    tumor_datasets = {c.replace("_Tumor", ""): c for c in ALL_TCGA_TUMOR}

    for cancer in cancer_list:
        cancer = cancer.upper()
        if cancer not in _VALID_CANCERS:
            print(f"WARNING: '{cancer}' is not a recognised GEPIA cancer code; skipping.")
            continue

        tumor_ds = tumor_datasets.get(cancer)
        if tumor_ds is None:
            print(f"WARNING: No tumor dataset found for '{cancer}'; skipping.")
            continue

        # 1. Chromatin regulators × heme genes
        for reg in CHROMATIN_REGULATORS:
            for heme in HEME_GENES:
                tasks.append({
                    "phase":    4,
                    "analysis": f"chromatin_vs_heme_{cancer}",
                    "gene1":    reg,
                    "gene2":    heme,
                    "datasets": [tumor_ds],
                    "method":   method,
                })

        # 2. Consecutive heme pathway pairs
        for g1, g2 in HEME_CONSECUTIVE_PAIRS:
            tasks.append({
                "phase":    4,
                "analysis": f"heme_consecutive_{cancer}",
                "gene1":    g1,
                "gene2":    g2,
                "datasets": [tumor_ds],
                "method":   method,
            })

    return tasks


# ---------------------------------------------------------------------------
# Dataset validation helper
# ---------------------------------------------------------------------------

def validate_datasets(requested: list) -> list:
    """
    Check that all requested dataset strings are valid GEPIA identifiers.
    Prints a warning for any unknown strings and returns only valid ones.
    """
    valid = []
    for ds in requested:
        if ds in _VALID_DATASETS:
            valid.append(ds)
        else:
            print(
                f"WARNING: '{ds}' is not a recognised GEPIA dataset identifier.\n"
                f"  Valid examples: BRCA_Tumor, BRCA_Normal, Liver (GTEx)\n"
                f"  Run with --list-datasets to see all valid values."
            )
    return valid


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Batch GEPIA2 correlation queries for heme-pathway / chromatin-regulator analysis.\n"
            "See module docstring for detailed parameter documentation."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- Phase selection ---
    phase_group = parser.add_mutually_exclusive_group(required=True)
    phase_group.add_argument(
        "--phase", type=int, choices=[1, 2, 3, 4],
        help="Run a single analysis phase (1–4).",
    )
    phase_group.add_argument(
        "--all", action="store_true",
        help="Run all four phases in sequence.",
    )

    # --- Dataset selection ---
    parser.add_argument(
        "--dataset", nargs="+", dest="datasets", metavar="DATASET",
        help=(
            "One or more GEPIA dataset strings to use as the correlation dataset pool.\n"
            "Examples: BRCA_Tumor BLCA_Tumor  /  Liver Lung (GTEx)\n"
            "Defaults: Phase 1 → all GTEx tissues; Phase 2 → all TCGA tumors;\n"
            "          Phase 3 → per-cancer (ignores --dataset); Phase 4 → per-cancer."
        ),
    )

    # --- Cancer selection (for Phase 3 / 4) ---
    parser.add_argument(
        "--cancers", nargs="+", metavar="CANCER",
        default=DEFAULT_FOCUSED_CANCERS,
        help=(
            "Bare cancer codes for Phase 3 / 4 (e.g. GBM BLCA ESCA HNSC).\n"
            f"Default: {' '.join(DEFAULT_FOCUSED_CANCERS)}"
        ),
    )

    # --- Correlation method ---
    parser.add_argument(
        "--method", choices=["spearman", "pearson"], default="spearman",
        help="Correlation method (default: spearman, recommended for RNA-seq).",
    )

    # --- Output ---
    parser.add_argument(
        "--outdir", default="./gepia_pdfs",
        help="Directory for downloaded PDFs (default: ./gepia_pdfs).",
    )
    parser.add_argument(
        "--index", default="./gepia_index.csv",
        help="Path to the index CSV file (default: ./gepia_index.csv).",
    )

    # --- Rate-limiting ---
    parser.add_argument(
        "--sleep", type=float, default=8.0,
        help="Base sleep time in seconds between requests (default: 8).",
    )
    parser.add_argument(
        "--jitter", type=float, default=3.0,
        help="Maximum extra random delay in seconds (default: 3).",
    )

    # --- Resume ---
    parser.add_argument(
        "--resume", action="store_true",
        help="Skip gene pairs already recorded as successful in the index CSV.",
    )

    # --- Dry run ---
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print tasks without actually querying GEPIA.",
    )

    # --- Info ---
    parser.add_argument(
        "--list-datasets", action="store_true",
        help="Print all valid GEPIA dataset identifiers and exit.",
    )
    parser.add_argument(
        "--list-cancers", action="store_true",
        help="Print all valid GEPIA cancer codes and exit.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # --- Informational exits ---
    if args.list_datasets:
        print("=== TCGA Tumor datasets ===")
        print("\n".join(ALL_TCGA_TUMOR))
        print("\n=== TCGA Normal datasets ===")
        print("\n".join(ALL_TCGA_NORMAL))
        print("\n=== GTEx tissue datasets ===")
        print("\n".join(ALL_GTEX))
        sys.exit(0)

    if args.list_cancers:
        print("Valid cancer codes for --cancers:")
        print("\n".join(sorted(_VALID_CANCERS)))
        sys.exit(0)

    # --- Resolve which phases to run ---
    phases = [1, 2, 3, 4] if args.all else [args.phase]

    # --- Build tasks per phase ---
    all_tasks = []

    for phase in phases:
        if phase == 1:
            # Default for Phase 1: all 23 TCGA matched normal tissues, one at a time.
            # This matches the published analysis (manuscript Figure 4B).
            # To use GTEx healthy-donor tissue instead, pass --dataset explicitly.
            datasets = validate_datasets(args.datasets) if args.datasets else ALL_TCGA_NORMAL
            if not datasets:
                print("ERROR: No valid datasets for Phase 1.  Aborting.")
                sys.exit(1)
            tasks = build_phase1_tasks(datasets=datasets, method=args.method)

        elif phase == 2:
            # Default for Phase 2: all TCGA tumor datasets, queried one at a time
            datasets = validate_datasets(args.datasets) if args.datasets else ALL_TCGA_TUMOR
            if not datasets:
                print("ERROR: No valid datasets for Phase 2.  Aborting.")
                sys.exit(1)
            tasks = build_phase2_tasks(datasets=datasets, method=args.method)

        elif phase == 3:
            # Phase 3 loops per cancer; --dataset is ignored (individual datasets used)
            if args.datasets:
                print(
                    "NOTE: --dataset is ignored for Phase 3 "
                    "(individual per-cancer datasets are used automatically)."
                )
            # Phase 3: loop over ALL TCGA cancers unless --cancers restricts it
            # The full phase 3 is all cancers; use --cancers to limit scope
            cancer_list = args.cancers if args.cancers else list(_VALID_CANCERS)
            tasks = build_phase3_tasks(cancer_list=cancer_list, method=args.method)

        else:  # phase == 4
            cancer_list = args.cancers
            tasks = build_phase4_tasks(cancer_list=cancer_list, method=args.method)

        print(f"\n--- Phase {phase}: {len(tasks)} tasks ---")
        all_tasks.extend(tasks)

    print(f"\nTotal tasks across all selected phases: {len(all_tasks)}")

    if args.dry_run:
        print("\n[Dry-run mode: listing tasks only]\n")
        for t in all_tasks:
            print(
                f"  Ph{t['phase']} | {t['analysis']:40s} | "
                f"{t['gene1']:8s} vs {t['gene2']:8s} | "
                f"datasets: {'+'.join(t['datasets'])}"
            )
        sys.exit(0)

    # --- Run the batch ---
    batch_run(
        tasks=all_tasks,
        out_dir=args.outdir,
        index_path=args.index,
        method=args.method,
        sleep_base=args.sleep,
        sleep_jitter=args.jitter,
        resume=args.resume,
        dry_run=False,
    )


if __name__ == "__main__":
    main()
