#!/usr/bin/env python3
"""
Batch GEPIA correlation PDF runner for HEM/chromatin/control pipeline.

Automatically extracts R and p-value from each PDF result and writes them
directly into the index CSV alongside the file path.

Usage:
    python run_gepia_pdf_batch.py --verify
    python run_gepia_pdf_batch.py --phase 1 --sleep 10 --resume
    python run_gepia_pdf_batch.py --phase 2 --sleep 12 --resume
    python run_gepia_pdf_batch.py --phase 3 --sleep 15 --resume
    python run_gepia_pdf_batch.py --phase 5 --cancers GBM_Tumor DLBC_Tumor LAML_Tumor --sleep 12 --resume
"""

from __future__ import annotations

import argparse
import concurrent.futures
import threading
import csv
import random
import re
import shutil
import string
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import gepia
from pdfminer.high_level import extract_text

# ---------------------------------------------------------------------------
# Gene lists
# ---------------------------------------------------------------------------

HEM_GENES = ["ALAS1", "ALAD", "HMBS", "UROS", "UROD", "CPOX", "PPOX", "FECH"]
CHROMATIN_REGULATORS = ["H2AFZ", "ATAD2", "BRD4", "BRD2"]
PROLIFERATION_CONTROLS = ["MKI67", "PCNA"]
GENERAL_CHROMATIN_CONTROLS = ["HDAC1"]
SPECIFICITY_CONTROLS = ["FASN", "ACACA", "TFRC", "FTH1", "ACTB", "GAPDH"]

CONSECUTIVE_HEM_PAIRS = [
    ("ALAS1", "ALAD"),
    ("ALAD",  "HMBS"),
    ("HMBS",  "UROS"),
    ("UROS",  "UROD"),
    ("UROD",  "CPOX"),
    ("CPOX",  "PPOX"),
    ("PPOX",  "FECH"),
]
NONCONSECUTIVE_HEM_PAIRS = [
    ("HMBS", "FECH"),
    ("ALAS1", "FECH"),
    ("HMBS", "PPOX"),
]
POSITIVE_CONTROLS = [
    ("E2F1",  "CCNE1"),
    ("CDK2",  "CCNE1"),
]

# ---------------------------------------------------------------------------
# Dataset labels
# ---------------------------------------------------------------------------

NORMAL_TISSUE_DATASET = gepia.TCGANormal  # 23 matched TCGA normals — only valid normal option for correlation API
TCGA_NORMAL_DATASET   = gepia.TCGANormal # 23 matched TCGA normals (alternative baseline)
PANCANCER_DATASET = gepia.TCGATumor  # all 33 TCGA tumor types — from gepia package

DEFAULT_CANCER_DATASETS = [
    "ACC_Tumor",  "BLCA_Tumor", "BRCA_Tumor", "CESC_Tumor", "CHOL_Tumor",
    "COAD_Tumor", "DLBC_Tumor", "ESCA_Tumor", "GBM_Tumor",  "HNSC_Tumor",
    "KICH_Tumor", "KIRC_Tumor", "KIRP_Tumor", "LAML_Tumor", "LGG_Tumor",
    "LIHC_Tumor", "LUAD_Tumor", "LUSC_Tumor", "MESO_Tumor", "OV_Tumor",
    "PAAD_Tumor", "PCPG_Tumor", "PRAD_Tumor", "READ_Tumor", "SARC_Tumor",
    "SKCM_Tumor", "STAD_Tumor", "TGCT_Tumor", "THCA_Tumor", "THYM_Tumor",
    "UCEC_Tumor", "UCS_Tumor",  "UVM_Tumor",
]

TCGA_INDIVIDUAL_NORMALS = [
    "BLCA_Normal",  "BRCA_Normal",  "CESC_Normal",  "CHOL_Normal",
    "COAD_Normal",  "ESCA_Normal",  "HNSC_Normal",  "KICH_Normal",
    "KIRC_Normal",  "KIRP_Normal",  "LIHC_Normal",  "LUAD_Normal",
    "LUSC_Normal",  "PAAD_Normal",  "PCPG_Normal",  "PRAD_Normal",
    "READ_Normal",  "SARC_Normal",  "SKCM_Normal",  "STAD_Normal",
    "THCA_Normal",  "THYM_Normal",  "UCEC_Normal",
]

MAX_QUERIES_PER_SESSION = 50

# ---------------------------------------------------------------------------
# CSV columns
# ---------------------------------------------------------------------------

CSV_FIELDS = [
    "phase", "analysis", "gene1", "gene2", "dataset",
    "method", "status", "R", "p_value", "path_or_message",
]

# ---------------------------------------------------------------------------
# PDF extraction
# ---------------------------------------------------------------------------

_UNICODE_MINUS = "\u2212"

_R_PAT = re.compile(
    r'R\s*=\s*([' + _UNICODE_MINUS + r'-]?\d+(?:\.\d+)?(?:[eE]['
    + _UNICODE_MINUS + r'+-]\d+)?)'
)
_P_PAT = re.compile(
    r'p[' + _UNICODE_MINUS + r'-]value\s*=\s*([^\s]+)'
)


def _norm(s: str) -> str:
    return s.replace(_UNICODE_MINUS, "-")


def extract_r_and_p(pdf_path: str) -> Tuple[Optional[str], Optional[str]]:
    """Extract R and p-value from a GEPIA PDF. Returns normalised ASCII strings."""
    try:
        text = extract_text(pdf_path)
    except Exception:
        return None, None
    r_m = _R_PAT.search(text)
    p_m = _P_PAT.search(text)
    return (_norm(r_m.group(1)) if r_m else None,
            _norm(p_m.group(1)) if p_m else None)


# ---------------------------------------------------------------------------
# Task builder
# ---------------------------------------------------------------------------

def build_tasks(phases: Iterable[str], cancers=None) -> List[Dict]:
    tasks: List[Dict] = []
    phases = set(str(p) for p in phases)

    if "1" in phases:
        for g1, g2 in CONSECUTIVE_HEM_PAIRS + NONCONSECUTIVE_HEM_PAIRS:
            tasks.append(dict(phase="1A", analysis="normal_within_heme",
                              gene1=g1, gene2=g2, use_normal=True))
        for reg in CHROMATIN_REGULATORS:
            for hem in HEM_GENES:
                tasks.append(dict(phase="1B", analysis="normal_chromatin_vs_heme",
                                  gene1=reg, gene2=hem, use_normal=True))
        for ctrl in PROLIFERATION_CONTROLS + GENERAL_CHROMATIN_CONTROLS:
            for hem in HEM_GENES:
                tasks.append(dict(phase="1C", analysis="normal_confound_vs_heme",
                                  gene1=ctrl, gene2=hem, use_normal=True))
        for ctrl in SPECIFICITY_CONTROLS:
            tasks.append(dict(phase="1C", analysis="normal_h2afz_specificity",
                              gene1="H2AFZ", gene2=ctrl, use_normal=True))

    if "2" in phases:
        for g1, g2 in CONSECUTIVE_HEM_PAIRS + NONCONSECUTIVE_HEM_PAIRS:
            tasks.append(dict(phase="2A", analysis="pancancer_within_heme",
                              gene1=g1, gene2=g2))
        for reg in CHROMATIN_REGULATORS:
            for hem in HEM_GENES:
                tasks.append(dict(phase="2B", analysis="pancancer_chromatin_vs_heme",
                                  gene1=reg, gene2=hem))
        for ctrl in PROLIFERATION_CONTROLS + GENERAL_CHROMATIN_CONTROLS:
            for hem in HEM_GENES:
                tasks.append(dict(phase="2C", analysis="pancancer_confound_vs_heme",
                                  gene1=ctrl, gene2=hem))
        for ctrl in SPECIFICITY_CONTROLS:
            tasks.append(dict(phase="2C", analysis="pancancer_h2afz_specificity",
                              gene1="H2AFZ", gene2=ctrl))
        for g1, g2 in POSITIVE_CONTROLS:
            tasks.append(dict(phase="2D", analysis="positive_control",
                              gene1=g1, gene2=g2))

    if "3" in phases:
        selected = cancers or DEFAULT_CANCER_DATASETS
        for ds in selected:
            tasks.append(dict(phase="3A", analysis="per_cancer_h2afz_hmbs",
                              gene1="H2AFZ", gene2="HMBS", dataset_override=ds))
            tasks.append(dict(phase="3B", analysis="per_cancer_hmbs_uros",
                              gene1="HMBS", gene2="UROS", dataset_override=ds))

    if "5" in phases:
        selected = cancers or ["GBM_Tumor", "DLBC_Tumor", "LAML_Tumor", "BLCA_Tumor"]
        for ds in selected:
            for reg in CHROMATIN_REGULATORS:
                for hem in HEM_GENES:
                    tasks.append(dict(phase="5A", analysis="focused_chromatin_vs_heme",
                                      gene1=reg, gene2=hem, dataset_override=ds))
            for g1, g2 in CONSECUTIVE_HEM_PAIRS:
                tasks.append(dict(phase="5B", analysis="focused_within_heme",
                                  gene1=g1, gene2=g2, dataset_override=ds))
            for ctrl in ["MKI67", "HDAC1"]:
                tasks.append(dict(phase="5C", analysis="focused_confound_hmbs",
                                  gene1=ctrl, gene2="HMBS", dataset_override=ds))

    if "6" in phases:
        # Phase 6 — per individual normal tissue type (23 TCGA matched normals)
        # Mirrors Phase 1 structure but runs each normal tissue separately
        # allowing direct comparison with per-cancer Phase 5 data
        normals = cancers or TCGA_INDIVIDUAL_NORMALS
        for ds in normals:
            # 6A: within-pathway HEM correlations per normal tissue
            for g1, g2 in CONSECUTIVE_HEM_PAIRS + NONCONSECUTIVE_HEM_PAIRS:
                tasks.append(dict(phase="6A", analysis="normal_indiv_within_heme",
                                  gene1=g1, gene2=g2, dataset_override=ds))
            # 6B: chromatin regulators vs HEM genes per normal tissue
            for reg in CHROMATIN_REGULATORS:
                for hem in HEM_GENES:
                    tasks.append(dict(phase="6B", analysis="normal_indiv_chromatin_vs_heme",
                                      gene1=reg, gene2=hem, dataset_override=ds))
            # 6C: proliferation controls per normal tissue
            for ctrl in PROLIFERATION_CONTROLS + GENERAL_CHROMATIN_CONTROLS:
                for hem in HEM_GENES:
                    tasks.append(dict(phase="6C", analysis="normal_indiv_confound_vs_heme",
                                      gene1=ctrl, gene2=hem, dataset_override=ds))

    return tasks


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_name(s: str) -> str:
    keep = string.ascii_letters + string.digits + "_-"
    return "".join(ch if ch in keep else "_" for ch in s)



def dataset_label(ds: list) -> str:
    """
    Return a short human-readable label for use in filenames.
    Avoids the 255-char filename limit when ds contains all 33 cancer names.
    """
    import gepia as _gepia
    if ds == _gepia.TCGATumor:
        return "pancancer"
    if ds == _gepia.GTEx:
        return "GTEx_normal"
    if ds == _gepia.TCGANormal:
        return "TCGA_normal"
    if len(ds) == 1:
        return ds[0]
    # Small custom subset — join up to 3 names then summarise
    if len(ds) <= 3:
        return "_".join(ds)
    return "_".join(ds[:2]) + f"_and_{len(ds)-2}_more"

def resolve_dataset(task: Dict) -> List[str]:
    if "dataset_override" in task:
        return [task["dataset_override"]]
    if task.get("use_normal"):
        return NORMAL_TISSUE_DATASET
    return PANCANCER_DATASET


# ---------------------------------------------------------------------------
# GEPIA query
# ---------------------------------------------------------------------------

def run_one(gene1, gene2, dataset, method, out_dir, prefix):
    c = gepia.correlation()
    c.setOutDir(str(out_dir))
    c.setParams({
        "dataset":         dataset,
        "methodoption":    method,
        "signature1":      [gene1],
        "signature1_norm": "",
        "signature2":      [gene2],
        "signature2_norm": "",
    })
    result = c.query()
    if result is None:
        return "failed", "GEPIA returned None — check dataset labels."

    src = Path(str(result))
    if not src.is_absolute():
        for candidate in (Path.cwd() / src, out_dir / src.name):
            if candidate.exists():
                src = candidate
                break

    if not src.exists():
        return "unknown", f"GEPIA returned '{result}' but file not found."

    dest = out_dir / f"{safe_name(prefix)}__{gene1}_vs_{gene2}.pdf"
    if src.resolve() != dest.resolve():
        shutil.move(str(src), str(dest))
    return "ok", str(dest)


def run_one_with_retry(gene1, gene2, dataset, method, out_dir, prefix,
                       max_retries=3, base_sleep=8.0):
    for attempt in range(max_retries):
        try:
            status, msg = run_one(gene1, gene2, dataset, method, out_dir, prefix)
            if status == "ok":
                return status, msg
            wait = base_sleep * (2 ** attempt) + random.uniform(5, 15)
            print(f"    Attempt {attempt+1}/{max_retries} failed ({msg}), "
                  f"retrying in {wait:.0f}s ...")
            time.sleep(wait)
        except Exception as exc:
            wait = base_sleep * (2 ** attempt) + random.uniform(5, 15)
            print(f"    Exception on attempt {attempt+1}/{max_retries}: {exc}, "
                  f"retrying in {wait:.0f}s ...")
            time.sleep(wait)
    return "failed_all_retries", f"Failed after {max_retries} attempts"


# ---------------------------------------------------------------------------
# Sleep strategy
# ---------------------------------------------------------------------------

def compute_sleep(base_sleep, status, queries_done, batch_pause_every=20):
    if queries_done > 0 and queries_done % batch_pause_every == 0:
        pause = random.uniform(60, 120)
        print(f"    -- Batch pause after {queries_done} queries: {pause:.0f}s --")
        return pause
    if status in ("failed", "error", "failed_all_retries", "unknown"):
        wait = base_sleep * 3 + random.uniform(5, 15)
        print(f"    Error backoff: {wait:.0f}s")
        return wait
    return max(base_sleep + random.uniform(-2, 5), 5.0)


# ---------------------------------------------------------------------------
# Verify
# ---------------------------------------------------------------------------

def verify_connection(out_dir, method):
    print("=" * 60)
    print("Verifying GEPIA: H2AFZ vs HMBS on BRCA_Tumor ...")
    test_dir = out_dir / "_verify"
    test_dir.mkdir(parents=True, exist_ok=True)
    try:
        status, msg = run_one("H2AFZ", "HMBS", ["BRCA_Tumor"],
                              method, test_dir, "verify_test")
        if status == "ok":
            r_val, p_val = extract_r_and_p(msg)
            print(f"Connection OK  -->  R={r_val}  p={p_val}")
            print(f"PDF saved to:  {msg}")
            print("=" * 60)
            return True
        print(f"Query failed: {msg}")
    except Exception as exc:
        print(f"Exception: {exc}")
    print("Check dataset label format and GEPIA package version.")
    print("=" * 60)
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def extract_all_pdfs(out_dir: Path, method: str) -> int:
    """
    Post-processing: scan all PDFs in out_dir, extract R/p, write fresh CSV.
    Fast — no network queries. Run once after all sessions complete.
    """
    pdfs = sorted(out_dir.glob("*.pdf"))
    print(f"Extracting R/p from {len(pdfs)} PDFs in {out_dir.resolve()} ...")

    index_path = out_dir / "gepia_pdf_index.csv"
    ok = failed = 0

    with open(index_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()

        for pdf in pdfs:
            stem  = pdf.stem
            parts = stem.split("__")
            if len(parts) < 4 or "_vs_" not in parts[3]:
                continue
            phase    = parts[0]
            analysis = parts[1]
            dataset  = parts[2]
            gene1, gene2 = parts[3].split("_vs_", 1)

            r_val, p_val = extract_r_and_p(str(pdf))
            if r_val is not None:
                ok += 1
            else:
                failed += 1
                print(f"  FAIL: {pdf.name}")

            writer.writerow({
                'phase': phase, 'analysis': analysis,
                'gene1': gene1, 'gene2': gene2,
                'dataset': dataset, 'method': method,
                'status': 'ok' if r_val else 'extraction_failed',
                'R': r_val or '', 'p_value': p_val or '',
                'path_or_message': str(pdf),
            })

    print(f"Done. Extracted: {ok}  Failed: {failed}")
    print(f"CSV: {index_path}")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="Batch GEPIA runner — extracts R/p into index CSV on the fly."
    )
    ap.add_argument("--phase", action="append", choices=["1", "2", "3", "5", "6"])
    ap.add_argument("--all", action="store_true", help="Run all phases.")
    ap.add_argument("--verify", action="store_true", help="Test query and exit.")
    ap.add_argument("--extract-only", action="store_true",
                    help="Extract R/p from all PDFs in --out dir and write CSV. No queries.")
    ap.add_argument("--cancers", nargs="*", default=None)
    ap.add_argument("--method", default="spearman", choices=["pearson", "spearman"])
    ap.add_argument("--out", default="gepia_pdf_results")
    ap.add_argument("--sleep", type=float, default=10.0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=1,
                    help="Parallel query workers (default: 1). 3 workers recommended.")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--batch-pause-every", type=int, default=20)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.verify:
        return 0 if verify_connection(out_dir, args.method) else 1

    if getattr(args, 'extract_only', False):
        return extract_all_pdfs(out_dir, args.method)

    phases = ["1", "2", "3", "5", "6"] if args.all else (args.phase or ["2"])
    tasks  = build_tasks(phases, cancers=args.cancers)
    session_limit = args.limit if args.limit is not None else MAX_QUERIES_PER_SESSION

    print(f"Total tasks in plan : {len(tasks)}")
    print(f"Session query limit : {session_limit}  (skipped tasks don't count)")
    print(f"Base sleep          : {args.sleep}s + jitter")
    print(f"Batch pause every   : {args.batch_pause_every} queries")
    print(f"Correlation method  : {args.method}")
    print(f"Output directory    : {out_dir.resolve()}")
    print(f"Resume mode         : {'ON' if args.resume else 'OFF'}")
    if not args.resume:
        print("TIP: --resume skips already-completed queries safely.")
    print()

    index_path   = out_dir / "gepia_pdf_index.csv"
    write_header = not index_path.exists()
    queries_done = 0

    with open(index_path, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        if write_header:
            writer.writeheader()

        # Split tasks into skip (already done) and new
        skip_tasks, new_tasks = [], []
        for task in tasks:
            g1 = task['gene1']
            g2 = task['gene2']
            ds = resolve_dataset(task)
            pfx = f"{task['phase']}__{task['analysis']}__{dataset_label(ds)}"
            exp = out_dir / f"{safe_name(pfx)}__{g1}_vs_{g2}.pdf"
            if args.resume and exp.exists():
                skip_tasks.append((task, exp))
            else:
                new_tasks.append((task, exp))

        # Handle skips — just record path, skip PDF parsing during run.
        # Use --extract-only afterwards to populate R/p values.
        if skip_tasks:
            print(f"Skipping {len(skip_tasks)} already-done tasks instantly (no PDF parsing).")
        for task, expected in skip_tasks:
            gene1 = task['gene1']
            gene2 = task['gene2']
            ds    = resolve_dataset(task)
            writer.writerow({
                'phase': task['phase'], 'analysis': task['analysis'],
                'gene1': gene1, 'gene2': gene2,
                'dataset': ';'.join(ds), 'method': args.method,
                'status': 'skipped',
                'R': '', 'p_value': '',
                'path_or_message': str(expected),
            })
        fh.flush()

        # Apply session limit to new tasks
        if len(new_tasks) > session_limit:
            print(f"\nCapping to session limit: {session_limit} of {len(new_tasks)} new tasks.")
            print("Run again with --resume to continue.")
            new_tasks = new_tasks[:session_limit]

        workers = getattr(args, 'workers', 1)
        print(f"\nSkipped: {len(skip_tasks)}  |  New: {len(new_tasks)}  |  Workers: {workers}\n")

        # Thread-safe write
        import threading
        import concurrent.futures
        csv_lock = threading.Lock()

        def process_task(idx_task):
            idx, (task, expected) = idx_task
            gene1  = task['gene1']
            gene2  = task['gene2']
            ds     = resolve_dataset(task)
            prefix = f"{task['phase']}__{task['analysis']}__{dataset_label(ds)}"
            total  = len(new_tasks)
            print(f"[{idx}/{total}] START  {task['phase']}  "
                  f"{gene1} vs {gene2}  on {dataset_label(ds)}")
            try:
                status, msg = run_one_with_retry(
                    gene1, gene2, ds, args.method,
                    out_dir, prefix, max_retries=3, base_sleep=args.sleep)
            except Exception as exc:
                status, msg = 'error', repr(exc)
            r_val, p_val = None, None
            if status == 'ok':
                r_val, p_val = extract_r_and_p(msg)
                print(f"[{idx}/{total}] DONE   {gene1} vs {gene2}  "
                      f"R={r_val}  p={p_val}  {Path(msg).name}")
            else:
                print(f"[{idx}/{total}] FAIL   {gene1} vs {gene2}  {status}: {msg}")
            with csv_lock:
                writer.writerow({
                    'phase': task['phase'], 'analysis': task['analysis'],
                    'gene1': gene1, 'gene2': gene2,
                    'dataset': ';'.join(ds), 'method': args.method,
                    'status': status,
                    'R': r_val or '', 'p_value': p_val or '',
                    'path_or_message': msg,
                })
                fh.flush()
            return status

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(process_task, enumerate(new_tasks, 1)))

        queries_done = len(new_tasks)
        failed_count = sum(1 for r in results if r != 'ok')
        print(f"\nAll workers finished. Failed: {failed_count}/{len(new_tasks)}")

    print(f"\nSession complete.")
    print(f"New queries submitted : {queries_done}")
    print(f"Index CSV             : {index_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
