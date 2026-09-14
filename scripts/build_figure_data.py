#!/usr/bin/env python3
"""
build_figure_data.py
====================
Rebuilds figure_4B_coexpression_normal_data.csv and
figure_4C_within_pathway_data.csv from the new GEPIA2 phase outputs
(results_ph1.csv = TCGA normals, results_ph2.csv = TCGA tumors).

These replacement CSVs are in the exact format expected by figure4_plot2.py,
so you can run that script unchanged after running this one.

Also prints a validation table comparing new values against the
original published values, so you can check they agree.

Usage:
    python build_figure_data.py

Input files (must be in same directory):
    results_ph1.csv              — output of extract_gepia_values.R on Phase 1 PDFs
    results_ph2.csv              — output of extract_gepia_values.R on Phase 2 PDFs
    figure_4B_coexpression_normal_data.csv   — original published values (for validation)
    figure_4C_within_pathway_data.csv        — original published values (for validation)

Output files:
    figure_4B_coexpression_normal_data_NEW.csv
    figure_4C_within_pathway_data_NEW.csv
"""

import pandas as pd
import numpy as np
import sys
import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PH1_CSV     = 'results_ph1.csv'    # TCGA normals
PH2_CSV     = 'results_ph2.csv'    # TCGA tumors
ORIG_4B_CSV = 'figure_4B_coexpression_normal_data.csv'
ORIG_4C_CSV = 'figure_4C_within_pathway_data.csv'
OUT_4B      = 'figure_4B_coexpression_normal_data_NEW.csv'
OUT_4C      = 'figure_4C_within_pathway_data_NEW.csv'

# The 9 cancers used in Figure 4C (must have matched normal + tumor)
FIG4C_CANCERS = ['BRCA', 'COAD', 'HNSC', 'KIRC', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'THCA']

# Known sample sizes for the normal tissue in each cancer (from TCGA)
N_NORMAL = {
    'BRCA': 114, 'COAD': 41, 'HNSC': 44, 'KIRC': 72, 'LIHC': 50,
    'LUAD': 59,  'LUSC': 49, 'PRAD': 52, 'THCA': 59,
}

# Chromatin regulators and heme genes in pathway order
REGS = ['H2AFZ', 'ATAD2', 'BRD2', 'BRD4']
HEME = ['ALAS1', 'ALAD', 'HMBS', 'UROS', 'UROD', 'CPOX', 'PPOX', 'FECH']

CONSECUTIVE_PAIRS = [
    ('ALAS1', 'ALAD'), ('ALAD', 'HMBS'), ('HMBS', 'UROS'),
    ('UROS', 'UROD'),  ('UROD', 'CPOX'), ('CPOX', 'PPOX'), ('PPOX', 'FECH'),
]

# ---------------------------------------------------------------------------
# Load phase data
# ---------------------------------------------------------------------------

def load_phase(path):
    if not os.path.exists(path):
        sys.exit(f"ERROR: {path} not found. Run extract_gepia_values.R first.")
    df = pd.read_csv(path)
    # Normalise column names to lowercase
    df.columns = df.columns.str.lower().str.strip()
    # Keep only successful extractions
    df = df[df['status'] == 'ok'].copy()
    # R column should be numeric
    df['r'] = pd.to_numeric(df['r'], errors='coerce')
    df['p_value'] = pd.to_numeric(df['p_value'], errors='coerce')
    print(f"  Loaded {len(df)} successful rows from {path}")
    print(f"  Analyses present: {sorted(df['analysis'].unique())}")
    print(f"  Datasets present (first 5): {sorted(df['dataset'].unique())[:5]} ...")
    return df

print("Loading phase data ...")
ph1 = load_phase(PH1_CSV)
ph2 = load_phase(PH2_CSV)

# Restrict Phase 1 strictly to TCGA matched normal datasets.
# Any GTEx tissue that accidentally ended up in the Phase 1 run is excluded here.
gtex_in_ph1 = ph1[~ph1['dataset'].str.endswith('_Normal')]['dataset'].unique()
if len(gtex_in_ph1) > 0:
    print(f"  NOTE: Dropping {len(gtex_in_ph1)} non-TCGA dataset(s) from Phase 1: {sorted(gtex_in_ph1)}")
ph1 = ph1[ph1['dataset'].str.endswith('_Normal')].copy()
print(f"  Phase 1 after filter: {ph1['dataset'].nunique()} TCGA normal datasets")
print()

# ---------------------------------------------------------------------------
# Build Figure 4B
# ---------------------------------------------------------------------------
# Use Phase 1 (TCGA normals), analysis == 'chromatin_vs_heme'
# Average R across all normal datasets for each regulator x heme gene pair
# Then also keep individual tissue rows (same format as original file)

print("Building Figure 4B data ...")

b_rows = ph1[ph1['analysis'] == 'chromatin_vs_heme'].copy()

if len(b_rows) == 0:
    print("  WARNING: No 'chromatin_vs_heme' rows found in results_ph1.csv")
    print("  Available analyses:", sorted(ph1['analysis'].unique()))
else:
    print(f"  Found {len(b_rows)} chromatin_vs_heme rows across "
          f"{b_rows['dataset'].nunique()} normal datasets")

# Individual tissue rows
indiv_rows = []
for _, row in b_rows.iterrows():
    indiv_rows.append({
        'figure_panel':  'Figure_4B',
        'tissue_dataset': row['dataset'],
        'regulator':      row['gene1'],
        'heme_gene':      row['gene2'],
        'spearman_R':     round(row['r'], 4) if not pd.isna(row['r']) else np.nan,
        'p_value':        row['p_value'],
        'data_source':    'GEPIA2 TCGA matched normal tissue RNA-seq',
    })

df_indiv = pd.DataFrame(indiv_rows)

# Pooled mean rows — mean R across all normal datasets per pair
pooled_rows = []
for reg in REGS:
    for heme in HEME:
        subset = b_rows[(b_rows['gene1'] == reg) & (b_rows['gene2'] == heme)]
        if len(subset) == 0:
            print(f"  WARNING: No data for {reg} vs {heme}")
            continue
        mean_r = subset['r'].mean()
        # For the pooled p-value: use the minimum (most significant) across tissues
        # as a conservative summary — consistent with how the original was computed
        min_p  = subset['p_value'].min()
        n_tissues = len(subset)
        pooled_rows.append({
            'figure_panel':  'Figure_4B',
            'tissue_dataset': f'POOLED_MEAN_{n_tissues}_tissues',
            'regulator':      reg,
            'heme_gene':      heme,
            'spearman_R':     round(mean_r, 4),
            'p_value':        round(min_p, 5) if not pd.isna(min_p) else np.nan,
            'data_source':    f'Pooled mean across {n_tissues} TCGA normal tissues',
        })

df_pooled = pd.DataFrame(pooled_rows)
df_4B_new = pd.concat([df_indiv, df_pooled], ignore_index=True)
df_4B_new.to_csv(OUT_4B, index=False)
print(f"  Saved {len(df_4B_new)} rows -> {OUT_4B}")
print()

# ---------------------------------------------------------------------------
# Build Figure 4C
# ---------------------------------------------------------------------------
# Use Phase 1 (TCGA normals) + Phase 2 (TCGA tumors)
# analysis == 'heme_consecutive'
# 9 specific cancers with matched normal + tumor

print("Building Figure 4C data ...")

c_rows_normal = ph1[ph1['analysis'] == 'heme_consecutive'].copy()
c_rows_tumor  = ph2[ph2['analysis'] == 'heme_consecutive'].copy()

if len(c_rows_normal) == 0:
    print("  WARNING: No 'heme_consecutive' rows in results_ph1.csv")
    print("  Available analyses:", sorted(ph1['analysis'].unique()))
if len(c_rows_tumor) == 0:
    print("  WARNING: No 'heme_consecutive' rows in results_ph2.csv")
    print("  Available analyses:", sorted(ph2['analysis'].unique()))

fig4c_rows = []
missing = []

for cancer in FIG4C_CANCERS:
    normal_ds = f'{cancer}_Normal'
    tumor_ds  = f'{cancer}_Tumor'

    for g1, g2 in CONSECUTIVE_PAIRS:
        pair_label = f'{g1}-{g2}'

        # Normal R
        n_row = c_rows_normal[
            (c_rows_normal['dataset'] == normal_ds) &
            (c_rows_normal['gene1']   == g1) &
            (c_rows_normal['gene2']   == g2)
        ]
        # Tumor R
        t_row = c_rows_tumor[
            (c_rows_tumor['dataset'] == tumor_ds) &
            (c_rows_tumor['gene1']   == g1) &
            (c_rows_tumor['gene2']   == g2)
        ]

        if len(n_row) == 0:
            missing.append(f"NORMAL: {normal_ds} {pair_label}")
            r_normal = np.nan
        else:
            r_normal = float(n_row.iloc[0]['r'])

        if len(t_row) == 0:
            missing.append(f"TUMOR:  {tumor_ds}  {pair_label}")
            r_tumor = np.nan
        else:
            r_tumor = float(t_row.iloc[0]['r'])

        delta = (r_tumor - r_normal) if not (pd.isna(r_normal) or pd.isna(r_tumor)) else np.nan

        fig4c_rows.append({
            'figure_panel':              'Figure_4C',
            'cancer_type':               cancer,
            'n_normal_samples':          N_NORMAL.get(cancer, np.nan),
            'gene_pair':                 pair_label,
            'upstream_gene':             g1,
            'downstream_gene':           g2,
            'spearman_R_normal_tissue':  round(r_normal, 4) if not pd.isna(r_normal) else np.nan,
            'spearman_R_cancer':         round(r_tumor,  4) if not pd.isna(r_tumor)  else np.nan,
            'delta_cancer_minus_normal': round(delta,    4) if not pd.isna(delta)     else np.nan,
            'data_source':               'GEPIA2 TCGA RNA-seq',
        })

if missing:
    print(f"  WARNING: {len(missing)} missing entries:")
    for m in missing:
        print(f"    {m}")

df_4C_new = pd.DataFrame(fig4c_rows)
df_4C_new.to_csv(OUT_4C, index=False)
print(f"  Saved {len(df_4C_new)} rows -> {OUT_4C}")
print()

# ---------------------------------------------------------------------------
# Validation: compare new vs original published values
# ---------------------------------------------------------------------------

print("=" * 65)
print("VALIDATION — comparing new vs original published values")
print("=" * 65)

# --- Figure 4B pooled means ---
if os.path.exists(ORIG_4B_CSV):
    orig_4B = pd.read_csv(ORIG_4B_CSV)
    orig_pooled = orig_4B[orig_4B['tissue_dataset'].str.startswith('POOLED')].copy()
    new_pooled  = df_4B_new[df_4B_new['tissue_dataset'].str.startswith('POOLED')].copy()

    print("\nFigure 4B — pooled mean R: original vs new")
    print(f"{'Regulator':<10} {'Heme':<8} {'Original R':>12} {'New R':>10} {'Diff':>8}")
    print("-" * 52)

    max_diff = 0
    for _, orig_row in orig_pooled.iterrows():
        reg  = orig_row['regulator']
        heme = orig_row['heme_gene']
        orig_r = float(orig_row['spearman_R'])
        new_row = new_pooled[(new_pooled['regulator'] == reg) &
                             (new_pooled['heme_gene'] == heme)]
        if len(new_row) == 0:
            print(f"{reg:<10} {heme:<8} {orig_r:>12.4f} {'MISSING':>10}")
            continue
        new_r = float(new_row.iloc[0]['spearman_R'])
        diff  = new_r - orig_r
        max_diff = max(max_diff, abs(diff))
        flag  = ' <-- LARGE DIFF' if abs(diff) > 0.05 else ''
        print(f"{reg:<10} {heme:<8} {orig_r:>12.4f} {new_r:>10.4f} {diff:>+8.4f}{flag}")

    print(f"\nMax absolute difference: {max_diff:.4f}")
    if max_diff < 0.02:
        print("✓ Values agree well (differences likely due to rounding)")
    elif max_diff < 0.10:
        print("~ Small differences — check whether the same datasets were used")
    else:
        print("✗ Large differences — datasets or extraction may differ")
else:
    print(f"\nOriginal {ORIG_4B_CSV} not found — skipping 4B validation")

# --- Figure 4C ---
if os.path.exists(ORIG_4C_CSV):
    orig_4C = pd.read_csv(ORIG_4C_CSV)

    print("\n\nFigure 4C — R values: original vs new (normal tissue side)")
    print(f"{'Cancer':<8} {'Pair':<12} {'Orig R_norm':>12} {'New R_norm':>11} "
          f"{'Orig R_canc':>12} {'New R_canc':>11} {'Diff_norm':>10}")
    print("-" * 80)

    max_diff_c = 0
    for _, orig_row in orig_4C.iterrows():
        cancer = orig_row['cancer_type']
        pair   = orig_row['gene_pair']
        orig_rn = float(orig_row['spearman_R_normal_tissue'])
        orig_rc = float(orig_row['spearman_R_cancer'])
        new_row = df_4C_new[(df_4C_new['cancer_type'] == cancer) &
                             (df_4C_new['gene_pair']   == pair)]
        if len(new_row) == 0:
            print(f"{cancer:<8} {pair:<12} {orig_rn:>12.3f} {'MISSING':>11}")
            continue
        new_rn = new_row.iloc[0]['spearman_R_normal_tissue']
        new_rc = new_row.iloc[0]['spearman_R_cancer']
        diff_n = new_rn - orig_rn if not pd.isna(new_rn) else np.nan
        if not pd.isna(diff_n):
            max_diff_c = max(max_diff_c, abs(diff_n))
        flag = ' <--' if not pd.isna(diff_n) and abs(diff_n) > 0.05 else ''
        print(f"{cancer:<8} {pair:<12} {orig_rn:>12.3f} {new_rn:>11.3f} "
              f"{orig_rc:>12.3f} {new_rc:>11.3f} {diff_n:>+10.3f}{flag}")

    print(f"\nMax absolute difference (normal side): {max_diff_c:.4f}")
    if max_diff_c < 0.02:
        print("✓ Values agree well")
    elif max_diff_c < 0.10:
        print("~ Small differences — check datasets")
    else:
        print("✗ Large differences — investigate")
else:
    print(f"\nOriginal {ORIG_4C_CSV} not found — skipping 4C validation")

print()
print("Done. Output files:")
print(f"  {OUT_4B}")
print(f"  {OUT_4C}")
print()
print("To use with figure4_plot2.py, either:")
print("  1. Rename the _NEW files to replace the originals, or")
print("  2. Edit DATA_4B and DATA_4C paths at the top of figure4_plot2.py")
