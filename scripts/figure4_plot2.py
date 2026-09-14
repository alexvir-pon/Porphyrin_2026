"""
figure4_plot.py
===============
Generates Figure 4 (panels A, B, C) for:

  "Chromatin regulator perturbation causes heme pathway imbalance
   in yeast and human cells"
  Mitkevich et al. 2026

Panels:
  A — Mammalian KD heatmap (H2AFZ and ATAD2 knockdown datasets)
  B — Chromatin regulator vs heme gene co-expression in normal tissue
  C — Within-pathway consecutive enzyme co-expression: normal vs cancer

Usage:
  python figure4_plot.py

Input files (must be in same directory or edit paths below):
  figure_4A_mammalian_KD_data.csv
  figure_4B_coexpression_normal_data.csv
  figure_4C_within_pathway_data.csv

Output:
  figure4_combined.png  (300 dpi, publication quality)
  figure4_combined.svg  (vector, for editing in Inkscape/Illustrator)

Requirements:
  pip install pandas numpy matplotlib openpyxl
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

# ── Configuration ─────────────────────────────────────────────────────────────
DATA_4A    = 'figure_4A_mammalian_KD_data.csv'
DATA_4B    = 'figure_4B_coexpression_normal_data_NEW.csv'
DATA_4C    = 'figure_4C_within_pathway_data_NEW.csv'
OUTPUT_PNG = 'figure4_combined.png'
OUTPUT_SVG = 'figure4_combined.svg'
DPI        = 300

# ── Font sizes ────────────────────────────────────────────────────────────────
FS_GENE      = 18   # gene names on x-axis (ALAS1, ALAD …)
FS_DATASET   = 16   # dataset labels on y-axis (WI38, U2OS …)
FS_REG       = 18   # regulator names on Panel B y-axis
FS_CELL      = 11   # numbers + stars inside heatmap cells (smaller — cells are tight)
FS_ND        = 10   # n.d. text
FS_CB_LABEL  = 13   # colourbar axis label
FS_CB_TICK   = 11   # colourbar tick numbers
FS_TITLE     = 18   # panel titles A, B, C
FS_YLABEL    = 16   # Panel C y-axis label
FS_XTICK_C   = 14   # Panel C junction labels (gene→gene)
FS_LEGEND    = 15   # Panel C legend

# ── Pathway order ─────────────────────────────────────────────────────────────
CORE = ['ALAS1', 'ALAD', 'HMBS', 'UROS', 'UROD', 'CPOX', 'PPOX', 'FECH']
REGS = ['H2AFZ', 'ATAD2', 'BRD2', 'BRD4']

# ── Colour map — blue → white → red ──────────────────────────────────────────
CMAP = LinearSegmentedColormap.from_list('bwr_custom',
    ['#2166ac', '#74add1', '#e0f3f8', '#f7f7f7', '#fee090', '#f46d43', '#a50026'])

# ── Regulator colours (matching snapshot) ────────────────────────────────────
REG_COLORS = {'H2AFZ': '#E07B00',   # orange  — Δhtz1
              'ATAD2': '#2166ac',   # blue    — Δyta7
              'BRD2':  '#008080',   # teal    — Δbdf1
              'BRD4':  '#008080'}

# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

df4A = pd.read_csv(DATA_4A)
df4B = pd.read_csv(DATA_4B)
df4C = pd.read_csv(DATA_4C)

# ── Panel A matrix ────────────────────────────────────────────────────────────
datasets_order = [
    ('H2AFZ', 'WI38',    'H2AFZ-siRNA\nWI38 (human)',       '#E07B00'),
    ('H2AFZ', 'U2OS',    'H2AFZ-siRNA\nU2OS (human)',       '#E07B00'),
    ('ATAD2', 'ES_cells','Atad2-siRNA\nES cells (mouse)',   '#2166ac'),
    ('ATAD2', 'C2C12',   'Atad2-shRNA\nC2C12 (mouse)',      '#2166ac'),
]

n_ds   = len(datasets_order)
n_gene = len(CORE)
lfc_mat  = np.full((n_ds, n_gene), np.nan)
padj_mat = np.full((n_ds, n_gene), np.nan)

for j, (reg, cell, label, col) in enumerate(datasets_order):
    sub = df4A[(df4A['regulator'] == reg) & (df4A['cell_line'] == cell)]
    for i, gene in enumerate(CORE):
        row = sub[sub['heme_gene_human'] == gene]
        if len(row):
            lfc_mat[j, i]  = float(row.iloc[0]['log2FC_KD_vs_ctrl'])
            padj_mat[j, i] = float(row.iloc[0]['padj']) if str(row.iloc[0]['padj']) != 'NA' else np.nan

# ── Panel B matrix ────────────────────────────────────────────────────────────
pooled = df4B[df4B['tissue_dataset'].str.startswith('POOLED_MEAN')].copy()

corr_mat = np.full((len(REGS), n_gene), np.nan)
pval_mat = np.full((len(REGS), n_gene), np.nan)

for i, reg in enumerate(REGS):
    for j, gene in enumerate(CORE):
        row = pooled[(pooled['regulator'] == reg) & (pooled['heme_gene'] == gene)]
        if len(row):
            corr_mat[i, j] = float(row.iloc[0]['spearman_R'])
            pval_mat[i, j] = float(row.iloc[0]['p_value'])

# ── Panel C data ──────────────────────────────────────────────────────────────
junctions = ['ALAS1-ALAD', 'ALAD-HMBS', 'HMBS-UROS', 'UROS-UROD',
             'UROD-CPOX', 'CPOX-PPOX', 'PPOX-FECH']

norm_vals   = [df4C[df4C['gene_pair'] == j]['spearman_R_normal_tissue'].mean() for j in junctions]
cancer_vals = [df4C[df4C['gene_pair'] == j]['spearman_R_cancer'].mean()        for j in junctions]

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE — make it taller to accommodate larger fonts
# ══════════════════════════════════════════════════════════════════════════════

fig = plt.figure(figsize=(20, 16))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(3, 1, figure=fig,
                       hspace=0.60,
                       left=0.14, right=0.96,
                       top=0.94, bottom=0.08,
                       height_ratios=[3, 3, 3.2])

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A
# ══════════════════════════════════════════════════════════════════════════════

ax_a   = fig.add_subplot(gs[0])
norm_a = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
im_a   = ax_a.imshow(lfc_mat, cmap=CMAP, norm=norm_a, aspect='auto')

ax_a.set_xticks(range(n_gene))
ax_a.set_xticklabels(CORE, fontsize=FS_GENE, fontweight='bold', color='#8B0000')
ax_a.xaxis.set_ticks_position('top')
ax_a.xaxis.set_label_position('top')

ax_a.set_yticks(range(n_ds))
ax_a.set_yticklabels([d[2] for d in datasets_order], fontsize=FS_DATASET)
for tk, (_, _, _, col) in zip(ax_a.get_yticklabels(), datasets_order):
    tk.set_color(col)
    tk.set_fontweight('bold')
ax_a.tick_params(length=0)

for j in range(n_ds):
    reg = datasets_order[j][0]
    for i in range(n_gene):
        lfc = lfc_mat[j, i]
        pv  = padj_mat[j, i]
        if np.isnan(lfc):
            ax_a.text(i, j, 'n.d.', ha='center', va='center',
                      fontsize=FS_ND, color='#bbb')
            continue
        if reg == 'H2AFZ':
            sig = '*'
        else:
            sig = ('***' if (not np.isnan(pv) and pv < 0.001) else
                   '**'  if (not np.isnan(pv) and pv < 0.01)  else
                   '*'   if (not np.isnan(pv) and pv < 0.05)  else '')
        txt = f'{lfc:+.2f}' + (f'\n{sig}' if sig else '')
        fc  = 'white' if abs(lfc) > 1.2 else '#1a1a1a'
        ax_a.text(i, j, txt, ha='center', va='center',
                  fontsize=FS_CELL, color=fc, fontweight='bold')

ax_a.axhline(1.5, color='#aaa', lw=1.5, ls='--', alpha=0.7)

cb_a = fig.colorbar(im_a, ax=ax_a, fraction=0.025, pad=0.01, shrink=0.9)
cb_a.set_label('log2FC\n(KD / ctrl)', fontsize=FS_CB_LABEL)
cb_a.ax.tick_params(labelsize=FS_CB_TICK)

ax_a.set_title('A   Heme pathway response to chromatin regulator knockdown',
               fontsize=FS_TITLE, fontweight='bold', color='#1a1a2e', pad=32, loc='left')
for sp in ax_a.spines.values():
    sp.set_linewidth(1.3)
    sp.set_color('#555')

# ══════════════════════════════════════════════════════════════════════════════
# PANEL B
# ══════════════════════════════════════════════════════════════════════════════

ax_b   = fig.add_subplot(gs[1])
norm_b = TwoSlopeNorm(vmin=-0.75, vcenter=0, vmax=0.75)
im_b   = ax_b.imshow(corr_mat, cmap=CMAP, norm=norm_b, aspect='auto')

ax_b.set_xticks(range(n_gene))
ax_b.set_xticklabels(CORE, fontsize=FS_GENE, fontweight='bold', color='#8B0000')
ax_b.xaxis.set_ticks_position('top')
ax_b.xaxis.set_label_position('top')

ax_b.set_yticks(range(len(REGS)))
ax_b.set_yticklabels(REGS, fontsize=FS_REG, fontweight='bold')
for tk, reg in zip(ax_b.get_yticklabels(), REGS):
    tk.set_color(REG_COLORS.get(reg, '#1a1a1a'))
ax_b.tick_params(length=0)

for i in range(len(REGS)):
    for j in range(n_gene):
        v = corr_mat[i, j]
        p = pval_mat[i, j]
        if np.isnan(v):
            continue
        sig = ('***' if (not np.isnan(p) and p < 0.001) else
               '**'  if (not np.isnan(p) and p < 0.01)  else
               '*'   if (not np.isnan(p) and p < 0.05)  else '')
        txt = f'{v:+.2f}' + (f'\n{sig}' if sig else '')
        fc  = 'white' if abs(v) > 0.45 else '#1a1a1a'
        ax_b.text(j, i, txt, ha='center', va='center',
                  fontsize=FS_CELL, color=fc, fontweight='bold')

cb_b = fig.colorbar(im_b, ax=ax_b, fraction=0.025, pad=0.01, shrink=0.9)
cb_b.set_label('Spearman R\n(TCGA matched normals)', fontsize=FS_CB_LABEL)
cb_b.ax.tick_params(labelsize=FS_CB_TICK)

ax_b.set_title('B   Chromatin regulator – heme pathway co-expression (TCGA matched normal tissues, mean of 23 types)',
               fontsize=FS_TITLE, fontweight='bold', color='#1a1a2e', pad=32, loc='left')
for sp in ax_b.spines.values():
    sp.set_linewidth(1.3)
    sp.set_color('#555')

# ══════════════════════════════════════════════════════════════════════════════
# PANEL C
# ══════════════════════════════════════════════════════════════════════════════

ax_c = fig.add_subplot(gs[2])
x = np.arange(len(junctions))
w = 0.32

ax_c.bar(x - w/2, norm_vals,   width=w, color='#4393c3', alpha=0.9,
         label='Normal tissue', edgecolor='#2166ac', linewidth=0.8)
ax_c.bar(x + w/2, cancer_vals, width=w, color='#f4a582', alpha=0.9,
         label='Cancer',        edgecolor='#b2182b', linewidth=0.8)

ax_c.axhline(0, color='#555', lw=0.8)
ax_c.set_xticks(x)
ax_c.set_xticklabels([f'{j.split("-")[0]}\n→\n{j.split("-")[1]}' for j in junctions],
                     fontsize=FS_XTICK_C, ha='center', color='#555')
ax_c.set_ylabel('Mean Spearman R\n(9 matched cancer types)', fontsize=FS_YLABEL)
ax_c.tick_params(axis='y', labelsize=FS_XTICK_C)
ax_c.set_ylim(-0.42, 0.92)
ax_c.legend(fontsize=FS_LEGEND, loc='upper right', framealpha=0.85)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)
ax_c.set_facecolor('white')

ax_c.set_title('C   Within-pathway co-expression: normal tissue vs cancer',
               fontsize=FS_TITLE, fontweight='bold', color='#1a1a2e', pad=8, loc='left')

# ── Save ──────────────────────────────────────────────────────────────────────
plt.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_SVG, bbox_inches='tight', facecolor='white')
plt.close()

print(f"Saved: {OUTPUT_PNG}")
print(f"Saved: {OUTPUT_SVG}")
