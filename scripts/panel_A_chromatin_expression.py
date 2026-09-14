"""
panel_A_chromatin_expression.py
--------------------------------
Generates Panel A: chromatin regulator expression changes across 33 TCGA
cancer types, shown as coloured arrows (red = significant up, blue =
significant down, grey = trend/ns, dash = no change).

Input:
    table_A_chromatin_expression.csv   (hand-corrected expression table)

Output:
    panel_A_chromatin_expression.png

Columns expected in CSV:
    cancer, in_main_figure,
    H2AFZ, H2AFZ_sign, ATAD2, ATAD2_sign, BRD4, BRD4_sign, BRD2, BRD2_sign

Values in expression columns:  'up' | 'down' | 'ns' | 'no' | 'NA'
Values in _sign columns:        ''   (significant, p<0.05)
                                'ns' (trend, not significant)

Usage:
    python panel_A_chromatin_expression.py
    python panel_A_chromatin_expression.py --csv my_table.csv --out my_figure.png
"""

import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_CSV = "table_A_chromatin_expression.csv"
DEFAULT_OUT = "panel_A_chromatin_expression.png"

# Cancer order: main-figure cancers first, then supplement
MAIN9  = ['BRCA','COAD','HNSC','KIRC','LIHC','LUAD','LUSC','PRAD','THCA']
REST24 = ['ACC','BLCA','CESC','CHOL','DLBC','ESCA','GBM','KICH','KIRP',
          'LAML','LGG','MESO','OV','PAAD','PCPG','READ','SARC','SKCM',
          'STAD','TGCT','THYM','UCEC','UCS','UVM']
ALL33  = MAIN9 + REST24

REGS  = ['H2AFZ','ATAD2','BRD4','BRD2']
YEAST = {'H2AFZ':'HTZ1','ATAD2':'YTA7','BRD4':'BDF1','BRD2':'BDF1'}

# Colours
UP_SIG   = '#c0392b'   # red  — significant upregulation
DN_SIG   = '#2166ac'   # blue — significant downregulation
GREY     = '#b8b8b8'   # grey — trend, not significant
NO_CHG_MAIN = '#c0c0c0'
NO_CHG_REST = '#dedede'

# Arrow geometry
SHAFT  = 0.28   # arrow length (data units)
HEAD_W = 0.28   # arrowhead width
HEAD_L = 0.16   # arrowhead length
LW     = 2.5    # shaft line width

# Spacing
DOT_SPACING = 1.05
GAP         = 1.8    # gap between main9 and rest24


def load_data(csv_path):
    rows = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return {r['cancer']: r for r in rows}


def draw_arrow(ax, x, y, direction, color):
    """Draw a short, fat arrow at (x, y) pointing up or down."""
    if direction == 'up':
        ax.annotate('',
            xy=(x, y + SHAFT), xytext=(x, y),
            arrowprops=dict(
                arrowstyle=f'->,head_width={HEAD_W},head_length={HEAD_L}',
                color=color, lw=LW))
    elif direction == 'down':
        ax.annotate('',
            xy=(x, y - SHAFT), xytext=(x, y),
            arrowprops=dict(
                arrowstyle=f'->,head_width={HEAD_W},head_length={HEAD_L}',
                color=color, lw=LW))


def make_figure(csv_path, out_path, dpi=180):
    data = load_data(csv_path)

    # x positions
    dot_xs_main = np.arange(9)  * DOT_SPACING
    dot_xs_rest = np.arange(24) * DOT_SPACING + 9 * DOT_SPACING + GAP
    dot_xs      = np.concatenate([dot_xs_main, dot_xs_rest])
    sep_x       = 9 * DOT_SPACING + GAP / 2 - 0.1
    x_max       = dot_xs[-1] + 1.0

    fig, ax = plt.subplots(figsize=(22, 7.5))
    fig.patch.set_facecolor('#f9f9f7')
    ax.set_facecolor('#f9f9f7')
    ax.set_xlim(-2.2, x_max + 0.3)
    ax.set_ylim(-2.2, 4.6)
    ax.axis('off')

    # ── Titles ────────────────────────────────────────────────────────────────
    ax.text(x_max / 2, 4.52,
            'Chromatin regulator expression changes in cancer vs. normal tissue',
            ha='center', fontsize=13, fontweight='bold', color='#1a1a2e')
    ax.text(x_max / 2, 4.25,
            '33 TCGA cancer types  ·  GEPIA2  ·  '
            'red/blue = significant (p<0.05)  ·  grey = trend (ns)',
            ha='center', fontsize=9, color='#888', style='italic')

    # ── Section labels ────────────────────────────────────────────────────────
    ax.text(dot_xs_main.mean(), 3.95, 'Main figure cancers',
            ha='center', fontsize=9, fontweight='bold', color='#444',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='#e8e8e8',
                      edgecolor='#aaa', linewidth=0.8))
    ax.text(dot_xs_rest.mean(), 3.95, 'Additional cancers (supplement)',
            ha='center', fontsize=8.5, color='#999', style='italic',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='#f2f2f2',
                      edgecolor='#ccc', linewidth=0.8))

    # ── Separator ─────────────────────────────────────────────────────────────
    ax.axvline(sep_x, color='#bbb', lw=1.5, ls='--', alpha=0.6,
               ymin=0.05, ymax=0.92)

    # ── One row per regulator ─────────────────────────────────────────────────
    for ri, reg in enumerate(REGS):
        y = 2.9 - ri * 1.0

        # Compute counts for summary label
        r_vals = [(data[c][reg].strip(), data[c][reg + '_sign'].strip())
                  for c in ALL33 if c in data]
        n_up_s = sum(1 for d, s in r_vals if d == 'up'   and s == '')
        n_up_g = sum(1 for d, s in r_vals if d == 'up'   and s == 'ns')
        n_dn_s = sum(1 for d, s in r_vals if d == 'down' and s == '')
        n_dn_g = sum(1 for d, s in r_vals if d == 'down' and s == 'ns')

        count = f'↑{n_up_s}' + (f'+{n_up_g}ns' if n_up_g else '')
        if n_dn_s or n_dn_g:
            count += f'  ↓{n_dn_s}' + (f'+{n_dn_g}ns' if n_dn_g else '')
        count += '/33'

        # Gene label (black)
        ax.text(-0.3, y + 0.10, f'{reg}  ({YEAST[reg]})',
                fontsize=11, fontweight='bold', color='#1a1a1a',
                va='center', ha='right')
        ax.text(-0.3, y - 0.28, count,
                fontsize=8.5, color='#555', va='center',
                ha='right', style='italic')

        # Arrows / dashes for each cancer
        for ci, cancer in enumerate(ALL33):
            if cancer not in data:
                continue
            x         = dot_xs[ci]
            direction = data[cancer][reg].strip()
            sign      = data[cancer][reg + '_sign'].strip()

            if direction == 'NA':
                ax.text(x, y, '?', ha='center', va='center',
                        fontsize=7, color='#ccc')
            elif direction == 'up':
                col = UP_SIG if sign == '' else GREY
                draw_arrow(ax, x, y, 'up', col)
            elif direction == 'down':
                col = DN_SIG if sign == '' else GREY
                draw_arrow(ax, x, y, 'down', col)
            else:
                col = NO_CHG_MAIN if ci < 9 else NO_CHG_REST
                ax.plot(x, y, '_', ms=8, color=col, lw=1.5, zorder=2)

    # ── Cancer name labels ────────────────────────────────────────────────────
    for ci, cancer in enumerate(ALL33):
        x  = dot_xs[ci]
        fc = '#1a1a1a' if ci < 9 else '#666666'
        fs = 8.5       if ci < 9 else 7.5
        ax.text(x, -0.60, cancer, fontsize=fs, ha='right', va='top',
                rotation=50, color=fc, fontweight='bold',
                fontfamily='monospace')

    # ── Legend ────────────────────────────────────────────────────────────────
    lx = dot_xs_main[0]
    ly = -1.72

    draw_arrow(ax, lx,       ly, 'up',   UP_SIG)
    ax.text(lx + 0.35,  ly + 0.10, 'Significant up',   fontsize=8.5, color='#555', va='center')

    draw_arrow(ax, lx + 3.5, ly, 'up',   GREY)
    ax.text(lx + 3.85,  ly + 0.10, 'Trend up (ns)',     fontsize=8.5, color='#888', va='center')

    ax.plot(lx + 8, ly, '_', ms=8, color='#c0c0c0', lw=1.5)
    ax.text(lx + 8.4,   ly,        'No change',         fontsize=8.5, color='#aaa', va='center')

    draw_arrow(ax, lx + 12.5, ly, 'down', DN_SIG)
    ax.text(lx + 12.85, ly - 0.10, 'Significant down',  fontsize=8.5, color='#555', va='center')

    draw_arrow(ax, lx + 17,   ly, 'down', GREY)
    ax.text(lx + 17.35, ly - 0.10, 'Trend down (ns)',   fontsize=8.5, color='#888', va='center')

    plt.tight_layout(pad=0.5)
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', facecolor='#f9f9f7')
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate Panel A — chromatin expression figure.')
    parser.add_argument('--csv', default=DEFAULT_CSV,
                        help=f'Input CSV (default: {DEFAULT_CSV})')
    parser.add_argument('--out', default=DEFAULT_OUT,
                        help=f'Output PNG (default: {DEFAULT_OUT})')
    parser.add_argument('--dpi', type=int, default=180,
                        help='Output DPI (default: 180)')
    args = parser.parse_args()
    make_figure(args.csv, args.out, args.dpi)
