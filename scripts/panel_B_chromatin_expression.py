"""
cancer_figure_combined.py
--------------------------
Three-panel combined cancer figure (horizontal layout):

  A  — Chromatin regulator expression changes across 33 TCGA cancer types
  B  — HEM pathway enzyme expression changes across 33 TCGA cancer types
  C  — Co-expression heatmaps (matched normal vs cancer, 9 tissue types):
         C-top:    HEM consecutive-enzyme co-expression (normal | cancer | delta)
         C-bottom: Chromatin regulator × HEM enzyme co-expression (normal | cancer)

Inputs (all expected in the same directory unless --data-dir is specified):
    table_A_chromatin_expression.csv   hand-corrected chromatin regulator table
    table_B_HEM_expression.csv         hand-corrected HEM enzyme table
    phase5.csv                         GEPIA2 correlations — cancer tissues
    phase6.csv                         GEPIA2 correlations — individual normals

Output:
    cancer_figure_combined.png

CSV column conventions (panels A and B):
    cancer, in_main_figure,
    GENE, GENE_sign, ...
  Values:  'up' | 'down' | 'ns' | 'no' | 'NA'
  _sign:   '' = significant (p<0.05)   'ns' = trend, not significant

Usage:
    python cancer_figure_combined.py
    python cancer_figure_combined.py --data-dir /path/to/data --out figure.png --dpi 200
"""

import argparse
import csv
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgridspec
import numpy as np
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

# ── Cancer sets ───────────────────────────────────────────────────────────────
MAIN9  = ['BRCA','COAD','HNSC','KIRC','LIHC','LUAD','LUSC','PRAD','THCA']
REST24 = ['ACC','BLCA','CESC','CHOL','DLBC','ESCA','GBM','KICH','KIRP',
          'LAML','LGG','MESO','OV','PAAD','PCPG','READ','SARC','SKCM',
          'STAD','TGCT','THYM','UCEC','UCS','UVM']
ALL33  = MAIN9 + REST24

NORMAL_N = {'BRCA':114,'COAD':41,'HNSC':44,'KIRC':72,'LIHC':50,
            'LUAD':59,'LUSC':49,'PRAD':52,'THCA':59}

# ── Gene sets ─────────────────────────────────────────────────────────────────
CHROM_REGS = ['H2AFZ','ATAD2','BRD4','BRD2']
HEM_GENES  = ['HMBS','UROS','UROD','CPOX','PPOX','FECH']
HEM_ALL    = ['ALAS1','ALAD','HMBS','UROS','UROD','CPOX','PPOX','FECH']

YEAST_C = {'H2AFZ':'HTZ1','ATAD2':'YTA7','BRD4':'BDF1','BRD2':'BDF1'}
YEAST_H = {'HMBS':'Hem3','UROS':'Hem4','UROD':'Hem12',
            'CPOX':'Hem13','PPOX':'Hem14','FECH':'Hem15'}

HEM_PAIRS   = [('ALAS1','ALAD'),('ALAD','HMBS'),('HMBS','UROS'),
               ('UROS','UROD'),('UROD','CPOX'),('CPOX','PPOX'),('PPOX','FECH')]
PAIR_LABELS = ['ALAS1\n→ALAD','ALAD\n→HMBS','HMBS\n→UROS','UROS\n→UROD',
               'UROD\n→CPOX','CPOX\n→PPOX','PPOX\n→FECH']

# ── Colours & arrow geometry ──────────────────────────────────────────────────
UP_SIG = '#c0392b'   # red  — significant up
DN_SIG = '#2166ac'   # blue — significant down
GREY   = '#b8b8b8'   # grey — trend, not significant
SHAFT  = 0.28
HEAD_W = 0.28
HEAD_L = 0.16
LW     = 2.5

# ── Dot-plot spacing ──────────────────────────────────────────────────────────
DOT_SP = 1.05
GAP    = 1.8


# ── Data loaders ─────────────────────────────────────────────────────────────
def load_expr(path):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return {r['cancer']: r for r in rows}


def load_corr(path):
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def get_r(rows, g1, g2, ds):
    for r in rows:
        if (r['gene1'] == g1 and r['gene2'] == g2
                and r['dataset'] == ds and r['status'] == 'ok'):
            try:
                return float(r['R'])
            except (ValueError, TypeError):
                return None
    return None


# ── Dot-plot helpers ──────────────────────────────────────────────────────────
def make_dot_xs():
    main = np.arange(9)  * DOT_SP
    rest = np.arange(24) * DOT_SP + 9 * DOT_SP + GAP
    return np.concatenate([main, rest]), main, rest


def draw_arrow(ax, x, y, direction, color):
    kw = dict(arrowstyle=f'->,head_width={HEAD_W},head_length={HEAD_L}',
              color=color, lw=LW)
    if direction == 'up':
        ax.annotate('', xy=(x, y + SHAFT), xytext=(x, y), arrowprops=kw)
    elif direction == 'down':
        ax.annotate('', xy=(x, y - SHAFT), xytext=(x, y), arrowprops=kw)


def draw_dot_row(ax, y, gene, data_dict, dot_xs, label='', sublabel='', note=''):
    sign_col = gene + '_sign'
    ax.text(-0.3, y + 0.10, label, fontsize=10.5, fontweight='bold',
            color='#1a1a1a', va='center', ha='right')
    ax.text(-0.3, y - 0.28, sublabel, fontsize=8, color='#555',
            va='center', ha='right', style='italic')
    if note:
        ax.text(-0.3, y - 0.56, note, fontsize=7.5, color=DN_SIG,
                va='center', ha='right', style='italic', fontweight='bold')
    for ci, cancer in enumerate(ALL33):
        if cancer not in data_dict:
            continue
        x = dot_xs[ci]
        d = data_dict[cancer][gene].strip().lower()
        s = data_dict[cancer].get(sign_col, '').strip()
        if d == 'na':
            ax.text(x, y, '?', ha='center', va='center', fontsize=6, color='#ccc')
        elif d == 'up':
            draw_arrow(ax, x, y, 'up',   UP_SIG if s == '' else GREY)
        elif d == 'down':
            draw_arrow(ax, x, y, 'down', DN_SIG if s == '' else GREY)
        else:
            col = '#c0c0c0' if ci < 9 else '#dedede'
            ax.plot(x, y, '_', ms=7, color=col, lw=1.4, zorder=2)


def cancer_labels(ax, dot_xs, y_pos):
    for ci, cancer in enumerate(ALL33):
        fc = '#1a1a1a' if ci < 9 else '#777'
        fs = 8.0       if ci < 9 else 7.0
        ax.text(dot_xs[ci], y_pos, cancer, fontsize=fs, ha='right', va='top',
                rotation=50, color=fc, fontweight='bold', fontfamily='monospace')


def section_boxes(ax, dot_xs_main, dot_xs_rest, y):
    ax.text(dot_xs_main.mean(), y, 'Main figure cancers',
            ha='center', fontsize=8.5, fontweight='bold', color='#444',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='#e8e8e8',
                      edgecolor='#aaa', linewidth=0.8))
    ax.text(dot_xs_rest.mean(), y, 'Additional cancers (supplement)',
            ha='center', fontsize=8, color='#999', style='italic',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='#f2f2f2',
                      edgecolor='#ccc', linewidth=0.8))


def vsep(ax, sep_x, ymin=0.04, ymax=0.93):
    ax.axvline(sep_x, color='#bbb', lw=1.5, ls='--', alpha=0.6,
               ymin=ymin, ymax=ymax)


def dot_legend(ax, lx, ly):
    draw_arrow(ax, lx,       ly, 'up',   UP_SIG)
    ax.text(lx + 0.35,  ly + 0.10, 'Sig. up',   fontsize=8, color='#555', va='center')
    draw_arrow(ax, lx + 3.0,  ly, 'up',   GREY)
    ax.text(lx + 3.35,  ly + 0.10, 'Trend up',  fontsize=8, color='#888', va='center')
    ax.plot(lx + 6.5, ly, '_', ms=7, color='#c0c0c0', lw=1.4)
    ax.text(lx + 6.9,   ly,        'No change', fontsize=8, color='#aaa', va='center')
    draw_arrow(ax, lx + 10.5, ly, 'down', DN_SIG)
    ax.text(lx + 10.85, ly - 0.10, 'Sig. down', fontsize=8, color='#555', va='center')
    draw_arrow(ax, lx + 14.0, ly, 'down', GREY)
    ax.text(lx + 14.35, ly - 0.10, 'Trend dn',  fontsize=8, color='#888', va='center')


def count_str(gene, data_dict):
    r_vals = [(data_dict[c][gene].strip(),
               data_dict[c].get(gene + '_sign', '').strip())
              for c in ALL33 if c in data_dict]
    n_up_s = sum(1 for d, s in r_vals if d == 'up'   and s == '')
    n_up_g = sum(1 for d, s in r_vals if d == 'up'   and s == 'ns')
    n_dn_s = sum(1 for d, s in r_vals if d == 'down' and s == '')
    n_dn_g = sum(1 for d, s in r_vals if d == 'down' and s == 'ns')
    parts = []
    if n_up_s or n_up_g:
        parts.append(f'↑{n_up_s}' + (f'+{n_up_g}ns' if n_up_g else ''))
    if n_dn_s or n_dn_g:
        parts.append(f'↓{n_dn_s}' + (f'+{n_dn_g}ns' if n_dn_g else ''))
    if not parts:
        parts = ['ns all']
    return '  '.join(parts) + '/33'


# ── Heatmap helper ────────────────────────────────────────────────────────────
def draw_heatmap(ax, mat, cmap, norm, xlabels, ylabels, show_y,
                 title, tc, fs_cell=7, fs_x=7.5, fs_y=8.5):
    ax.imshow(mat, cmap=cmap, norm=norm, aspect='auto')
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=0, ha='center',
                        fontsize=fs_x, fontfamily='monospace', fontweight='bold')
    ax.set_yticks(range(len(MAIN9)))
    ax.set_yticklabels(ylabels if show_y else [], fontsize=fs_y, fontweight='bold')
    ax.tick_params(length=0)
    if title:
        ax.set_title(title, fontsize=10, fontweight='bold', color=tc, pad=5)
    for sp in ax.spines.values():
        sp.set_linewidth(1.6)
        sp.set_color(tc)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):
                continue
            fc = 'white' if abs(v) > 0.40 else '#1a1a1a'
            ax.text(j, i, f'{v:+.2f}', ha='center', va='center',
                    fontsize=fs_cell, color=fc, fontweight='bold')


# ── Correlation matrix builders ───────────────────────────────────────────────
def build_heme_mats(rows5, rows6):
    n = np.full((len(MAIN9), 7), np.nan)
    c = np.full((len(MAIN9), 7), np.nan)
    for i, cancer in enumerate(MAIN9):
        for j, (g1, g2) in enumerate(HEM_PAIRS):
            n[i, j] = get_r(rows6, g1, g2, f'{cancer}_Normal') or np.nan
            c[i, j] = get_r(rows5, g1, g2, f'{cancer}_Tumor')  or np.nan
    return n, c


def build_chrom_mat(reg, tissue, rows5, rows6):
    mat = np.full((len(MAIN9), 8), np.nan)
    for i, cancer in enumerate(MAIN9):
        for j, hem in enumerate(HEM_ALL):
            if tissue == 'normal':
                mat[i, j] = get_r(rows6, reg, hem, f'{cancer}_Normal') or np.nan
            else:
                mat[i, j] = get_r(rows5, reg, hem, f'{cancer}_Tumor')  or np.nan
    return mat


# ── Main figure function ──────────────────────────────────────────────────────
def make_figure(data_dir, out_path, dpi=150):
    # Load data
    chrom_data = load_expr(os.path.join(data_dir, 'table_A_chromatin_expression.csv'))
    hem_data   = load_expr(os.path.join(data_dir, 'table_B_HEM_expression.csv'))
    rows5 = load_corr(os.path.join(data_dir, 'phase5.csv'))
    rows6 = load_corr(os.path.join(data_dir, 'phase6.csv'))

    # Precompute correlation matrices
    heme_n, heme_c = build_heme_mats(rows5, rows6)
    heme_d  = heme_c - heme_n
    ylabels = [f"{c} (n={NORMAL_N[c]})" for c in MAIN9]

    # Colour maps
    cmap  = LinearSegmentedColormap.from_list('bwr2',
            ['#2166ac','#4393c3','#92c5de','#f7f7f7',
             '#f4a582','#d6604d','#b2182b'])
    cnorm = TwoSlopeNorm(vmin=-0.75, vcenter=0, vmax=0.75)
    dnorm = TwoSlopeNorm(vmin=-0.65, vcenter=0, vmax=0.65)

    # Dot positions
    dot_xs, dot_xs_main, dot_xs_rest = make_dot_xs()
    sep_x = 9 * DOT_SP + GAP / 2 - 0.1
    x_max = dot_xs[-1] + 1.0

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(48, 26))
    fig.patch.set_facecolor('#f9f9f7')

    main = mgridspec.GridSpec(1, 2,
                               width_ratios=[1.0, 1.0],
                               wspace=0.14,
                               left=0.05, right=0.98,
                               top=0.955, bottom=0.03)

    left = mgridspec.GridSpecFromSubplotSpec(2, 1,
        subplot_spec=main[0],
        height_ratios=[1, 1.5], hspace=0.10)

    # ── Panel A ───────────────────────────────────────────────────────────────
    ax_a = fig.add_subplot(left[0])
    ax_a.set_xlim(-2.2, x_max + 0.3)
    ax_a.set_ylim(-2.1, 4.5)
    ax_a.axis('off')

    ax_a.text(-0.04, 1.02, 'A', transform=ax_a.transAxes,
              fontsize=22, fontweight='bold', va='top', color='#222')
    ax_a.text(x_max / 2, 4.38,
              'Chromatin regulator expression in cancer vs. normal tissue  '
              '(33 TCGA types, GEPIA2)',
              ha='center', fontsize=11.5, fontweight='bold', color='#1a1a2e')

    section_boxes(ax_a, dot_xs_main, dot_xs_rest, 4.00)
    vsep(ax_a, sep_x)

    for ri, reg in enumerate(CHROM_REGS):
        y = 3.0 - ri * 0.93
        draw_dot_row(ax_a, y, reg, chrom_data, dot_xs,
                     label=f'{reg}  ({YEAST_C[reg]})',
                     sublabel=count_str(reg, chrom_data))

    cancer_labels(ax_a, dot_xs, -0.58)
    dot_legend(ax_a, dot_xs_main[0], -1.70)

    # ── Panel B ───────────────────────────────────────────────────────────────
    ax_b = fig.add_subplot(left[1])
    ax_b.set_xlim(-2.2, x_max + 0.3)
    ax_b.set_ylim(-2.1, 6.3)
    ax_b.axis('off')

    ax_b.text(-0.04, 1.02, 'B', transform=ax_b.transAxes,
              fontsize=22, fontweight='bold', va='top', color='#222')
    ax_b.text(x_max / 2, 6.18,
              'HEM pathway enzyme expression in cancer vs. normal tissue  '
              '(33 TCGA types, GEPIA2)',
              ha='center', fontsize=11.5, fontweight='bold', color='#1a1a2e')

    section_boxes(ax_b, dot_xs_main, dot_xs_rest, 5.80)
    vsep(ax_b, sep_x)

    ax_b.annotate('', xy=(dot_xs_main[-1] + 0.2, 5.52),
                  xytext=(dot_xs_main[0] - 0.1, 5.52),
                  arrowprops=dict(arrowstyle='->', color='#bbb', lw=1.1))
    ax_b.text((dot_xs_main[0] + dot_xs_main[-1]) / 2, 5.38,
              'heme biosynthesis order →',
              ha='center', fontsize=7.5, color='#ccc', style='italic')

    for gi, gene in enumerate(HEM_GENES):
        y = 4.6 - gi * 1.06
        note = 'never upregulated' if gene == 'FECH' else ''
        draw_dot_row(ax_b, y, gene, hem_data, dot_xs,
                     label=f'{gene}  ({YEAST_H[gene]})',
                     sublabel=count_str(gene, hem_data), note=note)

    cancer_labels(ax_b, dot_xs, -0.58)
    dot_legend(ax_b, dot_xs_main[0], -1.70)

    # ── Panel C ───────────────────────────────────────────────────────────────
    right = mgridspec.GridSpecFromSubplotSpec(2, 1,
        subplot_spec=main[1],
        height_ratios=[1, 2.2], hspace=0.12)

    c_pos = main[1].get_position(fig)
    fig.text(c_pos.x0 - 0.008, c_pos.y1 + 0.005, 'C',
             fontsize=22, fontweight='bold', color='#222', va='bottom')
    fig.text(c_pos.x0 + 0.008, c_pos.y1 + 0.005,
             'Co-expression analysis — matched normal vs. cancer  '
             '(9 tissue types, n≥40, Spearman R)',
             fontsize=11, fontweight='bold', color='#1a1a2e', va='bottom')

    # C-top: within-pathway
    ct     = mgridspec.GridSpecFromSubplotSpec(1, 3,
        subplot_spec=right[0], wspace=0.04, width_ratios=[7, 7, 7])
    ct_pos = right[0].get_position(fig)
    fig.text(ct_pos.x0, ct_pos.y1 + 0.003,
             'HEM consecutive-enzyme co-expression',
             fontsize=9.5, fontweight='bold', color='#555', va='bottom')

    for col, (mat, norm, title, tc) in enumerate([
            (heme_n, cnorm, 'Normal tissue',        '#2c7bb6'),
            (heme_c, cnorm, 'Cancer',               '#c0392b'),
            (heme_d, dnorm, 'Delta (Cancer−Normal)', '#555555')]):
        ax = fig.add_subplot(ct[col])
        draw_heatmap(ax, mat, cmap, norm, PAIR_LABELS,
                     ylabels if col == 0 else [], col == 0,
                     title, tc, fs_cell=8)

    cbar1 = fig.add_axes([ct_pos.x1 + 0.003,
                           ct_pos.y0 + ct_pos.height * 0.52,
                           0.005, ct_pos.height * 0.43])
    cb1 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=cnorm), cax=cbar1)
    cb1.set_label('Spearman R', fontsize=7)
    cb1.ax.tick_params(labelsize=6.5)

    cbar2 = fig.add_axes([ct_pos.x1 + 0.003,
                           ct_pos.y0 + 0.005,
                           0.005, ct_pos.height * 0.43])
    cb2 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=dnorm), cax=cbar2)
    cb2.set_label('ΔR', fontsize=7)
    cb2.ax.tick_params(labelsize=6.5)

    # C-bottom: chromatin × HEM
    cb_gs  = mgridspec.GridSpecFromSubplotSpec(4, 2,
        subplot_spec=right[1], hspace=0.06, wspace=0.04,
        width_ratios=[1, 1])
    cb_pos = right[1].get_position(fig)
    fig.text(cb_pos.x0, cb_pos.y1 + 0.003,
             'Chromatin regulator × HEM enzyme co-expression',
             fontsize=9.5, fontweight='bold', color='#555', va='bottom')

    for ri, reg in enumerate(CHROM_REGS):
        for col, (tissue, title, tc) in enumerate([
                ('normal', 'Normal tissue', '#2c7bb6'),
                ('tumor',  'Cancer',        '#c0392b')]):
            mat = build_chrom_mat(reg, tissue, rows5, rows6)
            ax  = fig.add_subplot(cb_gs[ri, col])
            draw_heatmap(ax, mat, cmap, cnorm, HEM_ALL,
                         ylabels if col == 0 else [], col == 0,
                         title if ri == 0 else '', tc,
                         fs_cell=6.5, fs_x=7, fs_y=8)
            if col == 0:
                ax.set_ylabel(f'{reg}\n({YEAST_C[reg]})', fontsize=9,
                             fontweight='bold', color='#1a1a1a', labelpad=5)

    cbar3 = fig.add_axes([cb_pos.x1 + 0.003,
                           cb_pos.y0 + 0.01,
                           0.005, cb_pos.height * 0.92])
    cb3 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=cnorm), cax=cbar3)
    cb3.set_label('Spearman R', fontsize=7)
    cb3.ax.tick_params(labelsize=6.5)

    # ── Main title ────────────────────────────────────────────────────────────
    fig.text(0.5, 0.972,
             'Chromatin dysregulation and heme pathway imbalance in cancer',
             ha='center', fontsize=16, fontweight='bold',
             fontfamily='DejaVu Serif', color='#1a1a2e')

    plt.savefig(out_path, dpi=dpi, bbox_inches='tight', facecolor='#f9f9f7')
    plt.close()
    print(f"Saved: {out_path}")


# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate combined three-panel cancer figure.')
    parser.add_argument('--data-dir', default='.',
                        help='Directory containing all input CSV files (default: .)')
    parser.add_argument('--out', default='cancer_figure_combined.png',
                        help='Output PNG path (default: cancer_figure_combined.png)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Output DPI (default: 150)')
    args = parser.parse_args()
    make_figure(args.data_dir, args.out, args.dpi)
