# Porphyrin_2026

Computational analysis code and data accompanying:

> **Chromatin regulator perturbation causes heme pathway imbalance in yeast and human cells**  
> Mitkevich OV, Kushnirov VV, Agaphonov MO, Alexandrov AI  
> *The FEBS Journal* (2026), manuscript ID FJ-26-0956

---

## What this repository contains

This repository provides everything needed to reproduce the GEPIA2-based co-expression
analyses underlying Figure 4 of the manuscript (panels B and C), and to regenerate
all figure panels.

**Scripts** (`scripts/`)

| File | Purpose |
|------|---------|
| `run_gepia_pdf_batch3.py` | Batch GEPIA2 correlation queries — downloads one PDF per gene pair × dataset |
| `build_figure_data.py` | Aggregates raw query results into the figure-ready CSVs |
| `figure4_plot2.py` | Generates Figure 4 (panels A, B, C) |
| `panel_A_chromatin_expression.py` | Generates Supplemental Figure panel: chromatin regulator expression changes across 33 TCGA cancers |
| `panel_B_chromatin_expression.py` | Generates Supplemental Figure: combined chromatin + HEM expression + co-expression heatmap panels |
| `render_fig5.py` | Renders the Figure 5 model schematic (SVG → PNG/PDF) |

**Note:** PDF value extraction from GEPIA2 output was performed with `extract_gepia_values.R`
(using the `pdftools` R package). This script is referenced in `build_figure_data.py`
and was used to produce the `results_ph*.csv` files provided here. The R script itself
is not included in this archive — the pre-extracted CSVs are provided instead, so
re-running the R step is not required to reproduce the figures.

**Raw GEPIA2 query results** (pre-extracted, ready for figure generation)

| File | Contents |
|------|---------|
| `results_ph1.csv` | Phase 1: Spearman R and p-values for all chromatin regulator × heme gene pairs across 23 TCGA matched normal tissue types |
| `results_ph2.csv` | Phase 2: Same pairs across all 33 TCGA tumor datasets (pan-cancer) |
| `results_ph3.csv` | Phase 3: Per-cancer correlation analysis for individual tumor types |
| `phase5.csv` | Duplicate of `results_ph2.csv` (retained for legacy compatibility with `panel_B_chromatin_expression.py`) |
| `phase6.csv` | Duplicate of `results_ph1.csv` (retained for legacy compatibility with `panel_B_chromatin_expression.py`) |

**Figure-ready data files**

| File | Used in |
|------|---------|
| `figure_4A_mammalian_KD_data.csv` | Figure 4A — manually curated from published mammalian knockdown datasets |
| `figure_4B_coexpression_normal_data_NEW.csv` | Figure 4B — chromatin regulator × heme gene Spearman R, averaged across 23 normal tissues |
| `figure_4C_within_pathway_data_NEW.csv` | Figure 4C — within-pathway consecutive-enzyme correlations, normal vs cancer |
| `table_A_chromatin_expression.csv` | Supplemental — chromatin regulator expression direction per cancer type |
| `table_B_HEM_expression.csv` | Supplemental — HEM enzyme expression direction per cancer type |

---

## Analysis design

Data were obtained from [GEPIA2](https://gepia2.cancer-pku.cn) (Tang et al. 2019), which
serves RNA-seq data from TCGA. All analyses used log2(TPM+1) values and the Spearman
correlation method.

Four analysis phases were run:

| Phase | Dataset | Purpose |
|-------|---------|---------|
| 1 | 23 TCGA matched normal tissue types | Normal-tissue baseline co-expression (Figure 4B) |
| 2 | All 33 TCGA tumor datasets | Pan-cancer co-expression |
| 3 | All TCGA tumor datasets, per-cancer | Per-cancer correlation maps |
| 4 | Focused cancers (GBM, BLCA, ESCA, HNSC) | Deep per-cancer analysis |

**Genes analysed:**

- Chromatin regulators: H2AFZ, ATAD2, BRD4, BRD2
- Heme biosynthesis pathway: ALAS1, ALAD, HMBS, UROS, UROD, CPOX, PPOX, FECH
- Confound controls: MKI67, HDAC1 (included in Phases 1 and 2)
- Positive controls: E2F1/CCNE1, CDK2/CCNE1 (included in Phase 2)

**Normal tissue baseline:** TCGA matched normals (biopsy samples from the non-tumour
tissue of cancer patients), not GTEx. Only tissue types with ≥40 samples were included
(23 types). Figure 4B values are unweighted arithmetic means of Spearman R across all
23 tissue types; the minimum p-value across tissues is reported.

---

## How to reproduce the figures (quick start)

The figure-ready CSVs are already in the repository, so you can go directly to
step 3 without re-querying GEPIA2.

### Step 1 — Install dependencies

```bash
pip install gepia pandas numpy matplotlib openpyxl
```

For Figure 5 rendering only:
```bash
pip install cairosvg
# Linux: sudo apt-get install fonts-liberation   (metric-compatible Arial substitute)
```

### Step 2 — (Optional) Re-run GEPIA2 queries

> ⚠️ This requires a working network connection to the GEPIA2 server
> (`gepia2.cancer-pku.cn` / `118.190.148.166`). **This must be run on your own
> machine or cluster**, not in sandboxed environments such as Claude.ai.
> The `gepia` package queries one dataset per request; multi-dataset pooled queries
> are not supported by the server.

```bash
# Phase 1 (TCGA matched normals):
python scripts/run_gepia_pdf_batch3.py --phase 1 --sleep 8

# Phase 2 (all TCGA tumors):
python scripts/run_gepia_pdf_batch3.py --phase 2 --sleep 8

# All phases at once:
python scripts/run_gepia_pdf_batch3.py --all --sleep 10

# Resume an interrupted run:
python scripts/run_gepia_pdf_batch3.py --all --resume --sleep 8
```

Downloaded PDFs are saved to `./gepia_pdfs/`. After downloading, extract values
with `extract_gepia_values.R` (requires R and the `pdftools` package), then run
`build_figure_data.py` to rebuild the `*_NEW.csv` files.

### Step 3 — Generate Figure 4

```bash
python scripts/figure4_plot2.py
```

Outputs: `figure4_combined.png` (300 dpi) and `figure4_combined.svg` (vector).

Input files expected in the working directory:
- `figure_4A_mammalian_KD_data.csv`
- `figure_4B_coexpression_normal_data_NEW.csv`
- `figure_4C_within_pathway_data_NEW.csv`

### Step 4 — Generate Supplemental cancer expression figures

```bash
# Panel A (chromatin regulators):
python scripts/panel_A_chromatin_expression.py

# Combined supplemental figure (chromatin + HEM + co-expression):
python scripts/panel_B_chromatin_expression.py
```

Input files: `table_A_chromatin_expression.csv`, `table_B_HEM_expression.csv`,
`phase5.csv`, `phase6.csv` (all in working directory).

### Step 5 — Render Figure 5 model schematic

```bash
python scripts/render_fig5.py
```

Expects `fig5_chromatin_iron_heme_model.svg` in the working directory.

---

## Data provenance

- **GEPIA2 co-expression values** were retrieved programmatically via the `gepia`
  Python package (version ≥0.3). Extracted R and p-values were spot-checked against
  the GEPIA2 web interface for a subset of gene pairs to confirm extraction accuracy.
- **Mammalian knockdown data** (Figure 4A) were manually curated from published
  GEO datasets: GSE131579 (H2AFZ knockdown, Lamaa et al. 2020), GSE70314 (Atad2
  siRNA in mouse ESCs, Morozumi et al. 2016), and GSE328851 (Atad2 shRNA in C2C12
  myoblasts, Sedraoui et al. 2025).
- **Cancer/normal expression direction tables** (`table_A/B_*.csv`) are hand-corrected
  summaries of GEPIA2 differential expression output.

---

## AI disclosure

Python and R scripts for GEPIA2 querying, PDF-based value extraction, and figure
generation were developed with the assistance of **Claude Sonnet 4.6** (Anthropic),
accessed via the claude.ai web interface. Three coding tasks were performed:
(1) automated GEPIA2 API querying, (2) extraction of numerical values from output PDFs,
and (3) figure generation. All analytical logic, gene and dataset selections, and
statistical interpretation were specified and verified by the authors. Prompts used
are provided in Supplemental File 1 of the manuscript.

---

## Citation

If you use this code or data, please cite:

> Mitkevich OV, Kushnirov VV, Agaphonov MO, Alexandrov AI.
> Chromatin regulator perturbation causes heme pathway imbalance in yeast and human cells.
> *The FEBS Journal* (2026). Manuscript ID: FJ-26-0956.

---

## Contact

Correspondence: aleksandr.aleksandrov@unige.ch
