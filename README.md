# Chromatin regulation of the heme biosynthesis pathway

Analysis code and data for:

**"Chromatin regulator perturbation causes heme pathway imbalance in yeast and human cells
"**

O.V. Mitkevich, V.V. Kushnirov, M.O. Agaphonov, A.I. Alexandrov

*[Journal and year to be updated upon publication]*

---

## Data files

| File | Figure | Description |
|------|--------|-------------|
| `all_figure_data.xlsx` | All | Master file with all panels as sheets plus metadata |
| `figure_4A_mammalian_KD_data.csv` | Fig 4A | Heme gene log2FC upon H2AFZ/ATAD2 knockdown in mammalian cells |
| `figure_4B_coexpression_normal_data.csv` | Fig 4B | Chromatin regulator vs heme gene Spearman R in normal human tissues |
| `figure_4C_within_pathway_data.csv` | Fig 4C | Within-pathway enzyme pair co-expression, normal vs cancer |
| `supplemental_pancancer_correlations.csv` | Suppl | Pan-cancer chromatin regulator vs heme gene correlations (all TCGA) |

---

## Data sources

All raw data is publicly available:

- **GEPIA2** — gepia2.cancer-pku.cn (Tang et al. 2019, *Nucleic Acids Res*)
  RNA-seq from TCGA (33 cancer types) and GTEx (normal tissues)
- **GSE70314** — Atad2 siRNA knockdown, mouse embryonic stem cells
- **GSE328851** — Atad2 shRNA knockdown, C2C12 mouse myoblasts
- **eLife-53375** — H2AFZ siRNA knockdown, WI38 and U2OS human cells
- **GSE172069** — ATAD2 overexpression, human melanocytes
  (Baggiolini et al. 2021, *Science*)

---

## Analysis scripts

The scripts used to query GEPIA2 and process the mammalian datasets
are available in the full project repository. Contact the corresponding
author (A.I. Alexandrov) for access, or see the Methods section of the
paper for a complete description of the analysis pipeline.

---

## Citation

[To be updated upon publication]

## Funding

Ministry of Science and Higher Education of the Russian Federation
(Federal Scientific and Technical Program for the Development of
Genetic Technologies for 2019–2030, agreement #075-15-2025-470)
