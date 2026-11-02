# Handoff: MaizExpressD

Plan for a development counterpart of MaizExpressN (nitrogen) and MaizExpressC (cold): look up a maize gene
and see where and when it is expressed across tissues and developmental stages, how much that varies among
~300 inbreds, whether the variation is genetic (eQTL, tissue by tissue) and whether it tracks agronomic traits.
Owner: Rubén Rellán-Álvarez (GEMMA Lab, NC State). Started 29 Oct 2026 (Design B: requantify the raw reads and
build a sibling of N and C).

## 1. Data

- **Expression: Kremling2018.** Kremling et al. 2018, Nature 555:520 (doi:10.1038/nature25966). 3′ mRNA-seq of
  the 282 diversity panel; SRA PRJNA383416 / SRP115041: 1,960 single-end runs, 532 GB, ~5.5 M reads each.
- **Results layer: Lei2025.** Lei et al. 2025, G3 15:jkaf140 (doi:10.1093/g3journal/jkaf140), a reanalysis of the
  same data (1,771 libraries, 298 lines, B73 v3 gene models). Supplementary Tables S1–S12: lines and
  subpopulations, tissue-specific genes, eQTLs (v3 and v5 positions), TWAS gene–tissue–trait hits, published QTLs
  and their colocalization with eQTLs.

| Stage (order) | Tissues (code, Lei code) |
|---|---|
| 1 Germination | germinating root (GRoot, RT), germinating shoot (GShoot, SH) |
| 2 Seedling, leaf 3 | 2 cm from the base (L3Base, LB), middle (L3Mid, 26 lines only), 2 cm from the tip (L3Tip, LT) |
| 3 Adult plant, after flowering | mature mid-leaf at mid-day (LMAD, LD) and mid-night (LMAN, LN) |
| 4 Kernel | developing kernels, 350 GDD after pollination (Kern, KN) |

## 2. Milestones (⏸ = review checkpoint: stop for Rubén)

| | Milestone | Reuses from MaizExpressN/C | New for development |
|---|---|---|---|
| M1 ⏸ | Sample sheet: runs, libraries, lines matched to Lei's 298 | `runs.tsv` layout | library-name parser, tissue table, line matching; checkpoint: `needs_review.tsv` |
| M2 | Download and quantification on Hazel | downloader, nf-core/rnaseq head job, B73 v5 salmon index (from N) | eight units by tissue; 3′-tag settings from the read-structure check (29 Oct 2026); storage paths approved at M1 |
| M3 ⏸ | Library QC | strict rules, tissue-marker check | B73 libraries as batch controls; tissue-identity check across 8 tissues |
| M4 ⏸ | Analysis | DESeq2 paired diversity-panel fit (C) | expression per tissue across lines, tissue specificity (Lei's method and τ), paired developmental contrasts, variation among lines; checkpoint: `contrasts.tsv` |
| M5 | Annotation and Lei results | `annotation/` as is | Lei tables lifted to v5 gene IDs |
| M6 | Co-expression partners | `coexpr.R`, `prepare_partners.R` | across lines within each tissue |
| M7 ⏸ | Shiny app (first slice for review) | `app/` (gene card, layout) | Traits & eQTL page, Lines page |
| M8 | Tests, repository, deployment | testthat suite, Posit Connect Cloud route | GitHub repo once Rubén chooses visibility |

## 3. Proposed developmental contrasts (to settle at M4)

All paired by line (`~ genotype + tissue`); sign log2(second / first):

- leaf maturation: L3Tip vs L3Base (and L3Mid where present);
- day/night: LMAN vs LMAD, within sampling group;
- seedling vs adult leaf: LMAD vs L3Tip;
- organ: GShoot vs GRoot.

## 4. Layout

Same as MaizExpressN and C: `catalogue/`, `metadata/`, `pipeline/`, `analysis/`, `annotation/`, `results/`,
`app/`, `tests/`, plus `private/` (gitignored).
