# Sample-sheet spec (M1, development)

MaizExpressD has one expression study, Kremling2018 (BioProject PRJNA383416): 3′ RNA-seq of the maize 282
diversity panel in seven tissues, plus leaf-3 middle. All 1,960 runs share one BioSample, so the design comes
from the library names, `[plate_][submission_]<tissue>_<line>_<barcode>`, for example
`KAKRNA23_10346932_L3Base_A556_ACAGAT` or `10343264_LMAN26_CI21E_ATGAAC`. `build_metadata.py` parses them;
nothing in `metadata/runs.tsv`, `samples.tsv` or `lines.tsv` is edited by hand.

Hand-curated inputs, all in `metadata/build/`:

| File | Content |
|---|---|
| `tissues.tsv` | One row per tissue: the library-name labels it covers, Lei et al. 2025 code, label, group, detail, Plant Ontology term, stage, stage class and order, sampling time, growth system, source |
| `line_aliases.tsv` | Library-name spellings that match a different name in Lei's Table S1, with the kind of match (`spelling`, `alias`, `probable`) and the evidence |
| `label_fixes.tsv` | Tissue labels in library names that the expression contradicts (`run_accession`, `label_in_name`, `label_used`, `evidence`); the sample takes the tissue of `label_used` but keeps the `sample_id` built from the name, since the quantified tables are keyed by it |
| `source/Lei2025_TableS1_lines.tsv` | The 298 lines of Lei et al. 2025 with origin, pedigree and subpopulation (`extract_lei_tables.R`) |
| `ena_runs/PRJNA383416.tsv` | ENA run report (FASTQ URLs, MD5s, read counts, library names, submitted file names) |

Open questions go in `metadata/needs_review.tsv` (`study_id`, `item`, `question`, `proposed_answer`,
`evidence`), written by hand.

## `metadata/samples.tsv`

One row per library (one run each).

| Column | Content |
|---|---|
| sample_id | `Kremling2018_<tissue>_<genotype>_R<replicate>` |
| study_id | `Kremling2018` |
| run_accessions | the run (SRR) |
| library_name | as submitted |
| genotype | canonical line name: the library-name spelling, merged across spellings of one line (the one with most libraries) |
| genotype_sra | the spelling in this library's name |
| genotype_lei | the line's name in Lei Table S1, empty if none matches |
| genotype_class | `B73`, `other inbred`, `landrace` (the MR series), `NR` (B73-10) |
| subpopulation | from Lei Table S1: `ss`, `nss`, `ts`, `sweet`, `popcorn`, `mixed`, `landrace`; `NR` when missing or unmatched |
| origin | from Lei Table S1 |
| tissue | `GRoot`, `GShoot`, `L3Base`, `L3Mid`, `L3Tip`, `LMAD`, `LMAN`, `Kern`, `LMid` (unclear, see review) |
| tissue_lei | Lei's code: `RT`, `SH`, `LB`, `LT`, `LD`, `LN`, `KN`; empty for L3Mid and LMid |
| tissue_label, tissue_group, tissue_detail, po_term | from `tissues.tsv`; groups `Root`, `Shoot`, `Leaf`, `Kernel` |
| stage, stage_class, stage_order | from `tissues.tsv`; classes `germination` (1), `seedling` (2), `adult` (3), `kernel` (4) |
| sampling_time | `mid-day` or `mid-night` for adult leaf, else empty |
| sampling_group | `8` or `26` for adult leaf (from LMAD8, LMAN26...), else empty |
| growth_system | `NR` until filled from Kremling et al. 2018 |
| seq_submission | numeric sequencing-submission ID in the library name, if any |
| plate | `KAKRNA<n>` in the library name, if any |
| flowcell | from the submitted file name |
| barcode | index sequence at the end of the library name |
| replicate | 1, 2, ... within genotype × tissue, by run accession |
| include | `TRUE` for every library at M1; library QC decides later |
| curation_note | library name, how the line was matched, flags |

## `metadata/lines.tsv`

One row per genotype in the library names, then one per Lei line with no libraries: `genotype`,
`genotype_sra_spellings`, `genotype_lei`, `lei_match` (`direct`, `spelling`, `alias`, `probable`, `none`),
`genotype_class`, `subpopulation`, `origin`, `pedigree`, `n_libraries`, `tissues_with_data` (tissue:count),
`adult_leaf_group`, `note`.

## `metadata/runs.tsv`

Same columns as in MaizExpressN and MaizExpressC (so the pipeline scripts read it unchanged), plus
`library_name`. FASTQ URLs are ENA HTTPS links with their MD5s.
