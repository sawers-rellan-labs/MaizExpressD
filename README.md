# MaizExpressD

Look up a maize gene and see where and when it is expressed across tissues and developmental stages
(germinating root and shoot, the base, middle and tip of seedling leaf 3, adult leaf by day and by night,
developing kernels) in about 300 inbred lines of the maize 282 diversity panel, how much that expression varies
among lines, and how it links to eQTLs and agronomic traits. GEMMA Lab, NC State (Rubén Rellán-Álvarez).

Sibling of [MaizExpressN](https://github.com/sawers-rellan-labs/MaizExpressN) (nitrogen) and
[MaizExpressC](https://github.com/sawers-rellan-labs/MaizExpressC) (cold), whose pipeline, gene annotation and app
it reuses. The plan is in `HANDOFF.md`; working notes and every decision taken are in `CLAUDE.md`.

**Status:** milestones 1-8 complete (sample sheet, quantification, library QC, analysis, annotation and Lei
results, co-expression partners, app, tests and deployment). Live app:
<https://01a1ab72-cb13-5f3b-83ac-eed2f6967ef8.share.connect.posit.cloud/>.

Large files (Parquet tables, the per-tissue gene count tables, annotation sources) are stored with
[Git LFS](https://git-lfs.com). Install it before cloning (`brew install git-lfs`, then `git lfs install`), or clone
with `GIT_LFS_SKIP_SMUDGE=1` to get only the code and pointer files. Two large outputs are not in the repository:
the per-line CPM table `results/expression_line.parquet` (rebuilt by `analysis/expression_tissue.R`; the app carries
the same values in `app/data/expr/`) and salmon's TPM tables (not used: CPM comes from the gene counts).

## In numbers

One study (Kremling et al. 2018), 1,960 3′ mRNA-seq libraries requantified on B73 RefGen_v5; 1,836 in use after
library QC, from 296 lines in nine tissues. 12 tissue contrasts (8 paired across lines with edgeR, 4 on B73's
replicate libraries with DESeq2), tissue specificity by Lei et al.'s method and by tau, variation among lines,
co-expression partners across lines in each tissue, and Lei et al. 2025's eQTL, TWAS and QTL colocalization
results on v5 genes.

## Data

- Expression: Kremling et al. 2018, *Nature* 555:520–523 ([doi:10.1038/nature25966](https://doi.org/10.1038/nature25966)),
  3′ mRNA-seq reads from SRA [PRJNA383416](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA383416), requantified here
  on B73 RefGen_v5.
- eQTL, TWAS and QTL colocalization results: Lei et al. 2025, *G3* 15:jkaf140
  ([doi:10.1093/g3journal/jkaf140](https://doi.org/10.1093/g3journal/jkaf140)).

## Run the app

```bash
Rscript -e 'shiny::runApp("app", launch.browser = TRUE)'
```

Run from the project root (the renv library activates there). The app reads `app/data/`, built by
`app/prepare_data.R`; set `MAIZEXPRESSD_DATA` to use a data folder elsewhere. Pages: Gene lookup, Traits & eQTL,
Lines, Co-expression, Rank genes, Gene list, Data, Methods.

## Deploy to Posit Connect Cloud

Live app: <https://01a1ab72-cb13-5f3b-83ac-eed2f6967ef8.share.connect.posit.cloud/> (public link; managed at
<https://connect.posit.cloud/gemmalab/content/01a1ab72-cb13-5f3b-83ac-eed2f6967ef8>, account `gemmalab`).

Same route as MaizExpressN and C. From an R session started in the project root:

```r
rsconnect::connectCloudUser()                                   # once per computer: authorise in the browser
rsconnect::deployApp(appDir = "app", appName = "maizexpressd")   # upload app/ (about 350 MB) and publish
```

`app/data` is about 345 MB (per-line CPM 165, partner lists 125, DE 19, atlas 16), under the 1 GiB bundle limit
of Connect Cloud's free plan. For publishing from GitHub, refresh `app/manifest.json` with
`rsconnect::writeManifest(appDir = "app")` and commit it with the `app/` folder. Re-run `deployApp()` after
`Rscript app/prepare_data.R` to update the app.

## Layout

| Folder | What is in it |
|---|---|
| `catalogue/` | `studies.tsv`: Kremling2018 for expression, Lei2025 for the eQTL, TWAS and QTL results |
| `metadata/` | `runs.tsv`, `samples.tsv`, `lines.tsv` (generated), `contrasts.tsv`, `needs_review.tsv`, `library_qc.tsv`, `library_exclusions.tsv`; builders and curated inputs (`tissues.tsv`, `line_aliases.tsv`, `label_fixes.tsv`) in `metadata/build/` |
| `pipeline/` | Download and quantification on NC State's Hazel cluster (Slurm, nf-core/rnaseq, salmon without length correction for 3′ tags); see `pipeline/README.md` |
| `results/` | Gene count tables by sequencing batch (`results/Kremling2018/<batch>/`), `de.parquet` and `de_summary.tsv`, expression per tissue, tissue specificity, variation among lines, `b73_unstable.tsv`, co-expression partners (`partners/`, `coexpr_strongest_pairs.tsv`) |
| `annotation/` | Gene annotation (as in MaizExpressN and C): sources in `annotation/source/`, gene, alias, GO, MaizeMine, literature and atlas tables; Lei et al. 2025's tables on v5 genes (`lei_*`) |
| `analysis/` | `common.R` (libraries in use, CPM), `de_tissue.R`, `expression_tissue.R`, `tissue_specificity.R`, `line_variation.R`, `coexpr_partners.R` |
| `app/` | The Shiny app (`app.R`, `R/`, `www/`, `methods.md`) and `prepare_data.R` |
| `tests/` | testthat suite: `Rscript tests/testthat.R` |

## Rebuild everything after a change

From the project root, in this order (each step reads the outputs of the one before):

```bash
Rscript metadata/build/extract_lei_tables.R          # Lei Tables S1-S12 -> metadata/build/source/
python3 -I metadata/build/build_metadata.py --write
python3 -I metadata/build/build_contrasts.py
python3 -I metadata/build/validate_contrasts.py
# quantification and pipeline/library_qc.py on Hazel (pipeline/README.md); gene counts to results/Kremling2018/
Rscript metadata/build/library_exclusions.R
Rscript analysis/de_tissue.R <contrast_id>           # one run per contrast (the 180-243-pair fits take hours)
Rscript analysis/de_tissue.R --merge
Rscript analysis/expression_tissue.R
Rscript analysis/tissue_specificity.R
Rscript analysis/line_variation.R
Rscript annotation/build_parquet.R                   # after the annotation fetchers, if the sources change
Rscript annotation/build_lei.R
Rscript analysis/coexpr_partners.R                   # ~16 min
Rscript app/prepare_data.R
Rscript tests/testthat.R
```

## Conventions

- Gene key: B73 RefGen_v5 gene IDs (`Zm00001eb…`); v4, v3, symbols and Arabidopsis names are search aliases.
  Lei et al.'s v3 genes are matched to v5 through MaizeGDB's cross-reference, every match kept and labelled.
- Expression is CPM from salmon's gene counts; tissue contrasts are log2(first tissue / second tissue).
- Panel contrasts pair each line with itself (day and night also within sampling group); B73 contrasts use B73's
  replicate libraries.
- Library QC is strict (`metadata/build/library_exclusions.R`).

## License

The code is released under the [MIT License](LICENSE) (© 2026 Rubén Rellán-Álvarez and the GEMMA Lab).
The expression data come from published studies whose raw reads are in public repositories; please cite them
(listed with their DOIs in `catalogue/studies.tsv`) when using the results.
