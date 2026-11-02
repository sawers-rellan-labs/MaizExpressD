# Shared setup for the MaizExpressD tests (run with Rscript tests/testthat.R from the project root).
suppressPackageStartupMessages({ library(data.table); library(arrow) })
ROOT <- normalizePath(file.path(testthat::test_path(), "..", ".."))
pq <- function(...) as.data.table(read_parquet(file.path(ROOT, ...)))
tsv <- function(...) fread(file.path(ROOT, ...), colClasses = "character", na.strings = NULL)
split_ids <- function(x) setdiff(strsplit(x, ";", fixed = TRUE)[[1]], "")
GENE <- c(pep1 = "Zm00001eb383680", pdk1 = "Zm00001eb287770", ssu1 = "Zm00001eb197410", pcna2 = "Zm00001eb192050",
          cca1 = "Zm00001eb172450", lhy1 = "Zm00001eb415760", o2 = "Zm00001eb301570", zp27 = "Zm00001eb313800",
          sh2 = "Zm00001eb159060", nii1 = "Zm00001eb255880", nrt2 = "Zm00001eb209670")
PILOT_FC <- "HHV55BGXX"            # pilot leaf-3 libraries: L3Base from this flowcell is not used (analysis/common.R)
LEI_TISSUES <- c("GRoot", "GShoot", "L3Base", "L3Tip", "LMAD", "LMAN", "Kern")

skip_if_no_app <- function() testthat::skip_if_not(file.exists(file.path(ROOT, "app", "data", "BUILT")), "app data not built")

# the app's data helpers (app/R/data.R, plots.R, coexpr.R) on the built app data, loaded once
.app <- new.env()
app_env <- function() {
  skip_if_no_app()
  if (is.null(.app$e)) {
    e <- new.env()
    old <- Sys.getenv("MAIZEXPRESSD_DATA", NA)
    Sys.setenv(MAIZEXPRESSD_DATA = file.path(ROOT, "app", "data"))
    on.exit(if (is.na(old)) Sys.unsetenv("MAIZEXPRESSD_DATA") else Sys.setenv(MAIZEXPRESSD_DATA = old))
    for (f in c("data.R", "plots.R", "coexpr.R")) sys.source(file.path(ROOT, "app", "R", f), envir = e)
    .app$e <- e
  }
  .app$e
}
