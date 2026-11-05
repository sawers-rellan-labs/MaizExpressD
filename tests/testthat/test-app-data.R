# App data (app/prepare_data.R): the files hold exactly the project results, only libraries in use reach the app,
# and nothing in it says why a library or line was left out.

md5 <- function(f) unname(tools::md5sum(f))

test_that("the app's result tables and partner buckets are the project's own files", {
  skip_if_no_app()
  pairs <- c(de.parquet = "results/de.parquet", expr_tissue.parquet = "results/expression_tissue.parquet",
             specificity.parquet = "results/tissue_specificity.parquet", line_variation.parquet = "results/line_variation.parquet")
  for (f in names(pairs)) expect_equal(md5(file.path(ROOT, "app", "data", f)), md5(file.path(ROOT, pairs[[f]])), label = f)
  for (f in Sys.glob(file.path(ROOT, "annotation", "lei_*.parquet")))
    expect_equal(md5(file.path(ROOT, "app", "data", basename(f))), md5(f), label = basename(f))
  a <- Sys.glob(file.path(ROOT, "app", "data", "partners", "*.parquet"))
  b <- Sys.glob(file.path(ROOT, "results", "partners", "*.parquet"))
  expect_equal(basename(a), basename(b))
  expect_equal(md5(a), md5(b))
})

test_that("the per-line CPM buckets hold every row of results/expression_line.parquet", {
  skip_if_no_app()
  f <- file.path(ROOT, "results", "expression_line.parquet")
  skip_if_not(file.exists(f), "results/expression_line.parquet is not in the repository (analysis/expression_tissue.R rebuilds it)")
  n <- sum(vapply(Sys.glob(file.path(ROOT, "app", "data", "expr", "*.parquet")),
                  function(b) nrow(read_parquet(b, col_select = "gene_id")), numeric(1)))
  expect_equal(n, nrow(read_parquet(f, col_select = "gene_id")))
})

test_that("a gene lookup is quick and returns every line of every tissue", {
  e <- app_env()
  t0 <- Sys.time(); x <- e$get_expr(GENE[["pep1"]]); p <- e$get_partners(GENE[["pep1"]])
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 1)
  expect_equal(nrow(x), nrow(e$LINE_TISSUE))                    # lines without reads appear with CPM 0
  expect_false(anyNA(x$cpm)); expect_false(anyNA(x$subpopulation))
  f <- file.path(ROOT, "results", "expression_line.parquet")
  if (file.exists(f)) {
    full <- as.data.table(read_parquet(f))[gene_id == GENE[["pep1"]]]
    expect_equal(x[cpm > 0, .N], nrow(full))
  }
  expect_gt(nrow(p), 0)
})

test_that("only libraries in use reach the app: QC passed, no B73-10, no pilot leaf-3 base", {
  skip_if_no_app()
  s <- pq("app", "data", "samples.parquet")
  expect_length(intersect(tsv("metadata", "library_exclusions.tsv")$sample_id, s$sample_id), 0)
  expect_false(any(s$genotype == "B73-10"))
  expect_false(any(s$tissue == "L3Base" & s$flowcell == PILOT_FC))
  expect_setequal(s$sample_id, unlist(lapply(tsv("results", "expression_lines.tsv")$sample_ids, split_ids)))
  ln <- pq("app", "data", "lines.parquet"); ts <- pq("app", "data", "tissues.parquet")
  expect_setequal(ln$genotype, unique(s$genotype))
  expect_equal(ts$tissue, c("GRoot", "GShoot", "L3Base", "L3Mid", "LMid", "L3Tip", "LMAD", "LMAN", "Kern"))
  expect_equal(ts$n_libraries, s[, .N, by = tissue][match(ts$tissue, tissue), N])
  expect_equal(ts[tissue == "LMid", tissue_label], "Leaf 3 middle (pilot)")
})

test_that("nothing in the app says why a library or line was left out", {
  skip_if_no_app()
  own <- c("samples", "lines", "line_tissue", "tissues", "contrasts")       # the app's own tables (not annotation)
  rx <- "exclu|fail|\\bQC\\b|contamin|swap|outlier|left out|mislabel|curation"
  for (t in own) {
    d <- pq("app", "data", paste0(t, ".parquet"))
    expect_false(any(grepl("reason|exclu|curation|_qc$", names(d))), label = t)
    for (cl in names(d)[vapply(d, is.character, TRUE)]) expect_false(any(grepl(rx, d[[cl]], ignore.case = TRUE)), label = paste(t, cl))
  }
  expect_false("excluded_libraries" %in% names(pq("app", "data", "contrasts.parquet")))
  for (f in c(file.path(ROOT, "app", "methods.md"), Sys.glob(file.path(ROOT, "app", "*.R")), Sys.glob(file.path(ROOT, "app", "R", "*.R")))) {
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl("private/", txt, fixed = TRUE), label = basename(f))
    expect_false(grepl(PILOT_FC, txt, fixed = TRUE), label = basename(f))
  }
})

test_that("p-values that underflowed show as < 1e-300 and sort first", {
  e <- app_env()
  expect_true(any(pq("results", "de.parquet")$padj == 0, na.rm = TRUE))   # the strongest panel results reach exactly 0
  expect_equal(e$fmtp(c(0, 4e-317, 1e-300, 2e-5, 0.0123, NA)), c("< 1e-300", "< 1e-300", "1.0e-300", "2.0e-05", "0.0123", "–"))
  expect_equal(e$p_tab(c(0, 4e-317, 3.456e-30, 0.0123, NA)), c(0, 0, 3.46e-30, 0.0123, NA))
})
