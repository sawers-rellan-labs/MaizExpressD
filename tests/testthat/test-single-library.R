# Fold changes carried by one library (analysis/de_tissue.R, single_library; MaizExpressN rule, as in C): flagged,
# no p-value, log2FC kept for display.

# the rule itself, taken from analysis/de_tissue.R without running the script
rule <- local({
  env <- new.env()
  for (e in parse(file.path(ROOT, "analysis", "de_tissue.R")))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% c("SPIKE_SHARE", "SPIKE_PADJ", "is_single_library")) eval(e, env)
  env$is_single_library
})

test_that("the single-library rule flags spikes and leaves consistent responses alone", {
  m <- function(...) matrix(c(...), nrow = 1)
  expect_true(rule(m(0, 0, 96), m(0, 0, 0), 22, 1e-10))          # all reads in one library
  expect_true(rule(m(0, 0, 96), m(0, 0, 0), 22, NA))             # ... whatever the p-value
  expect_true(rule(m(1, 2, 300), m(3, 2, 1), 5, 0.2))            # one library > 90% of the counts
  expect_false(rule(m(30, 32, 35), m(0, 0, 0), 20, 1e-8))        # on/off, but in every library
  expect_true(rule(m(100, 5), m(6, 6), 3, 0.01))                 # significant, gone without the top library
  expect_false(rule(m(100, 5), m(6, 6), 3, 0.2))                 # ... the leave-one-out part needs padj < 0.05
  expect_true(rule(m(6, 6), m(100, 5), -3, 0.01))                # same, higher in the denominator
  expect_false(rule(m(40, 50, 60), m(10, 12, 9), 2.2, 1e-5))     # ordinary response
  expect_false(rule(m(0, 0), m(0, 0), NA, NA))                   # nothing to judge
})

test_that("single-library genes in results/de.parquet have no p-value and keep their log2FC", {
  de <- pq("results", "de.parquet")
  expect_false(anyNA(de$single_library))
  x <- de[single_library == TRUE]
  expect_gt(nrow(x), 0)
  expect_true(all(is.na(x$pvalue) & is.na(x$padj)))
  expect_true(all(is.finite(x$log2fc)))
  share <- de[, mean(single_library), by = contrast_id]$V1
  expect_true(all(share < 0.01))                                 # a small minority in every contrast
})
