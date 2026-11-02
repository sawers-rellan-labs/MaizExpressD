# Gene search (annotation/gene_aliases.parquet, as in N and C): any identifier or symbol, case-insensitive, every
# match returned.

al <- pq("annotation", "gene_aliases.parquet")
hits <- function(q) unique(al[alias_lower == tolower(q), gene_id])

test_that("v5, v4 and v3 IDs and the symbol all find pep1", {
  for (q in c("Zm00001eb383680", "Zm00001d046170", "GRMZM2G083841", "pep1")) expect_true(GENE[["pep1"]] %in% hits(q), info = q)
})

test_that("search is case-insensitive", {
  expect_identical(hits("PEP1"), hits("pep1"))
  expect_identical(hits("zm00001eb383680"), hits("Zm00001eb383680"))
})

test_that("Arabidopsis names find their maize orthologs, including dotted forms of colon names", {
  expect_true(GENE[["nii1"]] %in% hits("NIR1"))
  expect_true(GENE[["nii1"]] %in% hits("AT2G15620"))
  expect_true(GENE[["nrt2"]] %in% hits("NRT2.1"))
})

test_that("an identifier that maps to several genes returns all of them", {
  expect_gt(length(hits("psbA")), 5)                                # the plastid gene and its nuclear copies
  expect_gt(length(hits("NIA1")), 1)
})

test_that("every quantified gene is in the annotation", {
  tables <- Sys.glob(file.path(ROOT, "results", "Kremling2018", "*", "salmon.merged.gene_counts.tsv"))
  skip_if(length(tables) == 0, "no gene tables")
  g <- pq("annotation", "genes.parquet")
  counts <- fread(tables[1], select = "gene_id")
  expect_true(all(counts$gene_id %in% g$gene_id))
})
