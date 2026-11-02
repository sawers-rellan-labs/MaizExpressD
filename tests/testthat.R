# Run the MaizExpressD test suite from the project root:  Rscript tests/testthat.R
library(testthat)
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
