# DRAGON cross-language parity test.
#
# Loads shared inputs and gold-standard outputs that are byte-identical to the
# fixtures in netZooPy's tests/dragon_parity/. If both this test and the
# matching netZooPy test pass, the two implementations agree to 1e-5 (atol).
# Coverage: lambdas, shrunken covariance, precision matrix, partial correlation
# matrix. Kappa / p-values are not covered because they are unimplemented in
# netZooR.

context("DRAGON cross-language parity (vs netZooPy gold)")

fixture_dir <- "dragon_parity"

X1 <- as.matrix(read.csv(file.path(fixture_dir, "X1.csv"), header = FALSE))
X2 <- as.matrix(read.csv(file.path(fixture_dir, "X2.csv"), header = FALSE))

py_lambdas <- as.numeric(readLines(file.path(fixture_dir, "lambdas.txt")))
py_cov     <- as.matrix(read.csv(file.path(fixture_dir, "cov.csv"),  header = FALSE))
py_prec    <- as.matrix(read.csv(file.path(fixture_dir, "prec.csv"), header = FALSE))
py_ggm     <- as.matrix(read.csv(file.path(fixture_dir, "ggm.csv"),  header = FALSE))

res <- dragon(layer1 = X1, layer2 = X2, pval = FALSE, verbose = FALSE)

test_that("[DRAGON parity] lambdas match netZooPy gold", {
  expect_equal(as.numeric(res$lambdas), py_lambdas, tolerance = 1e-5)
})

test_that("[DRAGON parity] shrunken covariance matches netZooPy gold", {
  expect_equal(as.vector(res$cov), as.vector(py_cov), tolerance = 1e-5)
})

test_that("[DRAGON parity] precision matrix matches netZooPy gold", {
  expect_equal(as.vector(res$prec), as.vector(py_prec), tolerance = 1e-5)
})

test_that("[DRAGON parity] partial correlation matrix matches netZooPy gold", {
  expect_equal(as.vector(res$ggm), as.vector(py_ggm), tolerance = 1e-5)
})
