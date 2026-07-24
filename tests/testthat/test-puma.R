context("test PUMA")

test_that("puma() runs and returns panda object", {
  skip_on_cran()

  # Create minimal test data
  set.seed(42)
  tfs <- paste0("TF", 1:3)
  genes <- paste0("G", 1:5)
  mirs <- paste0("TF", 1:3)

  motif <- data.frame(
    tf = rep(tfs, each = 5),
    gene = rep(genes, 3),
    score = rbinom(15, 1, 0.5)
  )

  expr <- matrix(rnorm(50), nrow = 5,
                 dimnames = list(genes, paste0("S", 1:10)))

  result <- puma(motif, expr, ppi = NULL, mir_file = mirs,
                 alpha = 0.1, hamming = 0.5, progress = FALSE)

  expect_s4_class(result, "panda")
  expect_true(!is.null(result@regNet))
  expect_equal(nrow(result@regNet), length(tfs))
  expect_equal(ncol(result@regNet), length(genes))
})

test_that("puma() works without expression data", {
  skip_on_cran()

  set.seed(42)
  tfs <- paste0("TF", 1:3)
  genes <- paste0("G", 1:5)

  motif <- data.frame(
    tf = rep(tfs, each = 5),
    gene = rep(genes, 3),
    score = rbinom(15, 1, 0.5)
  )

  result <- expect_warning(
    puma(motif, expr = NULL, ppi = NULL, mir_file = tfs,
         alpha = 0.1, hamming = 0.5, progress = FALSE),
    "No expression data given")

  expect_s4_class(result, "panda")
})
