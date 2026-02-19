context("test SPIDER")

test_that("spider() runs with epifilter and returns panda object", {
  skip_on_cran()

  # Create minimal test data
  set.seed(42)
  tfs <- paste0("TF", 1:3)
  genes <- paste0("G", 1:5)

  motif <- data.frame(
    tf = rep(tfs, each = 5),
    gene = rep(genes, 3),
    score = rbinom(15, 1, 0.5)
  )

  epifilter <- motif
  nind <- sample(seq_len(nrow(epifilter)), 5)
  epifilter[nind, 3] <- 0

  expr <- matrix(rnorm(50), nrow = 5,
                 dimnames = list(genes, paste0("S", 1:10)))

  ppi <- data.frame(
    tf1 = c("TF1", "TF2"),
    tf2 = c("TF2", "TF3"),
    score = c(1, 1)
  )

  result <- spider(motif, expr, epifilter, ppi,
                   alpha = 0.1, hamming = 0.5, progress = FALSE)

  expect_s4_class(result, "panda")
  expect_true(!is.null(result@regNet))
  expect_equal(nrow(result@regNet), length(tfs))
  expect_equal(ncol(result@regNet), length(genes))
})

test_that("spider() works without epifilter", {
  skip_on_cran()

  set.seed(42)
  tfs <- paste0("TF", 1:3)
  genes <- paste0("G", 1:5)

  motif <- data.frame(
    tf = rep(tfs, each = 5),
    gene = rep(genes, 3),
    score = rbinom(15, 1, 0.5)
  )

  expr <- matrix(rnorm(50), nrow = 5,
                 dimnames = list(genes, paste0("S", 1:10)))

  result <- spider(motif, expr, epifilter = NULL, ppi = NULL,
                   alpha = 0.1, hamming = 0.5, progress = FALSE)

  expect_s4_class(result, "panda")
})
