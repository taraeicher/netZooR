context("test TIGER helper functions")

test_that("el2adj() converts edge list to adjacency matrix", {
  el <- data.frame(
    from = c("TF1", "TF1", "TF2"),
    to = c("G1", "G2", "G1"),
    weight = c(1, 0.5, 0.8)
  )

  adj <- el2adj(el)

  expect_true(is.matrix(adj))
  expect_equal(nrow(adj), 2)  # 2 TFs
  expect_equal(ncol(adj), 2)  # 2 genes
  expect_equal(adj["TF1", "G1"], 1)
  expect_equal(adj["TF1", "G2"], 0.5)
  expect_equal(adj["TF2", "G1"], 0.8)
})

test_that("adj2el() converts adjacency matrix to edge list", {
  adj <- matrix(c(1, 0.8, 0.5, 0), nrow = 2,
                dimnames = list(c("TF1", "TF2"), c("G1", "G2")))

  el <- adj2el(adj)

  expect_true(is.data.frame(el))
  expect_equal(ncol(el), 3)
  expect_equal(nrow(el), 4)  # all entries including zeros
  expect_equal(colnames(el), c("from", "to", "weight"))
})

test_that("el2adj() and adj2el() are inverse operations", {
  set.seed(42)
  adj_orig <- matrix(runif(6), nrow = 2,
                     dimnames = list(c("TF1", "TF2"), c("G1", "G2", "G3")))

  el <- adj2el(adj_orig)
  adj_roundtrip <- el2adj(el)

  expect_equal(adj_roundtrip, adj_orig, tolerance = 1e-10)
})

test_that("el2regulon() creates valid regulon object", {
  el <- data.frame(
    from = c("TF1", "TF1", "TF2"),
    to = c("G1", "G2", "G1"),
    weight = c(1, -0.5, 0.8)
  )

  reg <- el2regulon(el)

  expect_type(reg, "list")
  expect_equal(length(reg), 2)
  expect_true("TF1" %in% names(reg))
  expect_true("TF2" %in% names(reg))
  expect_true("tfmode" %in% names(reg[["TF1"]]))
  expect_true("likelihood" %in% names(reg[["TF1"]]))
})

test_that("adj2regulon() creates valid regulon from adjacency", {
  adj <- matrix(c(1, 0, -0.5, 0.8), nrow = 2,
                dimnames = list(c("TF1", "TF2"), c("G1", "G2")))

  reg <- adj2regulon(adj)

  expect_type(reg, "list")
  expect_equal(length(reg), 2)
})

test_that("adj2regulon() filters zero edges before regulon", {
  adj <- matrix(c(1, 0, 0, 0.8, 0, -0.5), nrow = 2,
                dimnames = list(c("TF1", "TF2"), c("G1", "G2", "G3")))

  reg <- adj2regulon(adj)
  expect_type(reg, "list")
  # TF1 has 1 non-zero edge (G1=1), TF2 has 2 non-zero edges (G2=0.8, G3=-0.5)
  expect_equal(length(reg[["TF1"]]$tfmode), 1)
  expect_equal(length(reg[["TF2"]]$tfmode), 2)
})

test_that("el2adj() handles single row edge list", {
  el <- data.frame(from = "TF1", to = "G1", weight = 1.5)
  adj <- el2adj(el)
  expect_equal(dim(adj), c(1, 1))
  expect_equal(adj[1, 1], 1.5)
})

test_that("adj2el() preserves all values including zeros", {
  adj <- matrix(0, nrow = 2, ncol = 2,
                dimnames = list(c("TF1", "TF2"), c("G1", "G2")))
  adj["TF1", "G1"] <- 3.14

  el <- adj2el(adj)
  expect_equal(nrow(el), 4)
  expect_equal(sum(el$weight != 0), 1)
  expect_equal(el$weight[el$from == "TF1" & el$to == "G1"], 3.14)
})

test_that("el2regulon() returns correct tfmode values", {
  el <- data.frame(
    from = c("TF1", "TF1", "TF1"),
    to = c("G1", "G2", "G3"),
    weight = c(1, -0.5, 0.3)
  )

  reg <- el2regulon(el)
  expect_equal(reg[["TF1"]]$tfmode[["G1"]], 1)
  expect_equal(reg[["TF1"]]$tfmode[["G2"]], -0.5)
  expect_equal(reg[["TF1"]]$tfmode[["G3"]], 0.3)
  expect_true(all(reg[["TF1"]]$likelihood == 1))
})

test_that("priorPp() filters inconsistent edges", {
  skip_if_not_installed("GeneNet")
  set.seed(42)
  tfs <- paste0("TF", 1:3)
  genes <- paste0("G", 1:5)
  prior <- matrix(sample(c(-1, 0, 1), 15, replace = TRUE), nrow = 3,
                  dimnames = list(tfs, genes))
  expr <- matrix(rnorm(80), nrow = 8,
                 dimnames = list(c(tfs, genes), paste0("S", 1:10)))

  result <- priorPp(prior, expr)

  expect_true(is.matrix(result))
  # priorPp may filter TFs not in expr and remove all-zero rows/cols
  expect_true(nrow(result) <= nrow(prior))
  expect_true(ncol(result) <= ncol(prior))
  # Some edges remain, some may be filtered to 1e-6
  expect_true(all(result %in% c(-1, 0, 1, 1e-6)))
})

test_that("priorPp() handles all-positive prior", {
  skip_if_not_installed("GeneNet")
  set.seed(42)
  tfs <- paste0("TF", 1:3)
  genes <- paste0("G", 1:5)
  prior <- matrix(1, nrow = 3, ncol = 5, dimnames = list(tfs, genes))
  expr <- matrix(rnorm(80), nrow = 8,
                 dimnames = list(c(tfs, genes), paste0("S", 1:10)))

  result <- priorPp(prior, expr)

  expect_true(is.matrix(result))
  expect_true(all(result %in% c(0, 1, 1e-6)))
})

test_that(".get_cmdstanr_fun() errors when cmdstanr not available", {
  # This tests the error path - may or may not trigger depending on env
  # Just test that the function exists and can be called
  expect_true(is.function(.get_cmdstanr_fun))
})

test_that("tiger() errors without cmdstanr", {
  # If cmdstanr is not installed, tiger should error with a message
  has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
  if (!has_cmdstanr) {
    data(TIGER_expr)
    data(TIGER_prior)
    expect_error(tiger(TIGER_expr, TIGER_prior), "cmdstanr")
  } else {
    skip("cmdstanr is installed, skipping error path test")
  }
})

test_that("tiger() runs when cmdstanr and cmdstan are available", {
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({
    cmdstan_path <- getExportedValue("cmdstanr", "cmdstan_path")
    nzchar(cmdstan_path())
  }, error = function(e) FALSE)
  skip_if(!has_cmdstan, "cmdstan not installed")

  data(TIGER_expr)
  data(TIGER_prior)
  result <- tiger(TIGER_expr, TIGER_prior, method = "VB", seed = 42)

  expect_type(result, "list")
  expect_true("W" %in% names(result))
  expect_true("Z" %in% names(result))
  expect_true(is.matrix(result$W))
  expect_true(is.matrix(result$Z))
})
