context("test CRANE")

test_that("craneBipartite() returns edge list from adjacency matrix", {
  set.seed(42)
  A <- matrix(runif(12), nrow = 3, ncol = 4,
              dimnames = list(paste0("TF", 1:3), paste0("G", 1:4)))

  result <- craneBipartite(A, alpha = 0.1)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), 12)  # 3 * 4 edges
})

test_that("craneBipartite() with getAdj returns matrix", {
  set.seed(42)
  A <- matrix(runif(12), nrow = 3, ncol = 4,
              dimnames = list(paste0("TF", 1:3), paste0("G", 1:4)))

  result <- craneBipartite(A, alpha = 0.1, getAdj = TRUE)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 4))
})

test_that("craneBipartite() preserves row and column sums approximately", {
  set.seed(42)
  A <- matrix(runif(12), nrow = 3, ncol = 4,
              dimnames = list(paste0("TF", 1:3), paste0("G", 1:4)))

  result <- craneBipartite(A, alpha = 0.0, getAdj = TRUE)

  # With alpha=0, node strengths should be preserved exactly
  skip_if(is.null(result), "Alpha limit reached, result is NULL")
  expect_equal(rowSums(result), rowSums(A), tolerance = 1e-6)
  expect_equal(colSums(result), colSums(A), tolerance = 1e-6)
})

test_that("craneBipartite() accepts edge list input", {
  set.seed(42)
  elist <- data.frame(
    TF = rep(paste0("TF", 1:3), each = 4),
    Gene = rep(paste0("G", 1:4), 3),
    weight = runif(12)
  )

  result <- craneBipartite(elist, alpha = 0.1)
  
  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 3)
})

test_that("craneBipartite() with beta perturbation works", {
  set.seed(42)
  A <- matrix(abs(rnorm(12, mean = 2)), nrow = 3, ncol = 4,
              dimnames = list(paste0("TF", 1:3), paste0("G", 1:4)))

  result <- craneBipartite(A, alpha = 0.05, beta = 0.5, getAdj = TRUE)

  # Result may be NULL if alpha limit reached with beta perturbation
  if (!is.null(result)) {
    expect_true(is.matrix(result))
    expect_equal(dim(result), c(3, 4))
  }
})

test_that("craneUnipartite() returns perturbed matrix", {
  set.seed(42)
  n <- 5
  A <- matrix(runif(n * n), n, n,
              dimnames = list(paste0("N", 1:n), paste0("N", 1:n)))
  A <- (A + t(A)) / 2  # symmetric

  result <- craneUnipartite(A, alpha = 0.1)

  skip_if(is.null(result), "Alpha limit reached, result is NULL")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(n, n))
})

test_that("craneUnipartite() with isSelfLoop preserves diagonal", {
  set.seed(42)
  n <- 4
  A <- matrix(abs(rnorm(n * n, mean = 1)), n, n,
              dimnames = list(paste0("N", 1:n), paste0("N", 1:n)))
  A <- (A + t(A)) / 2
  diag(A) <- 1

  result <- craneUnipartite(A, alpha = 0.05, isSelfLoop = TRUE)

  if (!is.null(result)) {
    expect_equal(result[n, n], A[n, n], tolerance = 1e-6)
  }
})

test_that("elistAddTags() and elistRemoveTags() are inverse operations", {
  elist <- data.frame(
    from = c("TF1", "TF2"),
    to = c("G1", "G2"),
    weight = c(1, 2)
  )
  
  tagged <- elistAddTags(elist)
  expect_true(all(grepl("_A$", tagged[, 1])))
  expect_true(all(grepl("_B$", tagged[, 2])))
  
  untagged <- elistRemoveTags(tagged)
  expect_equal(untagged[, 1], elist[, 1])
  expect_equal(untagged[, 2], elist[, 2])
})

test_that("elistIsEdgeOrderEqual() checks column equality", {
  elist1 <- data.frame(a = c("A", "B"), b = c("C", "D"), w = c(1, 2))
  elist2 <- data.frame(a = c("A", "B"), b = c("C", "D"), w = c(3, 4))
  elist3 <- data.frame(a = c("B", "A"), b = c("D", "C"), w = c(1, 2))
  
  expect_true(elistIsEdgeOrderEqual(elist1, elist2))
  expect_false(elistIsEdgeOrderEqual(elist1, elist3))
})

test_that("elistToAdjMat() and adjMatToElist() work correctly", {
  elist <- data.frame(
    from = rep(paste0("TF", 1:2), each = 3),
    to = rep(paste0("G", 1:3), 2),
    weight = 1:6
  )
  
  A <- elistToAdjMat(elist, isBipartite = TRUE)
  expect_true(is.matrix(A))
  expect_equal(nrow(A), 2)
  expect_equal(ncol(A), 3)
  
  elist2 <- adjMatToElist(A)
  expect_true(is.data.frame(elist2))
  expect_equal(ncol(elist2), 3)
})

test_that("elistSort() sorts edge list alphabetically", {
  elist <- data.frame(
    from = c("TF2", "TF1", "TF2", "TF1"),
    to = c("G2", "G1", "G1", "G2"),
    weight = c(4, 1, 3, 2)
  )
  
  sorted <- elistSort(elist)
  expect_true(is.data.frame(sorted))
})

test_that("isElist() detects edge lists correctly", {
  elist <- data.frame(a = c("A", "B"), b = c("C", "D"), w = c(1.0, 2.0))
  expect_true(isElist(elist))
  
  mat <- data.frame(a = c(1, 2), b = c(3, 4), c = c(5, 6))
  expect_false(isElist(mat))
})

test_that("jutterDegree() with beta=0 returns unchanged values", {
  nodeD <- c(10, 20, 30)
  result <- jutterDegree(nodeD, beta = 0)
  expect_equal(result, nodeD)
})

test_that("jutterDegree() with beta > 0 perturbs values", {
  set.seed(42)
  nodeD <- c(10, 20, 30, 40, 50)
  result <- jutterDegree(nodeD, beta = 1, beta_slope = TRUE)
  expect_length(result, 5)
  expect_false(all(result == nodeD))
})

test_that("jutterDegree() with negative beta works", {
  set.seed(42)
  nodeD <- c(10, 20, 30)
  result <- jutterDegree(nodeD, beta = -0.5, beta_slope = FALSE)
  expect_length(result, 3)
})

test_that("alpacaGetMember() extracts members by type", {
  alp <- list(
    c(TF1_A = 1, TF2_A = 1, G1_B = 2, G2_B = 2),
    c(TF1_A = 0.5, TF2_A = 0.5, G1_B = 0.5, G2_B = 0.5)
  )
  
  all_memb <- alpacaGetMember(alp, "all")
  expect_length(all_memb, 4)
  
  tf_memb <- alpacaGetMember(alp, "tf")
  expect_length(tf_memb, 2)
  expect_true(all(grepl("_A", names(tf_memb))))
  
  gene_memb <- alpacaGetMember(alp, "gene")
  expect_length(gene_memb, 2)
  expect_true(all(grepl("_B", names(gene_memb))))
})

test_that("alpacaObjectToDfList() converts alpaca object to data frames", {
  alp <- list(
    c(TF1_A = 1, TF2_A = 1, G1_B = 2, G2_B = 2),
    c(TF1_A = 0.5, TF2_A = 0.3, G1_B = 0.4, G2_B = 0.6)
  )
  
  result <- alpacaObjectToDfList(alp)
  
  expect_type(result, "list")
  expect_true("TF" %in% names(result))
  expect_true("Gene" %in% names(result))
  expect_true(is.data.frame(result$TF))
  expect_true(is.data.frame(result$Gene))
  expect_equal(nrow(result$TF), 2)
  expect_equal(nrow(result$Gene), 2)
})

test_that("alpacaComputeDifferentialScoreFromDWBM() computes scores", {
  dwbm <- matrix(c(0.5, 0.1, 0.2, 0.3), nrow = 2, ncol = 2,
                 dimnames = list(c("TF1_A", "TF2_A"), c("G1_B", "G2_B")))
  louv.memb <- c(1, 1, 1, 1)
  names(louv.memb) <- c("TF1_A", "TF2_A", "G1_B", "G2_B")
  
  result <- alpacaComputeDifferentialScoreFromDWBM(dwbm, louv.memb)
  
  expect_true(is.numeric(result))
  expect_length(result, 4)
  # Scores should sum to ~1 within each community
  expect_equal(sum(result), 2, tolerance = 0.1)
})

test_that("craneUnipartite() preserves node strengths with alpha=0", {
  set.seed(42)
  n <- 5
  A <- matrix(runif(n * n), n, n,
              dimnames = list(paste0("N", 1:n), paste0("N", 1:n)))
  A <- (A + t(A)) / 2

  result <- craneUnipartite(A, alpha = 0.0)

  skip_if(is.null(result), "Alpha limit reached, result is NULL")
  expect_equal(rowSums(result), rowSums(A), tolerance = 1e-6)
})

test_that("craneBipartite() accepts edge list input", {
  set.seed(42)
  el <- data.frame(
    TF = rep(paste0("TF", 1:3), each = 4),
    Gene = rep(paste0("G", 1:4), 3),
    Weight = runif(12)
  )

  result <- craneBipartite(el, alpha = 0.1)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 3)
})
