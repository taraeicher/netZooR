context("test ALPACA result")

test_that("ALPACA works", {
  
  suppressWarnings(load("./testDataset.RData"))
  
  simp.alp <- alpaca(simp.mat,NULL,verbose=F)
  
  expect_equal(as.vector(simp.alp[[1]]),simp.memb)
  
  
})

test_that("alpacaTopEnsembltoTopSym() converts gene identifiers", {
  mod.top <- list(c("ENSG001", "ENSG002", "ENSG003"),
                  c("ENSG002", "ENSG004"))
  annot.vec <- c(ENSG001 = "TP53", ENSG002 = "BRCA1", ENSG003 = "MYC")
  
  result <- alpacaTopEnsembltoTopSym(mod.top, annot.vec)
  
  expect_type(result, "list")
  expect_length(result, 2)
  expect_equal(unname(result[[1]]), c("TP53", "BRCA1", "MYC"))
  expect_equal(unname(result[[2]]), "BRCA1")
})

test_that("alpacaNodeToGene() removes tags", {
  expect_equal(alpacaNodeToGene("TP53_A"), "TP53")
  expect_equal(alpacaNodeToGene("BRCA1_B"), "BRCA1")
  expect_equal(alpacaNodeToGene("MYC"), "MYC")
})

test_that("alpacaComputeWBMmat() computes modularity matrix", {
  edge.mat <- data.frame(
    TF = c("TF1_A", "TF1_A", "TF2_A", "TF2_A"),
    Gene = c("G1_B", "G2_B", "G1_B", "G2_B"),
    weight = c(1, 2, 3, 4)
  )
  
  result <- alpacaComputeWBMmat(edge.mat)
  
  expect_true(is.matrix(result) || is.array(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(sum(result), 0, tolerance = 1e-10)
})

test_that("alpacaMetaNetwork() computes meta network", {
  J <- matrix(c(1, 0.5, 0.5, 2), 2, 2)
  S <- c(1, 1)
  result <- alpacaMetaNetwork(J, S)
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 1)
  
  S2 <- c(1, 2)
  result2 <- alpacaMetaNetwork(J, S2)
  expect_equal(nrow(result2), 2)
  expect_equal(ncol(result2), 2)
})

test_that("alpacaTidyConfig() renumbers communities", {
  S <- c(3, 3, 5, 5, 3)
  result <- alpacaTidyConfig(S)
  expect_equal(result, c(1, 1, 2, 2, 1))
  
  S2 <- c(1, 2, 3)
  expect_equal(alpacaTidyConfig(S2), c(1, 2, 3))
})


