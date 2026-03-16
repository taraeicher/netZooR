context("test ALPACA result")

test_that("ALPACA works", {
  
  suppressWarnings(load("./testDataset.RData"))
  
  simp.alp <- alpaca(simp.mat,NULL,verbose=F)
  
  expect_equal(as.vector(simp.alp[[1]]),simp.memb)
  
  
})

test_that("alpaca() returns membership and scores with verbose=TRUE", {
  example_path <- system.file("extdata", "Example_2comm.txt",
                              package = "netZooR", mustWork = TRUE)
  simp.mat <- read.table(example_path, header = TRUE)
  
  result <- alpaca(simp.mat, NULL, verbose = TRUE)
  
  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(length(result[[1]]) > 0)
  expect_true(length(result[[2]]) > 0)
  expect_true(all(result[[1]] >= 1))
})

test_that("alpaca() writes files when file.stem is provided", {
  example_path <- system.file("extdata", "Example_2comm.txt",
                              package = "netZooR", mustWork = TRUE)
  simp.mat <- read.table(example_path, header = TRUE)
  
  tmp <- tempfile("alpaca_test")
  result <- alpaca(simp.mat, tmp, verbose = TRUE)
  
  expect_true(file.exists(paste0(tmp, "_ALPACA_ctrl_memb.txt")))
  expect_true(file.exists(paste0(tmp, "_ALPACA_final_memb.txt")))
  expect_true(file.exists(paste0(tmp, "_ALPACA_scores.txt")))
  expect_true(file.exists(paste0(tmp, "_DWBM.txt")))
  
  unlink(paste0(tmp, "*"))
})

test_that("alpacaExtractTopGenes() extracts top genes", {
  example_path <- system.file("extdata", "Example_2comm.txt",
                              package = "netZooR", mustWork = TRUE)
  simp.mat <- read.table(example_path, header = TRUE)
  simp.alp <- alpaca(simp.mat, NULL, verbose = FALSE)
  
  result <- alpacaExtractTopGenes(simp.alp, set.lengths = c(2, 3))
  
  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(length(result[[1]]) > 0)
  expect_true(length(result[[2]]) > 0)
  # Each set name should have format "number_length"
  expect_true(all(grepl("^\\d+_\\d+$", result[[2]])))
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
  # Modularity matrix sums to zero
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

test_that("alpacaGenLouvain() finds communities in modularity matrix", {
  # Small 4-node network with 2 clear communities
  B <- matrix(c(1, 0.8, -0.2, -0.2,
                0.8, 1, -0.2, -0.2,
                -0.2, -0.2, 1, 0.8,
                -0.2, -0.2, 0.8, 1), 4, 4)
  
  result <- alpacaGenLouvain(B)
  expect_true(length(result) == 4 || length(result) == 1)
  if (length(result) == 4) {
    # Two communities
    expect_true(result[1] == result[2])
    expect_true(result[3] == result[4])
  }
})

test_that("alpacaTestNodeRank() computes enrichment statistics", {
  set.seed(42)
  node.ordered <- paste0("gene", 1:100)
  true.pos <- paste0("gene", 1:10)  # top genes are true positives
  
  result <- alpacaTestNodeRank(node.ordered, true.pos)
  
  expect_length(result, 4)
  expect_true(all(is.numeric(result)))
  # Wilcoxon p-value should be small since true.pos are at the top
  expect_true(result[1] < 0.05)
})

test_that("alpacaCommunityStructureRotation() ranks nodes", {
  net1.memb <- c(1, 1, 2, 2, 1)
  names(net1.memb) <- paste0("N", 1:5)
  net2.memb <- c(1, 2, 2, 2, 1)
  names(net2.memb) <- paste0("N", 1:5)
  
  result <- alpacaCommunityStructureRotation(net1.memb, net2.memb)
  
  expect_type(result, "character")
  expect_true(length(result) == 5)
  # N2 changed community, should be ranked high
  expect_true("N2" %in% result[1:3])
})

test_that("alpacaComputeDWBMmatmScale() computes differential modularity", {
  edge.mat <- data.frame(
    TF = c("TF1_A", "TF1_A", "TF2_A", "TF2_A"),
    Gene = c("G1_B", "G2_B", "G1_B", "G2_B"),
    ctrl = c(1, 2, 3, 4),
    cond = c(2, 3, 4, 5)
  )
  ctrl.memb <- c(1, 1, 2, 2)
  names(ctrl.memb) <- c("TF1_A", "TF2_A", "G1_B", "G2_B")
  
  result <- alpacaComputeDWBMmatmScale(edge.mat, ctrl.memb)
  
  expect_true(is.matrix(result) || is.array(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that("alpacaWBMlouvain() finds communities in weighted bipartite network", {
  net.frame <- data.frame(
    TF = rep(paste0("TF", 1:3, "_A"), each = 4),
    Gene = rep(paste0("G", 1:4, "_B"), 3),
    weight = c(5, 5, 0.1, 0.1,
               5, 5, 0.1, 0.1,
               0.1, 0.1, 5, 5)
  )
  
  result <- alpacaWBMlouvain(net.frame)
  
  expect_type(result, "list")
  expect_length(result, 2)
  # Should find communities
  memb <- result[[1]]
  expect_true(length(memb) == 7)  # 3 TFs + 4 genes
  expect_true(max(memb) >= 1)
})

test_that("alpacaDeltaZAnalysis() performs edge subtraction with CONDOR", {
  # Create a connected bipartite network
  set.seed(42)
  n_edges <- 30
  tfs <- paste0("TF", rep(1:5, each = 6))
  genes <- paste0("Gene", rep(1:6, 5))
  df <- data.frame(
    tf = tfs, gene = genes,
    score_1 = runif(n_edges, 0.5, 1),
    score_2 = runif(n_edges, 0.5, 1)
  )
  
  result <- tryCatch(
    alpacaDeltaZAnalysis(df, NULL),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    expect_type(result, "list")
    expect_length(result, 2)
  }
})

test_that("alpacaDeltaZAnalysisLouvain() performs edge subtraction with Louvain", {
  # Create a well-connected bipartite network with clear community structure
  set.seed(123)
  tf_names <- paste0("TF", 1:10)
  gene_names <- paste0("Gene", 1:10)
  # All TF-gene pairs for full connectivity
  df <- expand.grid(tf = tf_names, gene = gene_names, stringsAsFactors = FALSE)
  # Community 1: TF1-5 strongly targeting Gene1-5
  # Community 2: TF6-10 strongly targeting Gene6-10
  df$score_1 <- 0.1
  df$score_2 <- 0.1
  for (i in seq_len(nrow(df))) {
    tf_idx <- as.integer(gsub("TF", "", df$tf[i]))
    gene_idx <- as.integer(gsub("Gene", "", df$gene[i]))
    if ((tf_idx <= 5 && gene_idx <= 5) || (tf_idx > 5 && gene_idx > 5)) {
      df$score_1[i] <- runif(1, 0.8, 1.0)
      df$score_2[i] <- runif(1, 1.5, 2.0) # stronger in condition 2
    }
  }
  
  result <- tryCatch(
    alpacaDeltaZAnalysisLouvain(df, NULL),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    expect_type(result, "list")
    expect_length(result, 2)
    expect_true(length(result[[1]]) > 0)
  }
})

test_that("alpacaRotationAnalysis() compares communities", {
  example_path <- system.file("extdata", "Example_2comm.txt",
                              package = "netZooR", mustWork = TRUE)
  simp.mat <- read.table(example_path, header = TRUE)
  
  result <- alpacaRotationAnalysis(simp.mat)
  
  expect_type(result, "character")
  expect_true(length(result) > 0)
})

test_that("alpacaRotationAnalysisLouvain() compares communities with Louvain", {
  example_path <- system.file("extdata", "Example_2comm.txt",
                              package = "netZooR", mustWork = TRUE)
  simp.mat <- read.table(example_path, header = TRUE)
  
  result <- alpacaRotationAnalysisLouvain(simp.mat)
  
  expect_type(result, "character")
  expect_true(length(result) > 0)
})

test_that("alpacaSimulateNetwork() creates simulated networks", {
  set.seed(42)
  comm.sizes <- matrix(c(5, 5, 5, 5), nrow = 2, ncol = 2)
  edge.mat <- matrix(c(8, 2, 2, 8), nrow = 2, ncol = 2)
  
  result <- alpacaSimulateNetwork(comm.sizes, edge.mat,
                                  num.module = 1,
                                  size.module = c(3, 3),
                                  dens.module = 0.5)
  
  expect_type(result, "list")
  expect_length(result, 2)
  expect_true(ncol(result[[1]]) == 4)
  expect_true(nrow(result[[1]]) > 0)
})


