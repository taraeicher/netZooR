# unit-tests for UNAGI
context("test UNAGI functions")
test_that("[UNAGI] BreadthFirstSearchU() function yields expected results", {
  
  # Construct a starting network, which will be modified.
  # Here, we expect genes A and B to be connected after 1 hop via TF2, genes A and D
  # to be connected after 1 hop via TF3, and genes A and C to be connected after 2
  # hops via TF4.
  bindWithReversed <- function(dataFrame) {
    reversedDataFrame <- dataFrame
    reversedDataFrame$source <- dataFrame$target
    reversedDataFrame$target <- dataFrame$source
    return(rbind(dataFrame, reversedDataFrame))
  }
  startingNetwork <- bindWithReversed(data.frame(source = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4"),
                                                 target = c("geneB", "geneC", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD"),
                                                 score = c(3, 5, 4, 4, 5, 3, 1, 1)))
  rownames(startingNetwork) <- paste(startingNetwork$source, startingNetwork$target, sep = "__")
  
  # Test that errors are thrown when appropriate.
  expect_error(BreadthFirstSearchU(networks = startingNetwork,
                                   startingNodes = c("blob", "fish"),
                                   nodesToExclude = c()),
               "ERROR: Starting nodes do not overlap with network nodes")
  expect_error(BreadthFirstSearchU(networks = startingNetwork, 
                                   startingNodes = c("geneA", "geneB", "geneC"),
                                   nodesToExclude = c("blob", "fish")),
               "ERROR: List of nodes to exclude does not overlap with network nodes")
  expect_error(BreadthFirstSearchU(networks = startingNetwork, 
                                   startingNodes = c("geneA", "geneB", "geneC"),
                                   nodesToExclude = c("geneA", "geneB")),
               "ERROR: Starting nodes cannot overlap with nodes to exclude")
  
  # Ensure that, when starting from genes A, B, and C, we obtain the correct values.
  result1 <- BreadthFirstSearchU(networks = startingNetwork,
                                 startingNodes = c("geneA", "geneB", "geneC"),
                                 nodesToExclude = c())
  expect_setequal(rownames(result1), setdiff(rownames(startingNetwork), c("gene3__geneD", "gene4__geneD", "geneD__gene3", "geneD__gene4")))
  
  # Ensure that, when starting from genes gene2, gene3, and gene4 and removing genes A, B, and C, we obtain the correct values.
  result2 <- BreadthFirstSearchU(networks = startingNetwork,
                                 startingNodes = c("gene2", "gene3", "gene4"),
                                 nodesToExclude = c("geneA", "geneB", "geneC"))
  expect_setequal(rownames(result2), c("gene3__geneD", "gene4__geneD", "geneD__gene3", "geneD__gene4"))
  
  # Unipartite structure.
  startingNetworkU <- bindWithReversed(data.frame(source = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4", "gene1", "gene1", "geneA", "geneB"),
                                                  target = c("geneB", "geneC", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD", "gene2", "gene3", "geneD", "geneC"),
                                                  score = c(3, 5, 4, 4, 5, 3, 1, 1, 1, 1, 2, 2)))
  rownames(startingNetworkU) <- paste(startingNetworkU$source, startingNetworkU$target, sep = "__")
  
  # Ensure that, when starting from genes A, B, and C, we obtain the correct values.
  result1 <- BreadthFirstSearchU(networks = startingNetworkU,
                                 startingNodes = c("geneA", "geneB", "geneC"),
                                 nodesToExclude = c())
  expect_setequal(rownames(result1), setdiff(rownames(startingNetworkU), c("gene3__geneD", "gene4__geneD", "gene1__gene2", "gene1__gene3",
                                                                           "geneD__gene3", "geneD__gene4", "gene2__gene1", "gene3__gene1")))
  
  # Ensure that, when starting from genes gene2, gene3, and gene4 and removing genes A, B, and C, we obtain the correct values.
  result2 <- BreadthFirstSearchU(networks = startingNetworkU,
                                 startingNodes = c("gene2", "gene3", "gene4"),
                                 nodesToExclude = c("geneA", "geneB", "geneC"))
  expect_setequal(rownames(result2), c("gene3__geneD", "gene4__geneD", "gene1__gene2", "gene1__gene3",
                                       "geneD__gene3", "geneD__gene4", "gene2__gene1", "gene3__gene1"))
  
  # Test messaging.
  expect_message(BreadthFirstSearchU(networks = startingNetworkU,
                                     startingNodes = c("gene2", "gene3", "gene4"),
                                     nodesToExclude = c("geneA", "geneB", "geneC"), verbose = TRUE), "Retained 8 edges")
})
test_that("[UNAGI] FindEdgesForHopU() function yields expected results",{
  
  # Construct a starting network, which will be modified.
  # Here, we expect genes A and B to be connected after 1 hop via gene2, genes A and D
  # to be connected after 1 hop via TF3, and genes A and C to be connected after 2
  # hops via TF4.
  bindWithReversed <- function(dataFrame){
    reversedDataFrame <- dataFrame
    reversedDataFrame$source <- dataFrame$target
    reversedDataFrame$target <- dataFrame$source
    return(rbind(dataFrame, reversedDataFrame))
  }
  startingNetwork <- bindWithReversed(data.frame(source = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4"),
                                                 target = c("geneB", "geneC", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD"),
                                                 score = c(3, 5, 4, 4, 5, 3, 1, 1)))
  rownames(startingNetwork) <- paste(startingNetwork$source, startingNetwork$target, sep = "__")
  
  # Set up the subnetwork from the previous example.
  geneANetHop1 <- data.frame(source = c("gene2", "gene3", "geneA", "geneA"), target = c("geneA", "geneA", "gene2", "gene3"))
  rownames(geneANetHop1) <- paste(geneANetHop1$source, geneANetHop1$target, sep = "__")
  geneANetHop2 <- data.frame(source = c("gene2", "gene3", "geneB", "geneD"), target = c("geneB", "geneD", "gene2", "gene3"))
  rownames(geneANetHop2) <- paste(geneANetHop2$source, geneANetHop2$target, sep = "__")
  geneANetHop3 <- data.frame(source = c("gene1", "gene4", "geneB", "geneD"), target = c("geneB", "geneD", "gene1", "gene4"))
  rownames(geneANetHop3) <- paste(geneANetHop3$source, geneANetHop3$target, sep = "__")
  geneBNetHop1 <- data.frame(source = c("gene2", "gene1", "geneB", "geneB"), target = c("geneB", "geneB", "gene1", "gene2"))
  rownames(geneBNetHop1) <- paste(geneBNetHop1$source, geneBNetHop1$target, sep = "__")
  geneBNetHop2 <- data.frame(source = c("gene2", "gene1", "geneC", "geneA"), target = c("geneA", "geneC", "gene1", "gene2"))
  rownames(geneBNetHop2) <- paste(geneBNetHop2$source, geneBNetHop2$target, sep = "__")
  geneBNetHop3 <- data.frame(source = c("gene3", "gene4", "geneC", "geneA"), target = c("geneA", "geneC", "gene4", "gene3"))
  rownames(geneBNetHop3) <- paste(geneBNetHop3$source, geneBNetHop3$target, sep = "__")
  geneCNetHop1 <- data.frame(source = c("gene4", "gene1", "geneC", "geneC"), target = c("geneC", "geneC", "gene1", "gene4"))
  rownames(geneCNetHop1) <- paste(geneCNetHop1$source, geneCNetHop1$target, sep = "__")
  geneCNetHop2 <- data.frame(source = c("gene4", "gene1", "geneB", "geneD"), target = c("geneD", "geneB", "gene1", "gene4"))
  rownames(geneCNetHop2) <- paste(geneCNetHop2$source, geneCNetHop2$target, sep = "__")
  geneCNetHop3 <- data.frame(source = c("gene3", "gene2", "geneB", "geneD"), target = c("geneD", "geneB", "gene2", "gene3"))
  rownames(geneCNetHop3) <- paste(geneCNetHop3$source, geneCNetHop3$target, sep = "__")
  geneDNetHop1 <- data.frame(source = c("gene3", "gene4", "geneD", "geneD"), target = c("geneD", "geneD", "gene3", "gene4"))
  rownames(geneDNetHop1) <- paste(geneDNetHop1$source, geneDNetHop1$target, sep = "__")
  geneDNetHop2 <- data.frame(source = c("gene3", "gene4", "geneA", "geneC"), target = c("geneA", "geneC", "gene3", "gene4"))
  rownames(geneDNetHop2) <- paste(geneDNetHop2$source, geneDNetHop2$target, sep = "__")
  geneDNetHop3 <- data.frame(source = c("gene2", "gene1", "geneA", "geneC"), target = c("geneA", "geneC", "gene2", "gene1"))
  rownames(geneDNetHop3) <- paste(geneDNetHop3$source, geneDNetHop3$target, sep = "__")
  subnetworksFullU <- list(geneA = list(geneANetHop1, geneANetHop2, geneANetHop3),
                           geneB = list(geneBNetHop1, geneBNetHop2, geneBNetHop3),
                           geneC = list(geneCNetHop1, geneCNetHop2, geneCNetHop3),
                           geneD = list(geneDNetHop1, geneDNetHop2, geneDNetHop3))
  # Test method with all genes as input.
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneB", "geneC", "geneD"),
                               network = startingNetwork,
                               hopConstraint = 3)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[1]]), rownames(sigEdges$geneA[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[2]]), rownames(sigEdges$geneA[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[3]]), rownames(sigEdges$geneA[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[1]]), rownames(sigEdges$geneB[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[2]]), rownames(sigEdges$geneB[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[3]]), rownames(sigEdges$geneB[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[1]]), rownames(sigEdges$geneC[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[2]]), rownames(sigEdges$geneC[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[3]]), rownames(sigEdges$geneC[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneD[[1]]), rownames(sigEdges$geneD[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneD[[2]]), rownames(sigEdges$geneD[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneD[[3]]), rownames(sigEdges$geneD[[3]]))) == 0)
  
  # Test method with only genes A, B, C.
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneB", "geneC"),
                               network = startingNetwork, hopConstraint = 3)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[1]]), rownames(sigEdges$geneA[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[2]]), rownames(sigEdges$geneA[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[3]]), rownames(sigEdges$geneA[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[1]]), rownames(sigEdges$geneB[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[2]]), rownames(sigEdges$geneB[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[3]]), rownames(sigEdges$geneB[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[1]]), rownames(sigEdges$geneC[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[2]]), rownames(sigEdges$geneC[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[3]]), rownames(sigEdges$geneC[[3]]))) == 0)
  
  # Test method with only genes A and C.
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneC"),
                               network = startingNetwork,
                               hopConstraint = 3)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[1]]), rownames(sigEdges$geneA[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[2]]), rownames(sigEdges$geneA[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[3]]), rownames(sigEdges$geneA[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[1]]), rownames(sigEdges$geneC[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[2]]), rownames(sigEdges$geneC[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[3]]), rownames(sigEdges$geneC[[3]]))) == 0)
  
  # Test the same for the unipartite case.
  startingNetworkU <- bindWithReversed(data.frame(source = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4", "gene1", "gene1", "geneA", "geneB"),
                                                  target = c("geneB", "geneC", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD", "gene2", "gene3", "geneD", "geneC"),
                                                  score = c(3, 5, 4, 4, 5, 3, 1, 1, 1, 1, 2, 2)))
  rownames(startingNetworkU) <- paste(startingNetworkU$source, startingNetworkU$target, sep = "__")
  geneAnetHop1U <- data.frame(source = c("gene2", "gene3", "geneA", "geneA", "geneA", "geneD"), 
                              target = c("geneA", "geneA", "gene2", "gene3", "geneD", "geneA"))
  rownames(geneAnetHop1U) <- paste(geneAnetHop1U$source, geneAnetHop1U$target, sep = "__")
  geneAnetHop2U <- data.frame(source = c("gene2", "gene3", "geneB", "geneD", "gene2", "gene1", "gene3", "gene1", "geneD", "gene4"), 
                              target = c("geneB", "geneD", "gene2", "gene3", "gene1", "gene2", "gene1", "gene3", "gene4", "geneD"))
  rownames(geneAnetHop2U) <- paste(geneAnetHop2U$source, geneAnetHop2U$target, sep = "__")
  geneAnetHop3U <- data.frame(source = c("gene1", "geneB", "geneC", "gene1", "gene4", "geneC", "geneB", "geneC"), 
                              target = c("geneB", "gene1", "gene1", "geneC", "geneC", "gene4", "geneC", "geneB"))
  rownames(geneAnetHop3U) <- paste(geneAnetHop3U$source, geneAnetHop3U$target, sep = "__")
  geneBnetHop1U <- data.frame(source = c("gene2", "gene1", "geneB", "geneB", "geneB", "geneC"), 
                              target = c("geneB", "geneB", "gene1", "gene2", "geneC", "geneB"))
  rownames(geneBnetHop1U) <- paste(geneBnetHop1U$source, geneBnetHop1U$target, sep = "__")
  geneBnetHop2U <- data.frame(source = c("gene2", "gene1", "geneC", "geneA", "gene1", "gene1", "gene2", "gene3", "gene4", "geneC"), 
                              target = c("geneA", "geneC", "gene1", "gene2", "gene2", "gene3", "gene1", "gene1", "geneC", "gene4"))
  rownames(geneBnetHop2U) <- paste(geneBnetHop2U$source, geneBnetHop2U$target, sep = "__")
  geneBnetHop3U <- data.frame(source = c("gene3", "geneA", "gene3", "geneD", "gene4", "geneD", "geneA", "geneD"), 
                              target = c("geneA", "gene3", "geneD", "gene3", "geneD", "gene4", "geneD", "geneA"))
  rownames(geneBnetHop3U) <- paste(geneBnetHop3U$source, geneBnetHop3U$target, sep = "__")
  geneCnetHop1U <- data.frame(source = c("gene4", "gene1", "geneC", "geneC", "geneB", "geneC"), 
                              target = c("geneC", "geneC", "gene1", "gene4", "geneC", "geneB"))
  rownames(geneCnetHop1U) <- paste(geneCnetHop1U$source, geneCnetHop1U$target, sep = "__")
  geneCnetHop2U <- data.frame(source = c("gene4", "gene1", "geneB", "geneD", "gene1", "gene2", "gene1", "gene3", "gene2", "geneB"), 
                              target = c("geneD", "geneB", "gene1", "gene4", "gene2", "gene1", "gene3", "gene1", "geneB", "gene2"))
  rownames(geneCnetHop2U) <- paste(geneCnetHop2U$source, geneCnetHop2U$target, sep = "__")
  geneCnetHop3U <- data.frame(source = c("gene3", "geneD", "gene2", "geneA", "gene3", "geneA", "geneA", "geneD"), 
                              target = c("geneD", "gene3", "geneA", "gene2", "geneA", "gene3", "geneD", "geneA"))
  rownames(geneCnetHop3U) <- paste(geneCnetHop3U$source, geneCnetHop3U$target, sep = "__")
  geneDnetHop1U <- data.frame(source = c("gene3", "gene4", "geneD", "geneD", "geneA", "geneD"), 
                              target = c("geneD", "geneD", "gene3", "gene4", "geneD", "geneA"))
  rownames(geneDnetHop1U) <- paste(geneDnetHop1U$source, geneDnetHop1U$target, sep = "__")
  geneDnetHop2U <- data.frame(source = c("gene3", "gene4", "geneA", "geneC", "gene3", "gene1", "gene2", "geneA"), 
                              target = c("geneA", "geneC", "gene3", "gene4", "gene1", "gene3", "geneA", "gene2"))
  rownames(geneDnetHop2U) <- paste(geneDnetHop2U$source, geneDnetHop2U$target, sep = "__")
  geneDnetHop3U <- data.frame(source = c("gene1", "geneC", "gene1", "geneB", "gene1", "gene2", "gene2", "geneB", "geneB", "geneC"), 
                              target = c("geneC", "gene1", "geneB", "gene1", "gene2", "gene1", "geneB", "gene2", "geneC", "geneB"))
  rownames(geneDnetHop3U) <- paste(geneDnetHop3U$source, geneDnetHop3U$target, sep = "__")
  subnetworksFullU <- list(geneA = list(geneAnetHop1U, geneAnetHop2U, geneAnetHop3U),
                           geneB = list(geneBnetHop1U, geneBnetHop2U, geneBnetHop3U),
                           geneC = list(geneCnetHop1U, geneCnetHop2U, geneCnetHop3U),
                           geneD = list(geneDnetHop1U, geneDnetHop2U, geneDnetHop3U))
  
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneB", "geneC", "geneD"),
                               network = startingNetworkU,
                               hopConstraint = 3)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[1]]), rownames(sigEdges$geneA[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[2]]), rownames(sigEdges$geneA[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneA[[3]]), rownames(sigEdges$geneA[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[1]]), rownames(sigEdges$geneB[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[2]]), rownames(sigEdges$geneB[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneB[[3]]), rownames(sigEdges$geneB[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[1]]), rownames(sigEdges$geneC[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[2]]), rownames(sigEdges$geneC[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneC[[3]]), rownames(sigEdges$geneC[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneD[[1]]), rownames(sigEdges$geneD[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneD[[2]]), rownames(sigEdges$geneD[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFullU$geneD[[3]]), rownames(sigEdges$geneD[[3]]))) == 0)
  
  # Test messaging.
  expect_message(FindEdgesForHopU(geneSet = c("geneA", "geneB", "geneC", "geneD"),
                                  network = startingNetworkU,
                                  hopConstraint = 3, verbose = TRUE), "Evaluating hop 1 for gene geneA")
  
  # Test thresholding by p-value.
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneB", "geneC", "geneD"),
                               network = startingNetworkU,
                               hopConstraint = 3, topX = 0.1)
  expect_length(rownames(sigEdges$geneA[[1]]), 1)
})
test_that("[UNAGI] FindConnectionsForAllHopCountsU() function yields expected results", {
  
  # Set up the subnetwork from the previous example.
  bindWithReversed <- function(dataFrame) {
    reversedDataFrame <- dataFrame
    reversedDataFrame$source <- dataFrame$target
    reversedDataFrame$target <- dataFrame$source
    return(rbind(dataFrame, reversedDataFrame))
  }
  geneANetHop1 <- bindWithReversed(data.frame(source = c("geneA", "geneA"), 
                                              target = c("gene2", "gene3")))
  geneBNetHop1 <- bindWithReversed(data.frame(source = "geneB", target = "gene2"))
  geneCNetHop1 <- bindWithReversed(data.frame(source = "geneC", target = "gene4"))
  geneDNetHop1 <- bindWithReversed(data.frame(source = c("geneD", "geneD"), 
                                              target = c("gene3", "gene4")))
  geneANetHop2 <- bindWithReversed(data.frame(source = c("gene2", "gene3"), 
                                              target = c("geneB", "geneD")))
  geneBNetHop2 <- bindWithReversed(data.frame(source = "gene2", target = "geneA"))
  geneCNetHop2 <- bindWithReversed(data.frame(source = "gene4", target = "geneD"))
  geneDNetHop2 <- bindWithReversed(data.frame(source = c("gene3", "gene4"), 
                                              target = c("geneA", "geneC")))
  geneANetHop3 <- bindWithReversed(data.frame(source = "geneD", target = "gene4"))
  geneBNetHop3 <- bindWithReversed(data.frame(source = "geneA", target = "gene3"))
  geneCNetHop3 <- bindWithReversed(data.frame(source = "geneD", target = "gene3"))
  geneDNetHop3 <- bindWithReversed(data.frame(source = "geneA", target = "gene2"))
  
  subnetworks <- list(geneA = list(geneANetHop1, geneANetHop2, geneANetHop3),
                      geneB = list(geneBNetHop1, geneBNetHop2, geneBNetHop3),
                      geneC = list(geneCNetHop1, geneCNetHop2, geneCNetHop3),
                      geneD = list(geneDNetHop1, geneDNetHop2, geneDNetHop3))
  
  # Function to sort results alphanumerically by gene.
  sortAlphanumeric <- function(result){
    sortedName <- unlist(lapply(1:nrow(result), function(i){
      gene1 <- result[i,1]
      gene2 <- result[i,2]
      rowname <- paste0(gene1, "__", gene2)
      if(gene2 < gene1){
        rowname <- paste0(gene2, "__", gene1)
      }
      return(rowname)
    }))
    return(sortedName)
  }
  
  # Obtain the overlaps for 1, 2, and 3 hops.
  result <- FindConnectionsForAllHopCountsU(subnetworks)
  
  if (is.data.frame(result) || is.matrix(result)) {
    expect_setequal(sortAlphanumeric(result),
                    c("gene2__geneA", "gene2__geneB", "gene3__geneA", "gene3__geneD",
                      "gene4__geneC", "gene4__geneD"))
    
  } else {
    stop("Error: FindConnectionsForAllHopCountsU did not return a data frame or matrix.")
  }
  
  # Obtain the overlaps for only the first two hops, for only A and C.
  subnetworks <- list(geneA = list(geneANetHop1, geneANetHop2),
                      geneC = list(geneCNetHop1, geneCNetHop2))
  result <- FindConnectionsForAllHopCountsU(subnetworks)
  
  if (is.data.frame(result) || is.matrix(result)) {
    expect_equal(nrow(result), 0)
  } else {
    stop("Error: FindConnectionsForAllHopCountsU did not return a data frame or matrix.")
  }
  
  # Set up as a subnetwork from the previous example (unipartite case).
  startingNetworkU <- bindWithReversed(data.frame(source = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4", "gene1", "gene1", "geneA", "geneB"),
                                                  target = c("geneB", "geneC", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD", "gene2", "gene3", "geneD", "geneC"),
                                                  score = c(3, 5, 4, 4, 5, 3, 1, 1, 1, 1, 2, 2)))
  rownames(startingNetworkU) <- paste(startingNetworkU$source, startingNetworkU$target, sep = "__")
  
  # 2-hop connections
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneB", "geneC", "geneD"),
                               network = startingNetworkU,
                               hopConstraint = 1)
  result <- FindConnectionsForAllHopCountsU(sigEdges)
  expResult <- c("gene2__geneA", "gene2__geneB", "gene3__geneA", "gene3__geneD", "geneA__geneD",
                 "gene4__geneD", "gene1__geneB", "geneB__geneC", "gene1__geneC", "gene4__geneC")
  expect_equal(length(setdiff(sortAlphanumeric(result), expResult)), 0)
  expect_equal(length(setdiff(expResult, sortAlphanumeric(result))), 0)
  
  # 4-hop connections between just A and D.
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneD"),
                               network = startingNetworkU,
                               hopConstraint = 2)
  result <- FindConnectionsForAllHopCountsU(sigEdges)
  expResult <- c("gene3__geneA", "gene3__geneD", "gene1__gene2", "gene1__gene3", "gene2__geneA",
                 "geneA__geneD")
  expect_equal(length(setdiff(sortAlphanumeric(result), expResult)), 0)
  
  # 6-hop connections between just A and D.
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneD"),
                               network = startingNetworkU,
                               hopConstraint = 3)
  result <- FindConnectionsForAllHopCountsU(sigEdges)
  expResult <- c("gene3__geneA", "gene3__geneD", "gene1__gene2", "gene1__gene3", "gene2__geneA",
                 "gene1__geneB", "gene1__geneC", "geneB__geneC", "gene4__geneC", "gene2__geneB",
                 "gene4__geneD", "geneA__geneD")
  expect_equal(length(setdiff(sortAlphanumeric(result), expResult)), 0)
  expect_equal(length(setdiff(expResult, sortAlphanumeric(result))), 0)
  
  # 8-hop connections with no connection from C to anything.
  startingNetworkChain <- bindWithReversed(data.frame(source = c("gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4", "gene1", "gene1", "geneA", "geneC", "geneE"),
                                                      target = c("geneB", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD", "gene2", "gene3", "geneD", "geneE", "gene5"),
                                                      score = c(3, 5, 4, 4, 5, 3, 1, 1, 1, 1, 2, 2)))
  rownames(startingNetworkChain) <- paste(startingNetworkChain$source, startingNetworkChain$target, sep = "__")
  
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneD"),
                               network = startingNetworkChain,
                               hopConstraint = 4)
  result <- FindConnectionsForAllHopCountsU(sigEdges)
  expResult <- c("gene2__geneA", "gene1__gene2", "gene1__geneB", "gene1__gene3",
                 "gene3__geneD", "gene3__geneA", "gene2__geneB", "geneA__geneD")
  expect_equal(length(setdiff(sortAlphanumeric(result), expResult)), 0)
  expect_equal(length(setdiff(expResult, sortAlphanumeric(result))), 0)
  
  # Check when there are paths with exactly 3 hops.
  startingNetwork3 <- bindWithReversed(data.frame(source = c("gene1", "gene1", "gene2", "gene4", "gene4", "gene1", "gene1", "geneA"),
                                                  target = c("geneB", "geneC", "geneA", "geneC", "geneD", "gene2", "gene3", "geneD"),
                                                  score = c(3, 5, 4, 4, 5, 3, 1, 1)))
  rownames(startingNetwork3) <- paste(startingNetwork3$source, startingNetwork3$target, sep = "__")
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneC"),
                               network = startingNetwork3,
                               hopConstraint = 2)
  result <- FindConnectionsForAllHopCountsU(sigEdges)
  expResult <- c("gene2__geneA", "gene1__gene2", "gene1__geneC", "geneA__geneD", "gene4__geneD", "gene4__geneC")
  expect_equal(length(setdiff(sortAlphanumeric(result), expResult)), 0)
  expect_equal(length(setdiff(expResult, sortAlphanumeric(result))), 0)
  
  # Check when the maximum count to check is odd (case 1).
  result <- FindConnectionsForAllHopCountsU(sigEdges, odd = TRUE)
  expect_equal(length(setdiff(sortAlphanumeric(result), expResult)), 0)
  expect_equal(length(setdiff(expResult, sortAlphanumeric(result))), 0)
  
  # Check when the maximum count to check is odd (case 2).
  startingNetworkChain <- bindWithReversed(data.frame(source = c("gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4", "gene1", "gene1", "geneA", "geneC", "geneE"),
                                                      target = c("geneB", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD", "gene2", "gene3", "geneD", "geneE", "gene5"),
                                                      score = c(3, 5, 4, 4, 5, 3, 1, 1, 1, 1, 2, 2)))
  rownames(startingNetworkChain) <- paste(startingNetworkChain$source, startingNetworkChain$target, sep = "__")
  
  sigEdges <- FindEdgesForHopU(geneSet = c("geneA", "geneD"),
                               network = startingNetworkChain,
                               hopConstraint = 4)
  result <- FindConnectionsForAllHopCountsU(sigEdges, odd = TRUE)
  expResult <- c("gene2__geneA", "gene1__gene2", "gene1__geneB", "gene1__gene3",
                 "gene3__geneD", "gene3__geneA", "gene2__geneB", "geneA__geneD")
  expect_equal(length(setdiff(sortAlphanumeric(result), expResult)), 0)
  expect_equal(length(setdiff(expResult, sortAlphanumeric(result))), 0)
  
  # Test messaging.
  expect_message(FindConnectionsForAllHopCountsU(sigEdges, odd = TRUE, verbose = TRUE), 
                 "Hop 2 - 1 overlapped between geneA and geneD")
})
test_that("[UNAGI] RunUNAGI() function yields expected results",{
  
  # Create a dummy network.
  bindWithReversed <- function(dataFrame){
    reversedDataFrame <- data.frame(source = dataFrame$target, target = dataFrame$source)
    return(rbind(dataFrame, reversedDataFrame))
  }
  startingNetwork <- data.frame(source = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4"),
                                target = c("geneB", "geneC", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD"),
                                score = c(3, 5, 4, 4, 5, 3, 1, 1))
  rownames(startingNetwork) <- paste(startingNetwork$source, startingNetwork$target, sep = "__")
  
  
  # Check errors.
  expect_error(RunUNAGI(nodeSet = 4, network = startingNetwork, hopConstraint = 5),
               paste("Wrong input type! nodeSet must be a character vector. network must be a data frame.",
                     "hopConstraint must be a scalar numeric value."))
  expect_error(RunUNAGI(nodeSet = "g1", network = 67, hopConstraint = 5),
               paste("Wrong input type! nodeSet must be a character vector. network must be a data frame.",
                     "hopConstraint must be a scalar numeric value."))
  expect_error(RunUNAGI(nodeSet = "g1", network = list(), hopConstraint = 5),
               paste("Wrong input type! nodeSet must be a character vector. network must be a data frame.",
                     "hopConstraint must be a scalar numeric value."))
  expect_error(RunUNAGI(nodeSet = "g1", network = list(), hopConstraint = "hi"),
               paste("Wrong input type! nodeSet must be a character vector. network must be a data frame.",
                     "hopConstraint must be a scalar numeric value."))
  expect_error(RunUNAGI(nodeSet = "g1", network = 1, hopConstraint = 5),
               "network must be a data frame.")
  expect_error(RunUNAGI(nodeSet = "g1", network = data.frame(mushroom = NA, target = NA, score = NA, blah = NA), 
                        hopConstraint = 5),
               paste("Network must have source genes in the first column,",
                     "target genes in the second column, and scores in the third column."))
  expect_error(RunUNAGI(nodeSet = "g1", network = data.frame(source = NA, target = NA, score = NA), 
                        hopConstraint = -4), "hopConstraint must be at least 1.")
  expect_error(RunUNAGI(nodeSet = "gene1", network = startingNetwork, hopConstraint = 4),
               "Node set must contain at least 2 nodes.")
  expect_no_error(RunUNAGI(nodeSet = c("gene1", "gene2", "gene3"), network = startingNetwork, hopConstraint = 4))
  expect_no_error(RunUNAGI(nodeSet = c("gene1", "gene2", "gene3"), network = startingNetwork, hopConstraint = 7))
})
test_that("[UNAGI] PlotNetworkU() function runs without error",{
  
  # Create a dummy network.
  bindWithReversed <- function(dataFrame){
    reversedDataFrame <- data.frame(source = dataFrame$target, target = dataFrame$source)
    return(rbind(dataFrame, reversedDataFrame))
  }
  startingNetwork <- data.frame(source = c("gene1", "gene1", "gene2", "gene2", "gene3", "gene3", "gene4", "gene4"),
                                target = c("geneB", "geneC", "geneA", "geneB", "geneA", "geneD", "geneC", "geneD"),
                                score = c(3, 5, 4, 4, 5, 3, 1, 1))
  rownames(startingNetwork) <- paste(startingNetwork$source, startingNetwork$target, sep = "__")
  
  # Run UNAGI.
  unagi <- RunUNAGI(nodeSet = c("gene1", "gene2", "gene3"), network = startingNetwork, hopConstraint = 4)
  
  # Plot the network.
  expect_no_error(PlotNetworkU(network = unagi,
                               vertexLabels = c("gene1", "gene2", "gene3", "gene4", "geneA", "geneB", "geneC", "geneD"),
               geneColorMapping = data.frame(gene = c("gene1", "gene2", "gene3", "gene4", "geneA", "geneB", "geneC", "geneD"), 
                                             color = c("red", "green", "blue", "orange", "purple", "yellow", "black", "gray"))))
  expect_no_error(PlotNetworkU(network = unagi, vertexLabels = NA,
                               geneColorMapping = data.frame(gene = c("gene1", "gene2", "gene3", "gene4", "geneA", "geneB", "geneC", "geneD"), 
                                                             color = c("red", "green", "blue", "orange", "purple", "yellow", "black", "gray"))))
  expect_no_error(PlotNetworkU(network = unagi,  
                               vertexLabels = c("gene1", "gene2", "gene3", "gene4", "geneA", "geneB", "geneC", "geneD"),
                               geneColorMapping = NULL))
  expect_no_error(PlotNetworkU(network = unagi, vertexLabels = NA,
                               geneColorMapping = NULL))
})