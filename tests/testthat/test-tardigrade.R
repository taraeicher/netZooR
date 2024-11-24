# unit-tests for TARDIGRADE
context("test TARDIGRADE functions")

test_that("[TARDIGRADE] DoRaMPSetup() function yields expected results", {
  # Versions that don't exist return an error.
  #expect_error(DoRaMPSetup(version = "myFakeVersion"),
  #             "ERROR: myFakeVersion is not available from RaMP.")
  
  # Setup returns the correct structure when version is supplied.
  
  # Setup returns the correct structure when the latest (default) version is used.
  rampDB <- DoRaMPSetup()
})
test_that("[TARDIGRADE] BuildPathwayGraph() function yields expected results", {
  
  # Test when the data frame is the wrong format.
  errorMsgFormat <- paste("ERROR: analyteHasPathway must be a data frame",
                    "with character columns 'Analyte' and 'Pathway'")
  expect_error(BuildPathwayGraph(analyteHasPathway = ""), errorMsgFormat)
  expect_error(BuildPathwayGraph(analyteHasPathway = data.frame(apples = c("1","2","3"),
                                                                oranges = c("1","2","3"))), errorMsgFormat)
  expect_error(BuildPathwayGraph(analyteHasPathway = data.frame(Analyte = c("1","2","3"),
                                                                Pathway = c("1","2","3"),
                                                                SomethingElse = c("1","2","3"))), errorMsgFormat)
  expect_error(BuildPathwayGraph(analyteHasPathway = data.frame(Analyte = c(1,2,3),
                                                                Pathway = c(1,2,3))), errorMsgFormat)
  expect_error(BuildPathwayGraph(analyteHasPathway = data.frame(Analyte = c("1","2","3"),
                                                                Pathway = c("1","2","3")), pathwaySizeLimit = "woohoo"), 
               "ERROR: pathwaySizeLimit must be numeric")
  
  # Test when the data frame is empty. In this case, there should be no graph.
  emptyGraphMapping <- data.frame(Analyte = as.character(c()),
                                Pathway = as.character(c()))
  emptyGraph <- methods::new("TARDIGRADE_KnowledgeGraph", 
                             adjacencyMatrix = new("dgCMatrix"),
                             edgeDict = list(), clusters = list(), 
                             analytePathwayMapping = emptyGraphMapping,
                             clusterPathwayDict = list())
  resEmpty <- BuildPathwayGraph(analyteHasPathway = emptyGraphMapping)
  expect_equal(resEmpty@adjacencyMatrix, emptyGraph@adjacencyMatrix)
  expect_equal(resEmpty@edgeDict, emptyGraph@edgeDict)
  expect_equal(resEmpty@analytePathwayMapping, emptyGraph@analytePathwayMapping)
  
  # Test when there is only one mapping. In this case, there should be no graph.
  emptyGraph@adjacencyMatrix <- as(matrix(0, nrow = 1), "dgCMatrix")
  rownames(emptyGraph@adjacencyMatrix) <- "a1"
  colnames(emptyGraph@adjacencyMatrix) <- "a1"
  emptyGraph@analytePathwayMapping <- data.frame(Analyte = c("a1"),
                                                 Pathway = c("p1"))
  resEmpty <- BuildPathwayGraph(analyteHasPathway = data.frame(Analyte = c("a1"),
                                                               Pathway = c("p1"))) 
  expect_equal(resEmpty@adjacencyMatrix, emptyGraph@adjacencyMatrix)
  expect_equal(resEmpty@edgeDict, emptyGraph@edgeDict)
  expect_equal(resEmpty@analytePathwayMapping, emptyGraph@analytePathwayMapping)
  
  # Test when pathways are disjoint. In this case, there should be no graph.
  emptyGraph@analytePathwayMapping <- data.frame(Analyte = c("a1", "a2"),
                                                 Pathway = c("p1", "p2"))
  emptyGraph@adjacencyMatrix <- as(matrix(rep(0,4), nrow = 2), "dgCMatrix")
  rownames(emptyGraph@adjacencyMatrix) <- c("a1", "a2")
  colnames(emptyGraph@adjacencyMatrix) <- c("a1", "a2")
  resEmpty <- BuildPathwayGraph(analyteHasPathway = data.frame(Analyte = c("a1", "a2"),
                                                               Pathway = c("p1", "p2")))
  expect_equal(resEmpty@adjacencyMatrix, emptyGraph@adjacencyMatrix)
  expect_equal(resEmpty@edgeDict, emptyGraph@edgeDict)
  expect_equal(resEmpty@analytePathwayMapping, emptyGraph@analytePathwayMapping)
  
  # Test when the input is valid, but the pathway size limit is set to 0. In this
  # case, there should be no graph.
  overlappingAdjMat <- as(matrix(c(rep(1, 4)), nrow = 2), "dgCMatrix")
  diag(overlappingAdjMat) <- 0
  colnames(overlappingAdjMat) <- c("a1", "a2")
  rownames(overlappingAdjMat) <- colnames(overlappingAdjMat)
  overlappingEdgeDict <- list("a1|a2" = "p1", "a2|a1" = "p1")
  overlappingMapping <- data.frame(Analyte = c("a1", "a2"),
                                   Pathway = c("p1", "p1"))
  overlappingPathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                      adjacencyMatrix = overlappingAdjMat,
                                      edgeDict = overlappingEdgeDict, 
                                      clusters = list(),
                                      analytePathwayMapping = overlappingMapping,
                                      clusterPathwayDict = list())
  resOverlapping <- BuildPathwayGraph(analyteHasPathway = overlappingMapping, pathwaySizeLimit = 0)
  expect_equal(resOverlapping@adjacencyMatrix, new("dgCMatrix"))
  expect_equal(resOverlapping@edgeDict, emptyGraph@edgeDict)
  expect_equal(resOverlapping@analytePathwayMapping, emptyGraphMapping)
  
  # Test when the input is valid, but the pathway size limit is set to 1. In this
  # case, there should be no graph.
  overlappingAdjMat <- as(matrix(c(rep(1, 4)), nrow = 2), "dgCMatrix")
  diag(overlappingAdjMat) <- 0
  colnames(overlappingAdjMat) <- c("a1", "a2")
  rownames(overlappingAdjMat) <- colnames(overlappingAdjMat)
  overlappingEdgeDict <- list("a1|a2" = "p1", "a2|a1" = "p1")
  overlappingMapping <- data.frame(Analyte = c("a1", "a2"),
                                   Pathway = c("p1", "p1"))
  overlappingPathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                      adjacencyMatrix = overlappingAdjMat,
                                      edgeDict = overlappingEdgeDict, 
                                      clusters = list(),
                                      analytePathwayMapping = overlappingMapping,
                                      clusterPathwayDict = list())
  resOverlapping <- BuildPathwayGraph(analyteHasPathway = overlappingMapping, pathwaySizeLimit = 1)
  expect_equal(resOverlapping@adjacencyMatrix, new("dgCMatrix"))
  expect_equal(resOverlapping@edgeDict, emptyGraph@edgeDict)
  expect_equal(resOverlapping@analytePathwayMapping, emptyGraphMapping)
  
  # Test when the input is valid, but the pathway size limit is set to 2. In this
  # case, there should be no graph.
  overlappingAdjMat <- as(matrix(c(rep(1, 4)), nrow = 2), "dgCMatrix")
  diag(overlappingAdjMat) <- 0
  colnames(overlappingAdjMat) <- c("a1", "a2")
  rownames(overlappingAdjMat) <- colnames(overlappingAdjMat)
  overlappingEdgeDict <- list("a1|a2" = "p1", "a2|a1" = "p1")
  overlappingMapping <- data.frame(Analyte = c("a1", "a2"),
                                   Pathway = c("p1", "p1"))
  overlappingPathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                      adjacencyMatrix = overlappingAdjMat,
                                      edgeDict = overlappingEdgeDict, 
                                      clusters = list(),
                                      analytePathwayMapping = overlappingMapping,
                                      clusterPathwayDict = list())
  resOverlapping <- BuildPathwayGraph(analyteHasPathway = overlappingMapping, pathwaySizeLimit = 2)
  expect_equal(resOverlapping@adjacencyMatrix, new("dgCMatrix"))
  expect_equal(resOverlapping@edgeDict, emptyGraph@edgeDict)
  expect_equal(resOverlapping@analytePathwayMapping, emptyGraphMapping)
  
  # Test when there are analytes that both belong to the same pathway.
  overlappingAdjMat <- as(matrix(c(rep(1, 4)), nrow = 2), "dgCMatrix")
  diag(overlappingAdjMat) <- 0
  colnames(overlappingAdjMat) <- c("a1", "a2")
  rownames(overlappingAdjMat) <- colnames(overlappingAdjMat)
  overlappingEdgeDict <- list("a1|a2" = "p1", "a2|a1" = "p1")
  overlappingMapping <- data.frame(Analyte = c("a1", "a2"),
                                   Pathway = c("p1", "p1"))
  overlappingPathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                      adjacencyMatrix = overlappingAdjMat,
                                      edgeDict = overlappingEdgeDict, 
                                      clusters = list(),
                                      analytePathwayMapping = overlappingMapping,
                                      clusterPathwayDict = list())
  resOverlapping <- BuildPathwayGraph(analyteHasPathway = overlappingMapping)
  expect_equal(resOverlapping@adjacencyMatrix, overlappingPathways@adjacencyMatrix)
  expect_equal(resOverlapping@edgeDict, overlappingPathways@edgeDict[names(resOverlapping@edgeDict)])
  expect_equal(resOverlapping@analytePathwayMapping, overlappingPathways@analytePathwayMapping)
  
  # Test when two small pathways overlap with a large pathway.
  lgPathwaySize <- 100
  smallLargeAdjMat <- as(matrix(c(rep(1/(lgPathwaySize-1), (lgPathwaySize + 2)^2)), nrow = lgPathwaySize + 2), "dgCMatrix")
  diag(smallLargeAdjMat) <- 0
  colnames(smallLargeAdjMat) <- c(paste0(rep("a", lgPathwaySize), 1:lgPathwaySize), "a#", "a*")
  rownames(smallLargeAdjMat) <- colnames(smallLargeAdjMat)
  smSharedMat <- 1
  smallLargeAdjMat[,"a#"] <- 0
  smallLargeAdjMat["a#",] <- 0
  smallLargeAdjMat[,"a*"] <- 0
  smallLargeAdjMat["a*",] <- 0
  smallLargeAdjMat["a1","a#"] <- smSharedMat
  smallLargeAdjMat["a#","a1"] <- smSharedMat
  smallLargeAdjMat["a2","a*"] <- smSharedMat
  smallLargeAdjMat["a*","a2"] <- smSharedMat
  smallLargeAdjMat["a1","a*"] <- 0
  smallLargeAdjMat["a*","a1"] <- 0
  smallLargeAdjMat["a2","a#"] <- 0
  smallLargeAdjMat["a#","a2"] <- 0
  smallLargeEdgeDict <- as.list(rep("p3", lgPathwaySize^2))
  names(smallLargeEdgeDict) <- paste(paste0("a", rep(1:lgPathwaySize, lgPathwaySize)),
                                     unlist(lapply(1:lgPathwaySize, function(i){return(paste0("a", rep(i, lgPathwaySize)))})),
                                     sep = "|")
  smallLargeEdgeDict["a1|a#"] <- "p1"
  smallLargeEdgeDict["a#|a1"] <- "p1"
  smallLargeEdgeDict["a2|a*"] <- "p2"
  smallLargeEdgeDict["a*|a2"] <- "p2"
  smallLargeMapping <- data.frame(Analyte = c("a1", "a#", "a2", "a*", 
                                               paste0(rep("a", lgPathwaySize), 1:lgPathwaySize)),
                                   Pathway = c("p1", "p1", "p2", "p2",
                                               rep("p3", lgPathwaySize)))
  smallLargePathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                      adjacencyMatrix = smallLargeAdjMat,
                                      edgeDict = smallLargeEdgeDict, 
                                      clusters = list(),
                                      analytePathwayMapping = smallLargeMapping,
                                      clusterPathwayDict = list())
  resSmallLarge <- BuildPathwayGraph(analyteHasPathway = smallLargeMapping)
  expect_equal(resSmallLarge@adjacencyMatrix, 
               smallLargePathways@adjacencyMatrix[rownames(resSmallLarge@adjacencyMatrix),
                                                  colnames(resSmallLarge@adjacencyMatrix)])
  expect_equal(resSmallLarge@edgeDict, smallLargePathways@edgeDict[names(resSmallLarge@edgeDict)])
  expect_equal(resSmallLarge@analytePathwayMapping, smallLargePathways@analytePathwayMapping)
  
  # Test when two small pathways overlap with a large pathway, but the large pathway
  # is excluded from the analysis.
  lgPathwaySize <- 100
  smallLargeAdjMat <- as(matrix(c(rep(0, 16)), nrow = 4), "dgCMatrix")
  diag(smallLargeAdjMat) <- 0
  colnames(smallLargeAdjMat) <- c("a1", "a#", "a2", "a*")
  rownames(smallLargeAdjMat) <- colnames(smallLargeAdjMat)
  smSharedMat <- 1
  smallLargeAdjMat[,"a#"] <- 0
  smallLargeAdjMat["a#",] <- 0
  smallLargeAdjMat[,"a*"] <- 0
  smallLargeAdjMat["a*",] <- 0
  smallLargeAdjMat["a1","a#"] <- smSharedMat
  smallLargeAdjMat["a#","a1"] <- smSharedMat
  smallLargeAdjMat["a2","a*"] <- smSharedMat
  smallLargeAdjMat["a*","a2"] <- smSharedMat
  smallLargeAdjMat["a1","a*"] <- 0
  smallLargeAdjMat["a*","a1"] <- 0
  smallLargeAdjMat["a2","a#"] <- 0
  smallLargeAdjMat["a#","a2"] <- 0
  smallLargeEdgeDict <- list()
  smallLargeEdgeDict["a1|a#"] <- "p1"
  smallLargeEdgeDict["a#|a1"] <- "p1"
  smallLargeEdgeDict["a2|a*"] <- "p2"
  smallLargeEdgeDict["a*|a2"] <- "p2"
  smallLargeMapping <- data.frame(Analyte = c("a1", "a#", "a2", "a*"),
                                  Pathway = c("p1", "p1", "p2", "p2"))
  smallLargePathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                     adjacencyMatrix = smallLargeAdjMat,
                                     edgeDict = smallLargeEdgeDict, 
                                     clusters = list(),
                                     analytePathwayMapping = smallLargeMapping,
                                     clusterPathwayDict = list())
  resSmallLarge <- BuildPathwayGraph(analyteHasPathway = smallLargeMapping, pathwaySizeLimit = lgPathwaySize)
  expect_equal(resSmallLarge@adjacencyMatrix, 
               smallLargePathways@adjacencyMatrix[rownames(resSmallLarge@adjacencyMatrix),
                                                  colnames(resSmallLarge@adjacencyMatrix)])
  expect_equal(resSmallLarge@edgeDict, smallLargePathways@edgeDict[names(resSmallLarge@edgeDict)])
  expect_equal(resSmallLarge@analytePathwayMapping, smallLargePathways@analytePathwayMapping)
  
  # Test when there are two disjoint overlap groups, with one pair of analytes sharing
  # multiple pathways.
  disjointAdjMat <- as(matrix(c(rep(0, 36)), nrow = 6), "dgCMatrix")
  colnames(disjointAdjMat) <- paste0("a", 1:6)
  rownames(disjointAdjMat) <- colnames(disjointAdjMat)
  disjointAdjMat["a1","a2"] <- 1
  disjointAdjMat["a1","a3"] <- 1
  disjointAdjMat["a2","a3"] <- 1
  disjointAdjMat["a4","a5"] <- 1
  disjointAdjMat["a5","a6"] <- 1
  disjointAdjMat["a4","a6"] <- 1
  disjointAdjMat["a2","a1"] <- 1
  disjointAdjMat["a3","a1"] <- 1
  disjointAdjMat["a3","a2"] <- 1
  disjointAdjMat["a5","a4"] <- 1
  disjointAdjMat["a6","a5"] <- 1
  disjointAdjMat["a6","a4"] <- 1
  disjointEdgeDict <- as.list(c(rep(paste0("p", 1:5), 2), rep("p6;p7", 2)))
  names(disjointEdgeDict) <- c("a1|a2", "a2|a3", "a1|a3", "a4|a5", "a5|a6",
                               "a2|a1", "a3|a2", "a3|a1", "a5|a4", "a6|a5", "a4|a6", "a6|a4")
  disjointMapping <- data.frame(Analyte = c("a1", "a2", "a2", "a3", "a1", "a3", "a4", "a5", "a5", "a6", "a4", "a6", "a4", "a6"),
                         Pathway = c("p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4", "p5", "p5", "p6", "p6", "p7", "p7"))
  disjointPathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                     adjacencyMatrix = disjointAdjMat,
                                     edgeDict = disjointEdgeDict, 
                                     clusters = list(),
                                     analytePathwayMapping = disjointMapping,
                                     clusterPathwayDict = list())
  resDisjoint <- BuildPathwayGraph(analyteHasPathway = disjointMapping)
  expect_equal(resDisjoint@adjacencyMatrix, 
               disjointPathways@adjacencyMatrix[rownames(resDisjoint@adjacencyMatrix),
                                                  colnames(resDisjoint@adjacencyMatrix)])
  expect_equal(resDisjoint@edgeDict, disjointPathways@edgeDict[names(resDisjoint@edgeDict)])
  expect_equal(resDisjoint@analytePathwayMapping, disjointPathways@analytePathwayMapping)
})
test_that("[TARDIGRADE] FindClusters() function yields expected results",{
  
  # If the graph is empty or has no edges, an error is thrown.
  emptyGraphMapping <- data.frame(Analyte = as.character(c()),
                                  Pathway = as.character(c()))
  emptyGraph <- methods::new("TARDIGRADE_KnowledgeGraph", 
                             adjacencyMatrix = new("dgCMatrix"),
                             edgeDict = list(), clusters = list(), 
                             analytePathwayMapping = emptyGraphMapping,
                             clusterPathwayDict = list())
  expect_error(FindClusters(knowledgeGraph = emptyGraph, clusterCount = 2),
               "ERROR: The knowledge graph has no edges.")
  emptyGraph@adjacencyMatrix <- as(matrix(0, nrow = 1), "dgCMatrix")
  rownames(emptyGraph@adjacencyMatrix) <- "a1"
  colnames(emptyGraph@adjacencyMatrix) <- "a1"
  emptyGraph@analytePathwayMapping <- data.frame(Analyte = c("a1"),
                                                 Pathway = c("p1"))
  expect_error(FindClusters(knowledgeGraph = emptyGraph, clusterCount = 2),
               "ERROR: The knowledge graph has no edges.")
  
  # If the number of clusters is greater than the number of nodes, an 
  # error is thrown.
  overlappingAdjMat <- as(matrix(c(rep(1, 4)), nrow = 2), "dgCMatrix")
  diag(overlappingAdjMat) <- 0
  colnames(overlappingAdjMat) <- c("a1", "a2")
  rownames(overlappingAdjMat) <- colnames(overlappingAdjMat)
  overlappingEdgeDict <- list("a1|a2" = "p1", "a2|a1" = "p1")
  overlappingMapping <- data.frame(Analyte = c("a1", "a2"),
                                   Pathway = c("p1", "p1"))
  overlappingPathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                      adjacencyMatrix = overlappingAdjMat,
                                      edgeDict = overlappingEdgeDict, 
                                      clusters = list(),
                                      analytePathwayMapping = overlappingMapping,
                                      clusterPathwayDict = list())
  expect_error(FindClusters(knowledgeGraph = overlappingPathways, clusterCount = 3),
               "ERROR: The number of clusters must be smaller than the number of nodes.")
  
  # If the graph has only one edge, a single cluster is returned.
  clusters1Edge <- FindClusters(knowledgeGraph = overlappingPathways, clusterCount = 1)
  clusters <- list("1" = c("a1", "a2"))
  clusterPathwayDict <- list("1" = "p1")
  expect_equal(clusters1Edge@adjacencyMatrix, overlappingPathways@adjacencyMatrix)
  expect_equal(clusters1Edge@edgeDict, overlappingPathways@edgeDict[names(clusters1Edge@edgeDict)])
  expect_equal(clusters1Edge@analytePathwayMapping, overlappingPathways@analytePathwayMapping)
  expect_equal(clusters1Edge@clusters, clusters)
  expect_equal(clusters1Edge@clusterPathwayDict, clusterPathwayDict)
  
  # Disjoint cliques are returned correctly.
  disjointAdjMat <- as(matrix(c(rep(0, 14^2)), nrow = 14), "dgCMatrix")
  colnames(disjointAdjMat) <- c(paste0("a", 1:14))
  rownames(disjointAdjMat) <- colnames(disjointAdjMat)
  disjointAdjMat[1:5, 1:5] <- 1/4
  disjointAdjMat[6:9, 6:9] <- 1/3
  disjointAdjMat[10:14, 10:14] <- 1/4
  diag(disjointAdjMat) <- 0
  disjointEdgeDict <- list("a1|a2" = "p1;p2;p3", "a1|a3" = "p1;p2;p3", "a1|a4" = "p1;p2;p3",
                           "a1|a5" = "p1;p2;p3", "a2|a1" = "p1;p2;p3", "a2|a3" = "p1;p2;p3", 
                           "a2|a4" = "p1;p2;p3", "a2|a5" = "p1;p2;p3", "a3|a1" = "p1;p2;p3",
                           "a3|a2" = "p1;p2;p3", "a3|a4" = "p1;p2;p3", "a3|a5" = "p1;p2;p3",
                           "a4|a1" = "p1;p2;p3", "a4|a2" = "p1;p2;p3", "a4|a3" = "p1;p2;p3",
                           "a4|a5" = "p1;p2;p3", "a5|a1" = "p1;p2;p3", "a5|a2" = "p1;p2;p3",
                           "a5|a3" = "p1;p2;p3", "a5|a4" = "p1;p2;p3",
                           "a6|a7" = "p4;p5", "a6|a8" = "p4;p5", "a6|a9" = "p4;p5",
                           "a7|a6" = "p4;p5", "a7|a8" = "p4;p5","a7|a9" = "p4;p5",
                           "a8|a6" = "p4;p5", "a8|a7" = "p4;p5", "a8|a9" = "p4;p5",
                           "a9|a6" = "p4;p5", "a9|a7" = "p4;p5", "a9|a8" = "p4;p5",
                           "a10|a11" = "p6", "a10|a12" = "p6", "a10|a13" = "p6",
                           "a10|a14" = "p6", "a11|a10" = "p6", "a11|a12" = "p6", 
                           "a11|a13" = "p6", "a11|a14" = "p6", "a12|a10" = "p6",
                           "a12|a11" = "p6", "a12|a13" = "p6", "a12|a14" = "p6",
                           "a13|a10" = "p6", "a13|a11" = "p6", "a13|a12" = "p6",
                           "a13|a14" = "p6", "a14|a10" = "p6", "a14|a11" = "p6",
                           "a14|a12" = "p6", "a14|a13" = "p6")
  disjointMapping <- data.frame(Analyte = c(paste0("a", 1:5), paste0("a", 1:5),
                                               paste0("a", 1:5), paste0("a", 6:9),
                                               paste0("a", 6:9), paste0("a", 10:14)),
                                   Pathway = c(rep("p1", 5), rep("p2", 5), 
                                               rep("p3", 5), rep("p4", 4),
                                               rep("p5", 4), rep("p6", 5)))
  disjointPathways <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                      adjacencyMatrix = disjointAdjMat,
                                      edgeDict = disjointEdgeDict, 
                                      clusters = list(),
                                      analytePathwayMapping = disjointMapping,
                                      clusterPathwayDict = list())
  clustersDisjoint <- FindClusters(knowledgeGraph = disjointPathways, clusterCount = 3)
  clusters <- list("1" = paste0("a", 1:5), "2" = paste0("a", 6:9), "3" = paste0("a", 10:14))
  clusterPathwayDict <- list("1" = "p1;p2;p3", "2" = "p4;p5", "3" = "p6")
  expect_equal(clustersDisjoint@adjacencyMatrix, disjointPathways@adjacencyMatrix)
  expect_equal(clustersDisjoint@edgeDict, disjointPathways@edgeDict[names(clustersDisjoint@edgeDict)])
  expect_equal(clustersDisjoint@analytePathwayMapping, disjointPathways@analytePathwayMapping)
  expect_equal(clustersDisjoint@clusters, clusters)
  expect_equal(clustersDisjoint@clusterPathwayDict, clusterPathwayDict)
  
  # Almost-disjoint cliques are returned correctly.
  almostDisjointPathways <- disjointPathways
  almostDisjointAdj <- disjointPathways@adjacencyMatrix
  almostDisjointAdj[1:5, 6:9] <- 1/8
  almostDisjointAdj[6:9, 1:5] <- 1/8
  almostDisjointAdj[1:5, 10:14] <- 1/9
  almostDisjointAdj[10:14, 1:5] <- 1/9
  almostDisjointPathways@adjacencyMatrix <- almostDisjointAdj
  almostDisjointEdgeDict <- disjointPathways@edgeDict
  edgesInFirstLargePathway1 <- unlist(lapply(1:5, function(i){
    source <- paste0("a", i)
    return(unlist(lapply(setdiff(1:5, i), function(j){
      print(paste(i,j))
      return(paste(source, paste0("a", j), sep = "|"))
    })))
  }))
  almostDisjointEdgeDict[edgesInFirstLargePathway1] <- "p1;p2;p3;p7;p8"
  edgesInFirstLargePathway2 <- unlist(lapply(1:5, function(i){
    source <- paste0("a", i)
    return(unlist(lapply(6:9, function(j){
      print(paste(i,j))
      return(paste(source, paste0("a", j), sep = "|"))
    })))
  }))
  almostDisjointEdgeDict[edgesInFirstLargePathway2] <- "p7"
  edgesInFirstLargePathway3 <- unlist(lapply(6:9, function(i){
    source <- paste0("a", i)
    return(unlist(lapply(1:5, function(j){
      print(paste(i,j))
      return(paste(source, paste0("a", j), sep = "|"))
    })))
  }))
  almostDisjointEdgeDict[edgesInFirstLargePathway3] <- "p7"
  edgesInFirstLargePathway4 <- unlist(lapply(6:9, function(i){
    source <- paste0("a", i)
    return(unlist(lapply(setdiff(6:9, i), function(j){
      print(paste(i,j))
      return(paste(source, paste0("a", j), sep = "|"))
    })))
  }))
  almostDisjointEdgeDict[edgesInFirstLargePathway4] <- "p4;p5;p7"
  edgesInSecondLargePathway1 <- unlist(lapply(1:5, function(i){
    source <- paste0("a", i)
    return(unlist(lapply(setdiff(10:14, i), function(j){
      print(paste(i,j))
      return(paste(source, paste0("a", j), sep = "|"))
    })))
  }))
  almostDisjointEdgeDict[edgesInSecondLargePathway1] <- "p8"
  edgesInSecondLargePathway2 <- unlist(lapply(10:14, function(i){
    source <- paste0("a", i)
    return(unlist(lapply(setdiff(1:5, i), function(j){
      print(paste(i,j))
      return(paste(source, paste0("a", j), sep = "|"))
    })))
  }))
  almostDisjointEdgeDict[edgesInSecondLargePathway2] <- "p8"
  edgesInSecondLargePathway2 <- unlist(lapply(10:14, function(i){
    source <- paste0("a", i)
    return(unlist(lapply(setdiff(10:14, i), function(j){
      print(paste(i,j))
      return(paste(source, paste0("a", j), sep = "|"))
    })))
  }))
  almostDisjointEdgeDict[edgesInSecondLargePathway2] <- "p6;p8"
  almostDisjointPathways@edgeDict <- almostDisjointEdgeDict
  almostDisjointMapping <- rbind(disjointPathways@analytePathwayMapping,
                                 data.frame(Analyte = c(paste0("a", 1:5), paste0("a", 1:5),
                                                        paste0("a", 6:9), paste0("a", 10:14)),
                                            Pathway = c(rep("p7", 5), rep("p8", 5), 
                                                        rep("p7", 4), rep("p8", 5))))
  almostDisjointPathways@analytePathwayMapping <- almostDisjointMapping
  clusterPathwayDict <- list("1" = "p1;p2;p3;p7;p8", "2" = "p4;p5;p7", "3" = "p6;p8")
  almostDisjointPathways@clusterPathwayDict <- clusterPathwayDict
  clustersAlmostDisjoint <- FindClusters(knowledgeGraph = almostDisjointPathways, clusterCount = 3)
  expect_equal(clustersAlmostDisjoint@adjacencyMatrix, almostDisjointPathways@adjacencyMatrix)
  expect_equal(clustersAlmostDisjoint@edgeDict, almostDisjointPathways@edgeDict[names(clustersAlmostDisjoint@edgeDict)])
  expect_equal(clustersAlmostDisjoint@analytePathwayMapping, almostDisjointPathways@analytePathwayMapping)
  expect_equal(clustersAlmostDisjoint@clusters, clusters)
  expect_equal(clustersAlmostDisjoint@clusterPathwayDict, clusterPathwayDict)
})
test_that("[TARDIGRADE] ModifyGraph() function yields expected results",{
  
  # Create two disjoint overlap groups, with one pair of analytes sharing
  # multiple pathways.
  disjointMapping <- data.frame(Analyte = c("a1", "a2", "a2", "a3", "a1", "a3", "a4", "a5", "a5", "a6", "a4", "a6", "a4", "a6"),
                                Pathway = c("p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4", "p5", "p5", "p6", "p6", "p7", "p7"))
  resDisjoint <- BuildPathwayGraph(analyteHasPathway = disjointMapping)
  
  # If there is no overlap between the network and the knowledge graph, an
  # error is thrown.
  network <- data.frame(Source = c("m1", "m2", "m3"),
                        Target = c("g1", "g2", "g3"))
  expect_error(ModifyGraph(knowledgeGraph = resDisjoint, network = network,
                           sourceAnalytes = c("m1", "m2", "m3"), 
                           targetAnalytes = c("g1", "g2", "g3")),
               "ERROR: The input network and knowledge graph do not overlap")
  
  # If the network analytes are not a subset of the source or target analytes,
  # an error is thrown.
  additionalMapping <- data.frame(Analyte = c("g1", "g2", "g2", "g3", "m1", "m3", "m4", "m5", "m5", "m6", "m4", "m6", "m4", "m6"),
                                Pathway = c("p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4", "p5", "p5", "p6", "p6", "p7", "p7"))
  resAdditional <- BuildPathwayGraph(analyteHasPathway = additionalMapping)
  expect_error(ModifyGraph(knowledgeGraph = resAdditional, network = network,
                           sourceAnalytes = c("i1", "i2", "i3"), 
                           targetAnalytes = c("g1", "g2", "g3")),
               "ERROR: Source or target analytes do not match the input network")
  expect_error(ModifyGraph(knowledgeGraph = resAdditional, network = network,
                           sourceAnalytes = c("m1", "m2", "m3"), 
                           targetAnalytes = c("z1", "z2", "z3")),
               "ERROR: Source or target analytes do not match the input network")
  
  # If one of the inputs is empty, an error is thrown.
  network <- data.frame(Source = c("a1", "a2", "a3"),
                        Target = c("a4", "a5", "a6"))
  resEmpty <- resDisjoint
  resEmpty@adjacencyMatrix <- new("dgCMatrix")
  resEmpty@analytePathwayMapping <- data.frame(Analyte = c(), Pathway = c())
  networkEmpty <- data.frame(Source = c(), Target = c())
  sourceEmpty <- c()
  targetEmpty <- c()
  expect_error(ModifyGraph(knowledgeGraph = resEmpty, network = network,
                           sourceAnalytes = c("a1", "a2", "a3"), 
                           targetAnalytes = c("a4", "a5", "a6")),
               "ERROR: The input network and knowledge graph do not overlap")
  expect_error(ModifyGraph(knowledgeGraph = resDisjoint, network = networkEmpty,
                           sourceAnalytes = c("a1", "a2", "a3"), 
                           targetAnalytes = c("a4", "a5", "a6")),
               "ERROR: The input network is empty")
  expect_error(ModifyGraph(knowledgeGraph = resDisjoint, network = network,
                           sourceAnalytes = sourceEmpty, 
                           targetAnalytes = c("a4", "a5", "a6")),
               "ERROR: Source or target analytes do not match the input network")
  expect_error(ModifyGraph(knowledgeGraph = resDisjoint, network = network,
                           sourceAnalytes = c("a1", "a2", "a3"), 
                           targetAnalytes = targetEmpty),
               "ERROR: Source or target analytes do not match the input network")
  
  # If one of the inputs is of the wrong format, an error is thrown.
  expect_error(ModifyGraph(knowledgeGraph = "hi", network = network,
                           sourceAnalytes = c("a1", "a2", "a3"), 
                           targetAnalytes = c("a4", "a5", "a6")),
               "ERROR: Input formats must be: TARDIGRADE_KnowledgeGraph object (knowledgeGraph), data frame (network), and character (sourceAnalytes and targetAnalytes)")
  expect_error(ModifyGraph(knowledgeGraph = resDisjoint, network = "hi",
                           sourceAnalytes = c("a1", "a2", "a3"), 
                           targetAnalytes = c("a4", "a5", "a6")),
               "ERROR: Input formats must be: TARDIGRADE_KnowledgeGraph object (knowledgeGraph), data frame (network), and character (sourceAnalytes and targetAnalytes)")
  expect_error(ModifyGraph(knowledgeGraph = resDisjoint, network = network,
                           sourceAnalytes = 1, 
                           targetAnalytes = c("a4", "a5", "a6")),
               "ERROR: Input formats must be: TARDIGRADE_KnowledgeGraph object (knowledgeGraph), data frame (network), and character (sourceAnalytes and targetAnalytes)")
  expect_error(ModifyGraph(knowledgeGraph = resDisjoint, network = network,
                           sourceAnalytes = c("a1", "a2", "a3"), 
                           targetAnalytes = 1),
               "ERROR: Input formats must be: TARDIGRADE_KnowledgeGraph object (knowledgeGraph), data frame (network), and character (sourceAnalytes and targetAnalytes)")
  
  # If there are only edges added, they are added correctly, even with partial overlap.
  trueModifiedMatrix <- as(matrix(c(0,1,1,1,0,0,
                                    1,0,1,0,0,0,
                                    1,1,0,0,0,0,
                                    1,0,0,0,1,1,
                                    0,0,0,1,0,1,
                                    0,0,0,1,1,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a2", "a1", "a4", "a5", "a4", "a1"),
                                  targetAnalytes = c("a2", "a3", "a3", "a5", "a6", "a6", "a4"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resDisjoint, network = networkAdditional,
                               sourceAnalytes = c("a1", "a2", "a3", "a4", "a5"), 
                               targetAnalytes = c("a2", "a3", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a4", "a5", "a4", "a1"),
                                  targetAnalytes = c("a2", "a5", "a6", "a6", "a4"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resDisjoint, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  
  # If there are only edges removed, they are removed correctly, even with partial overlap.
  trueModifiedMatrix <- as(matrix(c(0,0,1,0,0,0,
                                    0,0,1,0,0,0,
                                    1,1,0,0,0,0,
                                    0,0,0,0,1,1,
                                    0,0,0,1,0,1,
                                    0,0,0,1,1,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkRemoval <- data.frame(sourceAnalytes = c("a2", "a1", "a4", "a5", "a4"),
                                  targetAnalytes = c("a3", "a3", "a5", "a6", "a6"))
  modifiedGraphRemoval <- ModifyGraph(knowledgeGraph = resDisjoint, network = networkRemoval,
                                         sourceAnalytes = c("a1", "a2", "a3", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a3", "a4", "a5", "a6"))
  expect_equal(modifiedGraphRemoval, trueModifiedMatrix)
  networkRemoval <- data.frame(sourceAnalytes = c("a4", "a5", "a4"),
                                  targetAnalytes = c("a5", "a6", "a6"))
  modifiedGraphRemoval <- ModifyGraph(knowledgeGraph = resDisjoint, network = networkRemoval,
                                         sourceAnalytes = c("a1", "a2", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a4", "a5", "a6"))
  expect_equal(modifiedGraphRemoval, trueModifiedMatrix)
  
  # If there are only edges with modified weights, they are modified correctly, even with partial overlap.
  lowWeightMapping <- data.frame(Analyte = c("a1", "a2", "a3", "a3", "a4", "a5", "a5", "a6", "a4", "a6", "a4", "a6"),
                                Pathway = c("p1", "p1", "p1", "p3", "p4", "p4", "p5", "p5", "p6", "p6", "p7", "p7"))
  resLowWeight <- BuildPathwayGraph(analyteHasPathway = lowWeightMapping)
  trueModifiedMatrix <- as(matrix(c(0,1,1,0,0,0,
                                    1,0,1,0,0,0,
                                    1,1,0,0,0,0,
                                    0,0,0,0,1,1,
                                    0,0,0,1,0,1,
                                    0,0,0,1,1,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a2", "a1", "a4", "a5", "a4"),
                                  targetAnalytes = c("a2", "a3", "a3", "a5", "a6", "a6"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resLowWeight, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a3", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a3", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  trueModifiedMatrix <- as(matrix(c(0,1,0.5,0,0,0,
                                    1,0,0.5,0,0,0,
                                    0.5,0.5,0,0,0,0,
                                    0,0,0,0,1,1,
                                    0,0,0,1,0,1,
                                    0,0,0,1,1,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a4", "a5", "a4"),
                                  targetAnalytes = c("a2", "a5", "a6", "a6"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resLowWeight, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  
  # If addition, removal, and modification are all done, they are done correctly, even with partial overlap.
  # If there are only edges with modified weights, they are modified correctly, even with partial overlap.
  trueModifiedMatrix <- as(matrix(c(0,1,1,1,0,0,
                                    1,0,1,0,0,0,
                                    1,1,0,0,0,0,
                                    1,0,0,0,1,1,
                                    0,0,0,1,0,0,
                                    0,0,0,1,0,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a2", "a1", "a4", "a4", "a1"),
                                  targetAnalytes = c("a2", "a3", "a3", "a5", "a6", "a4"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resLowWeight, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a3", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a3", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  trueModifiedMatrix <- as(matrix(c(0,1,0.5,1,0,0,
                                    1,0,0.5,0,0,0,
                                    0.5,0.5,0,0,0,0,
                                    1,0,0,0,1,1,
                                    0,0,0,1,0,0,
                                    0,0,0,1,0,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a4", "a4", "a1"),
                                  targetAnalytes = c("a2", "a5", "a6", "a4"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resLowWeight, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  
  # If all edges are modified, it is done correctly, even with partial overlap.
  lowWeightAllMapping <- data.frame(Analyte = c("a1", "a2", "a3", "a3", "a4", "a5", "a6"),
                                 Pathway = c("p1", "p1", "p1", "p3", "p4", "p4", "p4"))
  resLowWeightAll <- BuildPathwayGraph(analyteHasPathway = lowWeightAllMapping)
  trueModifiedMatrix <- as(matrix(c(0,1,1,0,0,0,
                                    1,0,1,0,0,0,
                                    1,1,0,0,0,0,
                                    0,0,0,0,1,1,
                                    0,0,0,1,0,1,
                                    0,0,0,1,1,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a2", "a1", "a4", "a5", "a4"),
                                  targetAnalytes = c("a2", "a3", "a3", "a5", "a6", "a6"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resLowWeightAll, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a3", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a3", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  trueModifiedMatrix <- as(matrix(c(0,1,0.5,0,0,0,
                                    1,0,0.5,0,0,0,
                                    0.5,0.5,0,0,0,0,
                                    0,0,0,0,1,1,
                                    0,0,0,1,0,1,
                                    0,0,0,1,1,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1", "a4", "a5", "a4"),
                                  targetAnalytes = c("a2", "a5", "a6", "a6"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resLowWeightAll, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
  
  # If all edges are removed, it is done correctly.
  trueModifiedMatrix <- as(matrix(c(0,0,0,1,0,0,
                                    0,0,0,0,0,0,
                                    0,0,0,0,0,0,
                                    1,0,0,0,0,0,
                                    0,0,0,0,0,0,
                                    0,0,0,0,0,0), nrow = 6), "dgCMatrix")
  rownames(trueModifiedMatrix) <- paste0("a", rep(1:6))
  colnames(trueModifiedMatrix) <- rownames(trueModifiedMatrix)
  networkAdditional <- data.frame(sourceAnalytes = c("a1"),
                                  targetAnalytes = c("a4"))
  modifiedGraphAdditional <- ModifyGraph(knowledgeGraph = resLowWeightAll, network = networkAdditional,
                                         sourceAnalytes = c("a1", "a2", "a3", "a4", "a5"), 
                                         targetAnalytes = c("a2", "a3", "a4", "a5", "a6"))
  expect_equal(modifiedGraphAdditional, trueModifiedMatrix)
})
test_that("[TARDIGRADE] PlotKnowledgeGraph() function yields expected results",{
  # Test that an error is thrown when invalid values are supplied.
  
  # Test that an invalid state is returned when a file is open.
  
  # Test that a plot exists when it should exist.
  
  # Test that "OKAY" status is returned when it should be.
})
test_that("[TARDIGRADE] CalculateFruchtermanReingoldLayout() function yields expected results",{
  # Test that an error is thrown when invalid values are supplied.
  
  # Test that the correct format is returned.
})
test_that("[TARDIGRADE] DiffuseHeat() function yields expected results",{
  
  # If the inputs are of an incorrect format, an error should be thrown.
  
  # If one of the knowledge graphs is empty, an error should be thrown.
  
  # If the analyte measurements are empty, an error should be thrown.
  
  # If the analyte measurements do not match the graphs, an error should be thrown.
  
  # A standard case (connected graph, inverse can be taken) is handled appropriately.
  
  # Disconnected graphs are handled appropriately.
  
  # Singular graphs are handled appropriately.
  
  # Graphs without full rank are handled appropriately.
})
# GetNormalizedLaplacian
# DoPathwayAnalysis
# GetNormalizedEnrichmentScore