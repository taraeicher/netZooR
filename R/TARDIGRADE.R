#'  TARDIGRADE_KnowledgeGraph class
#'  This object type represents the list of results in a given directory.
#'  
#'  @name TARDIGRADE_KnowledgeGraph-class
#'  @rdname TARDIGRADE_KnowledgeGraph-class
#'  @exportClass TARDIGRADE_KnowledgeGraph
#'  @slot adjacencyMatrix The adjacency matrix for the graph.
#'  @slot edgeDict The dictionary mapping edges to pathways. This is a list.
#'  @slot clusters A named list of clusters.
methods::setClass(
  Class="TARDIGRADE_KnowledgeGraph",
  representation(adjacencyMatrix = "matrix",
                 edgeDict = "list",
                 clusters = "list",
                 analytePathwayMapping = "data.frame",
                 clusterPathwayDict = "list")
)

#' Builds the pathway graph from RaMP and stores in in the file specified by the user.
#' @param analytePathwayMappingFile Where the analyte-pathway mapping (input file) is stored.
#' @param pathwayGraphFile Where to store the pathway graph.
#' @returns An object of type TARDIGRADE_KnowledgeGraph.
#' @export
BuildPathwayGraph <- function(analytePathwayMappingFile = NULL, pathwayGraphFile = NULL){
  
  # Get tables needed to generate RaMP graph.
  analyteHasPathway <- read.csv(analytePathwayMappingFile)
  
  # Build the graph and the matching dictionary.
  # Initialize the graph and dict.
  G <- igraph::graph()
  edgePathwayDict < - list()
  
  # Obtain weights for each pathway.
  pathwayWeight <- list()
  for(pathway in analyteHasPathway$Pathway){
    
    # Get all analytes in the pathway.
    analytesInPathway  <- analyteHasPathway[which(analyteHasPathway$Pathway == pathway), "Analyte"]
    
    # Set the pathway's weight.
    pathwayWeight[pathway] <- 1.0 / length(analytesInPathway)
    
    # Add an edge for each analyte in the pathway.
    edgeList <- expand.grid(analytesInPathway, analytesInPathway)
    edgeListUnique <- edgeList[which(edgeList[,1] != edgeList[,2]),]
    igraph::add.edges(G, edgeListUnique)
  }
  
  # Add edges to the graph and add pathways to the edge dict.
  for(edge in igraph::as_ids(igraph::E(G))){
    
    # Get the source and target nodes.
    source <- strsplit(edge, "|")[[1]][1]
    target <- strsplit(edge, "|")[[1]][2]
    
    # Find the pathways that include the source, target, and both.
    sourceAllPathways <- analyteHasPathway[which(analyteHasPathway$Analyte == source), "Pathway"]
    targetAllPathways <- analyteHasPathway[which(analyteHasPathway$Analyte == target), "Pathway"] 
    sourceTargetSharedPathways <- intersect(sourceAllPathways, targetAllPathways) 
    
    # Add all shared pathways to the dictionary.
    edgePathwayDict[edge] <- sourceTargetSharedPathways
    
    # Compute the edge weight.
    sourceSum <- sum(unlist(lapply(sourceAllPathways, function(pathway){
      return(pathwayWeight[pathway])
    }))
    targetSum <- sum(unlist(lapply(targetAllPathways, function(pathway){
      return(pathwayWeight[pathway])
    }))
    sharedSum <- sum(unlist(lapply(sourceTargetSharedPathways, function(pathway){
      return(pathwayWeight[pathway])
    }))
    edgeWeight <- sharedSum / (sourceSum + targetSum - sharedSum)
    
    # Modify the edge weight in the graph.
    igraph::add.edges(G, edge, attr = c("weight", edgeWeight))
  }
  
  # Create the knowledge graph.
  adj <- igraph::graph.adjacency(G, mode="undirected", weighted=TRUE)
  knowledgeGraph <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                 adjacencyMatrix = adj, 
                                 edgeDict = edgePathwayDict,
                                 clusters = list(),
                                 analytePathwayMapping = analyteHasPathway,
                                 clusterPathwayDict = list())
  
  # Save the file.
  if(pathwayGraphFile != NULL){
    saveRDS(knowledgeGraph, pathwayGraphFile)
  }
}

#' Finds clusters in the pathway knowledge graph using Louvain clustering.
#' @param knowledgeGraph An object of type TARDIGRADE_KnowledgeGraph.
#' @param clusterCount The number of clusters.
#' @returns An object of type TARDIGRADE_KnowledgeGraph.
#' @export
findClusters <- function(knowledgeGraph = NULL, clusterCount = 0){
  
  # Do clustering.
  G <- igraph::graph_from_adjacency_matrix(knowledgeGraph$adjacencyMatrix)
  clusters <- igraph::cluster_louvain(G)
  
  # Map clusters to pathways.
  pathways <- lapply(clusters, function(cluster){
    
    # Subset only edges where both nodes are in the cluster.
    clusterEdges <- E(G)[intersect(which(ends(G, E(G))[,1] %in% cluster),
                                   which(ends(G, E(G))[,2] %in% cluster))]
    
    # Find the pathways associated with each edge.
    pathwaysForCluster <- unique(unlist(lapply(clusterEdges, function(edge){
      return(knowledgeGraph$edgeDict[which(knowledgeGraph$edgeDict[,1] == edge)], 2)
    })))
    
    # Return local list of pathways.
    return(pathwaysForCluster)
  })
  names(pathways) <- names(clusters)
  
  # Return the new knowledge graph with the clusters.
  return(methods::new("TARDIGRADE_KnowledgeGraph", 
                      adjacencyMatrix = knowledgeGraph$adjacencyMatrix, 
                      edgeDict = knowledgeGraph$edgeDict,
                      clusters = clusters,
                      analytePathwayMapping = knowledgeGraph$analytePathwayMapping,
                      clusterPathwayDict = clusterToPathway))
}

#' Modifies the edges in the graph given a data-driven network, with the following rules:
#' - If an edge exists in both the knowledge graph and the data-driven network, increase
#'   the edge weight to 1 and make no other changes.
#' - If an edge exists in the knowledge graph and CANNOT be present in the data-driven
#'   network (i.e. the nodes do not exist or are not of the correct analyte types to be
#'   connected), make no changes.
#' - If an edge exists in the knowledge graph and COULD be present in the data-driven network
#'   (but isn't), remove it.
#' - If an edge does not exist in the knowledge graph but does exist in the data-driven
#'   network, add it with no pathways to the dict.
#' @param knowledgeGraph An object of type TARDIGRADE_KnowledgeGraph.
#' @param network A data frame representing binary edges in a network.
#' @param sourceAnalytes The list of source analytes (e.g., TFs) evaluated.
#' @param targetAnalytes The list of target analytes (e.g., genes) evaluated.
#' @returns An object of type TARDIGRADE_KnowledgeGraph.
#' @export
modifyGraph <- function(knowledgeGraph = NULL, network = NULL, sourceAnalytes = NULL,
                        targetAnalytes = NULL){
  
  # Get all possible edges connecting source to target analytes.
  sourceToTarget <- rbind(lapply(sourceAnalytes, function(source){
    return(rbind(lapply(targetAnalytes, function(target){
      return(data.frame(source = source, target = target))
    })))
  }))
  
  # Get all possible edges connecting source to source analytes.
  sourceToSource <- rbind(lapply(sourceAnalytes, function(source){
    return(rbind(lapply(setdiff(sourceAnalytes, source), function(target){
      return(data.frame(source = source, target = target))
    })))
  }))
  
  # Get all possible edges connecting target to target analytes.
  targetToTarget <- rbind(lapply(targetAnalytes, function(source){
    return(rbind(lapply(setdiff(targetAnalytes, source), function(target){
      return(data.frame(source = source, target = target))
    })))
  }))
  
  # Get all possible edges connecting source to non-existent.
  unmeasured <- setdiff(c(colnames(knowledgeGraph$adjacencyMatrix),
                          rownames(knowledgeGraph$adjacencyMatrix)),
                        c(sourceAnalytes, targetAnalytes))
  sourceToUnmeasured <- rbind(lapply(sourceAnalytes, function(source){
    return(rbind(lapply(unmeasured, function(target){
      return(data.frame(source = source, target = target))
    })))
  }))
  targetToUnmeasured <- rbind(lapply(targetAnalytes, function(source){
    return(rbind(lapply(unmeasured, function(target){
      return(data.frame(source = source, target = target))
    })))
  }))
  
  # Add the weight column to each.
  network$weight <- 1
  expandedNetwork$weight <- 0
  sourceToSource$weight <- -1
  targetToTarget$weight <- -1
  sourceToUnmeasured$weight <- -1
  targetToUnmeasured$weight <- -1
  
  # Expand the adjacency list.
  networkEdges <- paste(network[,1], network[,2], sep = "|")
  sourceToTargetEdges <- paste(sourceToTarget[,1], sourceToTarget[,2], sep = "|")
  sourceToSourceEdges <- paste(sourceToSource[,1], sourceToSource[,2], sep = "|")
  targetToTargetEdges <- paste(targetToTarget[,1], targetToTarget[,2], sep = "|")
  expandedNetwork <- rbind(list(network, 
                                sourceToTarget[which(sourceToTargetEdges %in% setdiff(sourceToTargetEdges, networkEdges))],
                                sourceToSource[which(sourceToSourceEdges %in% setdiff(sourceToSourceEdges, sourceToTargetEdges))],
                                targetToTarget[which(targetToTargetEdges %in% setdiff(targetToTargetEdges, sourceToTargetEdges))],
                                sourceToUnmeasured, targetToUnmeasured))
  
  # Build an adjacency matrix using the input network.
  networkIgraph <- igraph::graph_from_data_frame(expandedNetwork)
  networkAdj <- igraph::graph.adjacency(networkIgraph, mode="undirected", weighted = "TRUE")
  networkAdj[rownames(knowledgeGraph$adjacencyMatrix), colnames(knowledgeGraph$adjacencyMatrix)]
  
  # Modify the knowledge graph matrix as follows:
  # Where the weight in the network matrix is higher (1), set to the network matrix value.
  # Where the weight in the knowledge graph is higher (-1 in network matrix), set to the knowledge graph value.
  # Finally, where the weight in the network matrix is 0, set the knowledge graph matrix to 0.
  modifiedMatrix <- pmax(knowledgeGraph, networkAdj)
  modifiedMatrix[which(networkAdj == 0)] <- 0
  
  # Create and return the new graph.
  knowledgeGraph <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                 adjacencyMatrix = modifiedMatrix, 
                                 edgeDict = edgePathwayDict,
                                 clusters = list())
  return(knowledgeGraph)
}

#' Performs heat diffusion for a single network, given the knowledge graph.
#' This happens by multiplying the seed values in a network by the conductance
#' matrix, which is the inverse of a negative normalized graph Laplacian.
#' @param knowledgeGraphs An list of TARDIGRADE_KnowledgeGraphs, one modified
#' graph for each sample.
#' @param scaledAnalyteMeasurements A list of analyte measurement matrices (e.g. DNA methylation,
#' gene expression, metabolite abundance, protein abundance), with each analyte
#' scaled between 0 and 1.
#' @export
diffuseHeat <- function(knowledgeGraphs = list(), scaledAnalyteMeasurements = list()){
  
  # Create normalized graph Laplacian matrices.
  normLaplacians <- lapply(knowledgeGraphs, function(knowledgeGraph){
    return(getNormalizedLaplacian(knowledgeGraph))
  })
  
  # Create the conductance matrix.
  conductance <- lapply(normLaplacians, function(normLaplacian){
    return(MASS::pinv(-1 * normLaplacian))
  })
  
  # Diffuse heat.
  heat <- lapply(1:length(conductance), function(i){
    return(conductance[i] %*% scaledAnalyteMeasurements[i])
  })
  
  # Return heat conductance.
  return(heat)
}

#' Obtains the normalized Laplacian for a knowledge graph. This is a helper
#' function for diffuseHeat.
#' @param knowledgeGraphs An TARDIGRADE_KnowledgeGraph.
#' @returns The normalized Laplacian as a matrix.
getNormalizedLaplacian <- function(knowledgeGraph){
  
  # Extract adjacency matrix.
  A <- knowledgeGraph$adjacencyMatrix
  
  # Compute the degree matrix.
  D <- matrix(data = rep(0, length(A)), nrow = nrow(A))
  diag(D) <- psum(D)
  
  # Compute the Laplacian matrix.
  L <- D - A
  
  # Normalize the Laplacian matrix.
  reciprocalD <- MASS::pinv(D)
  reciprocalSquareRootD <- reciprocalD ^ (1/2)
  normL <- reciprocalSquareRootD %*% L %*% reciprocalSquareRootD
  
  # Return normalized Laplacian.
  return(normL)
}

#' Performs pathway cluster enrichment analysis using a GSEA-like procedure.
#' Specifically, given the heat values of each node for each sample for 
#' both original and permuted data, (1) Find differential heat score for
#' each node and then rank them, (2) Find the composite ranking score for each
#' pathway cluster, and (3) Find the significance of that score.
#' @param pvalue FDR-adjusted p-value cutoff for significance.
#' @param knowledgeGraphsControlsHeat A list of heat diffusion results for all
#' samples in the control group.
#' @param knowledgeGraphsCasesHeat A list of heat diffusion results for all samples
#' in the case group.
#' @param permCount Number of permutations.
#' @returns A data frame containing p-values for each pathway cluster.
#' @export
doPathwayAnalysis <- function(knowledgeGraphsControlsHeat = list(), 
                              knowledgeGraphsCasesHeat = list(), pvalue = 1,
                              permCount = 1000){
  
  # Compute point biserial correlations between heat scores and case/control
  # status.
  geneCaseCorrelations <- unlist(lapply(rownames(knowledgeGraphsControlsHeat), function(node){
    
    # Get the controls for this node.
    controlsForNode <- unlist(lapply(knowledgeGraphsControlsHeat, function(graph){
      return(graph[node])
    }))
    return(controlsForNode)

    # Get the cases for this node.
    casesForNode <- unlist(lapply(knowledgeGraphsCasesHeat, function(graph){
      return(graph[node])
    }))
    return(casesForNode)
    
    # Return the point biserial correlation.
    caseControlVector <- c(controlsForNode, casesForNode)
    caseControlIDVector <- c(rep(0, length(controlsForNode)),
                             rep(0, length(casesForNode)))
    return(stats::cor.test(caseControlVector, caseControlIDVector)$cor)
  })
  names(geneCaseCorrelations) <- rownames(knowledgeGraphsControlsHeat)
  
  # Rank the nodes by their correlations.
  correlationOrder <- [order(-1 * geneCaseCorrelations)]
  nodeCorrelationNames <- names(geneCaseCorrelations)[correlationOrder]
  nodeCorrelationRank <- geneCaseCorrelations[correlationOrder]
  names(nodeCorrelationRank) <- nodeCorrelationNames
  
  # Obtain the pathway cluster information needed to compute normalized enrichment.
  uniquePathwayClusters <- sort(unique(pathwayClusters[,2]))
  pathwayClusterSizes <- unlist(lapply(uniquePathwayClusters, function(cluster){
    return(nrow(knowledgeGraphsControlsHeat@clusters[which(knowledgeGraphsControlsHeat@clusters[,2] == cluster)]))
  }))
  
  # Get the normalized enrichment scores for each pathway.
  normEnrichmentScore <- getNormalizedEnrichmentScore(correlationRanking = nodeCorrelationRank,
                                                      pathwayClusterMapping = knowledgeGraphsControlsHeat@clusters,
                                                      pathwayClusterNames = uniquePathwayClusters,
                                                      pathwayClusterSizes = pathwayClusterSizes)
  
  # Get the normalized enrichment scores for permuted pathway assignments.
  set.seed(1)
  permutedNormEnrichmentScoreList <- lapply(1:permCount, function(permIdx){
    
    # Randomly re-order the pathway assignment.
    permPathwayClusterMapping <- sample(knowledgeGraphsControlsHeat@clusters[,1], replace = FALSE)
    permNormEnrichmentScore <- getNormalizedEnrichmentScore(correlationRanking = nodeCorrelationRank,
                                                        pathwayClusterMapping = permPathwayClusterMapping,
                                                        pathwayClusterNames = uniquePathwayClusters,
                                                        pathwayClusterSizes = pathwayClusterSizes)
  })
  
  # Reorder the normalized enrichment scores.
  permutedNormEnrichmentScores <- lapply(names(permNormEnrichmentScoreList[1]), function(pathway){
    return(unlist(lapply(1:permCount, function(permIdx){
      return(permutedNormEnrichmentScoreList[permIdx][pathway])
    })))
  }))
  names(permutedNormEnrichmentScores) <- names(permNormEnrichmentScoreList[1])
  
  # Calculate significance.
  pvalEnrichmentScores <- lapply(normEnrichmentScore, function(score){
    return(t.test(x = permutedNormEnrichmentScores, mu = score, alternative = "less"))
  })
  names(pvalEnrichmentScores) <- names(normEnrichmentScore)
  
  # Return significance.
  return(pvalEnrichmentScores)
}

#' Obtains the normalized enrichment score for each pathway cluster based on
#' the node correlation rankings and the cluster definitions.
#' @param nodeCorrelationRank A ranked vector of nodes (analytes) by correlation
#' with case/control status.
#' @param pathwayClusterMapping The mapping from gene to pathway cluster.
#' @param pathwayClusterNames A list of names for each pathway cluster.
#' @param pathwayClusterSizes The sizes of each pathway cluster.
#' @returns A normalized enrichment score for each pathway cluster.
getNormalizedEnrichmentScore <- function(nodeCorrelationRank = NULL, pathwayClusterMapping = NULL,
                                         pathwayClusterNames = NULL, pathwayClusterSizes = NULL){
  
  # Initialize the running sum for each pathway cluster.
  pathwayClusterSums <- rep(0, length(pathwayClusterNames))
  pathwayClusterMaxes <- rep(0, length(pathwayClusterMaxes))
  names(pathwayClusterSums) <- pathwayClusterNames
  names(pathwayClusterMaxes) <- pathwayClusterNames
  
  # Walk through the list of nodes and add to the sum. Formula is taken
  # from the GSEA appendix.
  for(i in 1:length(nodeCorrelationRank)){
    
    # Extract the rank.
    rank <- nodeCorrelationRank[i]
    
    # Map the gene to a pathway.
    pathway <- pathwayClusterMapping[names(nodeCorrelationRank)[i], 2]
    
    sizeFactorSubtraction <- sqrt(pathwaySizes[pathway] / (length(nodeCorrelationRank) - pathwaySizes[pathway]))
    
    # Add to the pathway to which the gene belongs. Use the addition factor and the correlation.
    pathwayClusterSums[pathway] <- pathwayClusterSums[pathway] + 
      rank * sqrt((length(nodeCorrelationRank) - pathwaySizes[pathway]) / pathwaySizes[pathway])
    
    # Subtract from the pathway to which the gene does not belong. Use the subtraction factor
    # and the correlation.
    otherPathways <- setdiff(pathwayClusterNames, pathway)
    pathwayClusterSums[otherPathways] <- pathwayClusterSums[otherPathways] - 
      rank * sqrt(pathwaySizes[otherPathways] / (length(nodeCorrelationRank) - pathwaySizes[otherPathways]))
    
    # Update the max, if applicable.
    if(pathwayClusterSums[pathway] > pathwayClusterMaxes[pathway]){
      pathwayClusterMaxes[pathway] <- pathwayClusterSums[pathway]
    }
  }
  
  # Return the maxes, which are the normalized enrichment scores.
  return(pathwayClusterMaxes)
}
