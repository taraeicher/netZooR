#'  TARDIGRADE_KnowledgeGraph class
#'  This object type represents the list of results in a given directory.
#'  
#'  @name TARDIGRADE_KnowledgeGraph-class
#'  @rdname TARDIGRADE_KnowledgeGraph-class
#'  @exportClass TARDIGRADE_KnowledgeGraph
#'  @slot adjacencyMatrix The adjacency matrix for the graph.
#'  @slot edgeDict The dictionary mapping edges to pathways. This is a list.
#'  @slot clusters A named list of clusters.
#'  @slot analytePathwayMapping The original data frame mapping analytes to pathways.
#'  @slot clusterPathwayDict The dictionary mapping clusters to pathways.
methods::setClass(
  Class="TARDIGRADE_KnowledgeGraph",
  representation(adjacencyMatrix = "dgCMatrix",
                 edgeDict = "list",
                 clusters = "list",
                 analytePathwayMapping = "data.frame",
                 clusterPathwayDict = "list")
)

#' Implements the setup for RaMP.
#' @param version The version to use for setup.
#' @returns A data frame mapping from Analyte (rampId:common_name) to Pathway
#' (type:pathwayName)
#' @export
DoRaMPSetup <- function(version = NULL){
  
  # Set the query.
  query_new_ramp <- "
    select
         source.rampId,
         group_concat(distinct source.commonName) as common_name,
         pathway.*
    from
         analytehaspathway,
         pathway,
         source
    where
         analytehaspathway.pathwayRampId = pathway.pathwayRampId
         and analytehaspathway.rampId = source.rampId
    group by source.rampId, pathway.pathwayRampId"
  
  # Get the latest version of the database.
  newRampDB <- RaMP()
  # newRampDB <- NULL
  # if(!is.null(version)){
  #   versionsAvailable <- listAvailableRaMPDbVersions()
  #   str(versionsAvailable)
  #   newRampDB <- RaMP(version=version)
  #   stop(paste("ERROR:", version, "is not available from RaMP."))
  # }else{
  #   newRampDB <- RaMP()
  # }
  
  # Run the query.
  res <- runQuery(sql = query_new_ramp, db = newRampDB)
  # Should we remove the pathwayCategory = "smpdb3"?
  
  # Return the mapping.
  analyteHasPathway <- data.frame(Analyte = paste(res$rampId, res$common_name, sep = ":"),
                                  Pathway = paste(res$type, res$pathwayName, sep = ":"))

  # Return the DB object.
  return(analyteHasPathway)
}

#' Builds the pathway graph from RaMP and stores in in the file specified by the user.
#'
#' @param analyteHasPathway A data frame mapping analytes (column "Analyte") to 
#' pathways (column "Pathway").
#' @param pathwaySizeLimit The maximum number of analytes a pathway can contain.
#' Any pathways that are larger will be excluded from analysis. If -1, all are included.
#' @returns An object of type TARDIGRADE_KnowledgeGraph.
#' @export
BuildPathwayGraph <- function(analyteHasPathway, pathwaySizeLimit = -1){
  
  # Check that the input is correct.
  if(!is(analyteHasPathway, "data.frame") || length(colnames(analyteHasPathway)) != 2
     || colnames(analyteHasPathway)[1] != "Analyte" || colnames(analyteHasPathway)[2] != "Pathway"
     || !is.character(analyteHasPathway[,1]) || !is.character(analyteHasPathway[,2])){
    stop("ERROR: analyteHasPathway must be a data frame with character columns 'Analyte' and 'Pathway'")
  }
  if(!is.numeric(pathwaySizeLimit)){
    stop("ERROR: pathwaySizeLimit must be numeric")
  }
  
  # Filter input to include the size limit.
  pathwaySizes <- table(analyteHasPathway$Pathway)
  if(pathwaySizeLimit > -1){
    pathwaysWithinLimit <- names(pathwaySizes)[which(pathwaySizes < pathwaySizeLimit)]
    analyteHasPathway <- analyteHasPathway[which(analyteHasPathway$Pathway %in% pathwaysWithinLimit),]
    message(paste("Retained", length(pathwaysWithinLimit), "out of", length(pathwaySizes), "pathways"))
  }
  
  # Build the graph and the matching dictionary.
  # Initialize the graph and dict.
  G <- igraph::make_empty_graph(directed = FALSE)
  G <- set_edge_attr(G, "weight", value = numeric(0))
  edgePathwayDict <- list()
  
  # Add the vertices (analytes).
  uniqueAnalytes <- unique(analyteHasPathway$Analyte)
  G <- igraph::add.vertices(G, length(uniqueAnalytes), name = uniqueAnalytes)

  # Obtain weights for each pathway.
  pathwayWeight <- list()
  pathwaySizesUpdated <- table(analyteHasPathway$Pathway)
  uniquePathways <- names(pathwaySizesUpdated[order(pathwaySizesUpdated)])
  pathwayToEdge <- data.frame(Pathway = c(), Edge = c())
  for(i in 1:length(uniquePathways)){
    
    # Get all analytes in the pathway.
    pathway <- uniquePathways[i]
    analytesInPathway  <- unique(analyteHasPathway[which(analyteHasPathway$Pathway == pathway), "Analyte"])
    
    # Set the pathway's weight to be the probability of a given analyte in the pathway
    # directly influencing any other analyte in the pathway.
    pathwayWeight[pathway] <- 1.0 / (length(analytesInPathway) - 1)
    
    # Add an edge for each analyte in the pathway.
    clique <- make_full_graph(length(analytesInPathway), directed = FALSE)
    V(clique)$name <- analytesInPathway
    edgeMat <- t(igraph::as_edgelist(clique))
    edgeVec <- as.vector(edgeMat)
    
    
    # Add the edge list if it is greater than 0.
    if(length(edgeVec) > 0){

      # Add the edge to the pathway dict.
      edgeVecFormatted <- unique(paste(t(edgeMat)[,1], t(edgeMat)[,2], sep = "|"))
      pathwayToEdge <- rbind(pathwayToEdge, data.frame(Pathway = rep(pathway, length(edgeVec)),
                                                       Edge = edgeVecFormatted))
      
      # Add the edges to the graph, where the weight is the maximum pathway weight,
      # and the pathway weight is 1 / (pathway size - 1).
      isInGraph <- which(edgeVecFormatted %in% igraph::as_ids(igraph::E(G)))
      edgesInGraph <- edgeVecFormatted[isInGraph]
      edgesNotInGraphMat <- edgeMat[,setdiff(1:ncol(edgeMat), isInGraph)]
      edgesNotInGraph <- as.vector(edgesNotInGraphMat)
      edgesNotInGraphFormatted <- edgeVecFormatted[setdiff(1:ncol(edgeMat), isInGraph)]
      pathwayWeight <- 1 / (length(analytesInPathway) - 1)
      # Add the edges not already in the graph and set their weights.
      G <- add_edges(G, edgesNotInGraph)
      E(G)[edgesNotInGraphFormatted]$weight <- pathwayWeight
      # For edges already in the graph, modify their weights.
      weights <- pmax(E(G)[edgesInGraph]$weight, pathwayWeight)
      E(G)[edgesInGraph]$weight <- weights
    }

    # Add all shared pathways to the dictionary.
    message(paste("Added edges for", i, "out of", length(uniquePathways), "pathways"))
  }
  G <- igraph::simplify(G, remove.multiple = TRUE)

  # Add edges to the graph pathways to the edge dict.
  uniqueEdges <- unique(pathwayToEdge$Edge)
  if(length(uniqueEdges) > 0){
    edgePathwayDict <- as.list(do.call(c, lapply(1:length(uniqueEdges), function(i){
      edge <- uniqueEdges[i]
      retval <- paste(unique(pathwayToEdge[which(pathwayToEdge$Edge == edge), "Pathway"]), collapse = ";")
      message(paste("Built pathway dictionary for", i, "out of", length(uniqueEdges), "edges"))
      return(retval)
    })))
    names(edgePathwayDict) <- uniqueEdges
  }
  
  # Create the knowledge graph.
  adj <- igraph::as_adjacency_matrix(G, attr = "weight", sparse = TRUE)
  knowledgeGraph <- methods::new("TARDIGRADE_KnowledgeGraph", 
                                 adjacencyMatrix = adj, 
                                 edgeDict = edgePathwayDict,
                                 clusters = list(),
                                 analytePathwayMapping = analyteHasPathway,
                                 clusterPathwayDict = list())
  
  # Return
  return(knowledgeGraph)
}

#' Calculates the Fruchterman-Reingold layout.
#'
#' @param knowledgeGraph An object of type TARDIGRADE_KnowledgeGraph.
#' @return A numeric matrix with two (dim=2) or three (dim=3) columns, 
#' and as many rows as the number of vertices, the x, y and potentially z 
#' coordinates of the vertices.
CalculateFruchtermanReingoldLayout <- function(knowledgeGraph){
  
  # Convert the knowledge graph to a data frame.
  kgEdges <- which(knowledgeGraph@adjacencyMatrix != 0, arr.ind = TRUE)
  weights <- knowledgeGraph@adjacencyMatrix[kgEdges]
  kgEdgesNamed <- data.frame(
    from = rownames(knowledgeGraph@adjacencyMatrix)[kgEdges[, 1]],
    to = colnames(knowledgeGraph@adjacencyMatrix)[kgEdges[, 2]],
    weight = weights,
    stringsAsFactors = FALSE
  )
  str(kgEdgesNamed)
  
  # Invert the edge weights. Add a small "fudge factor" to prevent weights from
  # going to 0.
  kgEdgesNamed$weight <- max(kgEdgesNamed$weight) + 0.001 - kgEdgesNamed$weight
  str(kgEdgesNamed)

  # Construct nodes.
  uniqueNodes <- unique(c(kgEdgesNamed$from, kgEdgesNamed$to))
  str(uniqueNodes)
  nodeAttrs <- data.frame(node = uniqueNodes)
  str(nodeAttrs)
  
  # Create a graph object.
  graph <- igraph::graph_from_data_frame(kgEdgesNamed, vertices = nodeAttrs, directed = FALSE)
  str(graph)
  
  # Calculate and return the layout.
  message("Calculating Fruchterman-Reingold layout...")
  fruchtermanReingold <- igraph::layout_with_fr(graph)
  rownames(fruchtermanReingold) <- V(graph)$name
  message("Fruchterman-Reingold layout complete!")
  return(fruchtermanReingold)
}

#' Plots the knowledge graph using a Kamada-Kawai spatial layout. Edges are not shown,
#' as the location represents the orientation of the node in the graph and the
#' knowledge graph is expected to be dense. Nodes are single-pixel with low alpha
#' values by default, to show overlap. Clusters are color-coded if present. Otherwise,
#' no color is shown.
#'
#' @param knowledgeGraph An object of type TARDIGRADE_KnowledgeGraph.
#' @param nodeSize The size of the node in pixels.
#' @param alpha The opacity of each node.
#' @param edgeWidth The width of the edges.
PlotKnowledgeGraph <- function(knowledgeGraph, nodeSize = 1, alpha = 0.1, edgeWidth = 0,
                               layout = NULL){
  
  # # Convert the knowledge graph to a data frame.
  # kgEdges <- which(knowledgeGraph@adjacencyMatrix != 0, arr.ind = TRUE)
  # weights <- knowledgeGraph@adjacencyMatrix[kgEdges]
  # kgEdgesNamed <- data.frame(
  #   from = rownames(knowledgeGraph@adjacencyMatrix)[kgEdges[, 1]],
  #   to = colnames(knowledgeGraph@adjacencyMatrix)[kgEdges[, 2]],
  #   weight = weights,
  #   stringsAsFactors = FALSE
  # )
  # 
  # # Set the node attributes.
  # vertexLabelSize = 0
  # vertexLabelOffset = 0
  # vertexLabelColor = "black"
  # uniqueNodes <- unique(c(kgEdgesNamed$from, kgEdgesNamed$to))
  # nodeAttrs <- data.frame(node = uniqueNodes,
  #                         color = rep(rgb(0, 0, 0, max = 255, alpha = alpha * 255), length(uniqueNodes)),
  #                         size = rep(nodeSize, length(uniqueNodes)),
  #                         frame.width = rep(0, length(uniqueNodes)),
  #                         label.color = vertexLabelColor, label.cex = vertexLabelSize,
  #                         label.dist = vertexLabelOffset)
  # rownames(nodeAttrs) <- uniqueNodes
  # 
  # # Add edge attributes.
  # kgEdgesNamed$width <- edgeWidth
  # 
  # # Create a graph object.
  # graph <- igraph::graph_from_data_frame(kgEdgesNamed, vertices = nodeAttrs, directed = FALSE)
  # str(graph)
  # 
  # # Plot.
  # igraph::plot.igraph(graph, layout = layout, vertex.label = labels)
  plot(x = layout[,1], y = layout[,2], pch = 20, 
       col = rgb(red = 0, blue = 0, green = 0, alpha = 0.1))
}

#' Finds clusters in the pathway knowledge graph using Louvain clustering.
#' @param knowledgeGraph An object of type TARDIGRADE_KnowledgeGraph.
#' @param clusterCount The number of clusters.
#' @returns An object of type TARDIGRADE_KnowledgeGraph.
#' @export
FindClusters <- function(knowledgeGraph, clusterCount){
  
  # Check input.
  if(length(knowledgeGraph@adjacencyMatrix) == 0 || max(knowledgeGraph@adjacencyMatrix, na.rm = TRUE) == 0){
    stop("ERROR: The knowledge graph has no edges.")
  }
  if(nrow(knowledgeGraph@adjacencyMatrix) < clusterCount){
    stop("ERROR: The number of clusters must be smaller than the number of nodes.")
  }
  
  # Do clustering.
  G <- igraph::graph_from_adjacency_matrix(knowledgeGraph@adjacencyMatrix, 
                                           mode = "undirected", weighted = TRUE)
  clusters <- igraph::cluster_louvain(G)
  uniqueClusters <- sort(unique(clusters$membership))
  formattedClusters <- lapply(uniqueClusters, function(cluster){
    return(clusters$names[which(clusters$membership == cluster)])
  })
  names(formattedClusters) <- uniqueClusters

  # Map clusters to pathways.
  edges <- igraph::as_ids(igraph::E(G))
  pathways <- lapply(formattedClusters, function(cluster){
    
    # Subset only edges where both nodes are in the cluster.
    clusterEdges <- edges[intersect(which(igraph::ends(G, edges)[,1] %in% cluster),
                                   which(igraph::ends(G, edges)[,2] %in% cluster))]
    
    # Find the pathways associated with each edge.
    pathwaysForCluster <- unique(unlist(lapply(clusterEdges, function(edge){
      return(knowledgeGraph@edgeDict[edge])
    })))

    # Return local list of pathways.
    return(pathwaysForCluster)
  })
  names(pathways) <- names(formattedClusters)

  # Return the new knowledge graph with the clusters.
  return(methods::new("TARDIGRADE_KnowledgeGraph", 
                      adjacencyMatrix = knowledgeGraph@adjacencyMatrix, 
                      edgeDict = knowledgeGraph@edgeDict,
                      clusters = formattedClusters,
                      analytePathwayMapping = knowledgeGraph@analytePathwayMapping,
                      clusterPathwayDict = pathways))
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
#' @returns An adjacency matrix.
#' @export
ModifyGraph <- function(knowledgeGraph = NULL, network = NULL, sourceAnalytes = NULL,
                        targetAnalytes = NULL){
  
  # Verify that the input types are valid.
  if(!is(knowledgeGraph, "TARDIGRADE_KnowledgeGraph") || !is(network, "data.frame")
     || !is.character(sourceAnalytes) || !is.character(targetAnalytes)){
    stop("ERROR: Input formats must be: TARDIGRADE_KnowledgeGraph object (knowledgeGraph), data frame (network), and character (sourceAnalytes and targetAnalytes)")
  }
  
  # Verify that there is a match in the analytes.
  if(nrow(network) == 0){
    stop("ERROR: The input network is empty")
  }
  if(length(setdiff(knowledgeGraph@analytePathwayMapping$Analyte,
                    c(network[,1], network[,2]))) == length(unique(knowledgeGraph@analytePathwayMapping$Analyte))){
    stop("ERROR: The input network and knowledge graph do not overlap")
  }
  if(length(setdiff(c(network[,1], network[,2]),
                    sourceAnalytes)) == length(unique(c(network[,1], network[,2])))){
    stop("ERROR: Source or target analytes do not match the input network")
  }
  if(length(setdiff(c(network[,1], network[,2]),
                    targetAnalytes)) == length(unique(c(network[,1], network[,2])))){
    stop("ERROR: Source or target analytes do not match the input network")
  }
  
  # Expand knowledge graph to include all nodes exclusive to the source and target
  # lists (using adjacency matrix).
  adjacencyMatrixExpanded <- knowledgeGraph@adjacencyMatrix
  oldNames <- rownames(adjacencyMatrixExpanded)
  analytesNotInKG <- setdiff(c(sourceAnalytes, targetAnalytes), rownames(adjacencyMatrixExpanded))
  adjacencyMatrixExpanded <- cbind(adjacencyMatrixExpanded, matrix(rep(0, nrow(adjacencyMatrixExpanded) *
                                                                         length(analytesNotInKG)),
                                                                   nrow = nrow(adjacencyMatrixExpanded),
                                                                   ncol = length(analytesNotInKG)))
  adjacencyMatrixExpanded <- rbind(adjacencyMatrixExpanded, matrix(rep(0, length(analytesNotInKG) * 
                                                                         ncol(adjacencyMatrixExpanded)),
                                                                   nrow = length(analytesNotInKG),
                                                                   ncol = ncol(adjacencyMatrixExpanded)))
  rownames(adjacencyMatrixExpanded) <- c(oldNames, analytesNotInKG)
  colnames(adjacencyMatrixExpanded) <- rownames(adjacencyMatrixExpanded)
  
  # Add all analytes in the source-target list but not in the network with weights of 0.
  networkAdj <- igraph::as_adjacency_matrix(igraph::graph_from_data_frame(network, directed = FALSE))
  oldNames <- rownames(networkAdj)
  notInNetwork <- setdiff(c(sourceAnalytes, targetAnalytes), rownames(networkAdj))
  networkAdjExpanded <- cbind(networkAdj, matrix(rep(0, nrow(networkAdj) *
                                                               length(notInNetwork)), 
                                                         nrow = nrow(networkAdj),
                                                         ncol = length(notInNetwork)))
  networkAdjExpanded <- rbind(networkAdjExpanded, matrix(rep(0, ncol(networkAdjExpanded) *
                                                               length(notInNetwork)),
                                                         ncol = ncol(networkAdjExpanded),
                                                         nrow = length(notInNetwork)))
  rownames(networkAdjExpanded) <- c(oldNames, notInNetwork)
  colnames(networkAdjExpanded) <- rownames(networkAdjExpanded)
  
  # Add all analytes exclusive to the knowledge graph to the network with weights of -1.
  oldNames <- rownames(networkAdjExpanded)
  notInSrcTarget <- setdiff(rownames(adjacencyMatrixExpanded), c(sourceAnalytes, targetAnalytes))
  networkAdjExpanded <- cbind(networkAdjExpanded, matrix(-1, nrow = nrow(networkAdjExpanded), 
                                                 ncol = length(notInSrcTarget)))
  networkAdjExpanded <- rbind(networkAdjExpanded, matrix(-1, ncol = ncol(networkAdjExpanded), 
                                                 nrow = length(notInSrcTarget)))
  rownames(networkAdjExpanded) <- c(oldNames, notInSrcTarget)
  colnames(networkAdjExpanded) <- rownames(networkAdjExpanded)
  
  # Rearrange the network nodes.
  networkAdjRearranged <- networkAdjExpanded[rownames(adjacencyMatrixExpanded), 
                                             colnames(adjacencyMatrixExpanded)]

  # Zero out the places in the knowledge graph where the network is 0. These represent
  # connections found not to be present in the sample.
  adjacencyMatrixExpanded[which(networkAdjRearranged == 0)] <- 0

  # For everything else, take the maximum, which is either the knowledge graph
  # value for edges exclusive to the knowledge graph, or the network value for
  # edges in the network.
  modifiedAdjMat <- adjacencyMatrixExpanded
  modifiedAdjMat <- pmax(adjacencyMatrixExpanded, networkAdjRearranged)

  # Return.
  return(modifiedAdjMat)
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
  }))
  names(geneCaseCorrelations) <- rownames(knowledgeGraphsControlsHeat)
  
  # Rank the nodes by their correlations.
  correlationOrder <- order(-1 * geneCaseCorrelations)
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
  })
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
