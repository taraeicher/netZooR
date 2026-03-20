#' Given a set of genes of interest, full unipartite networks with scores (one network for each sample), a significance
#' cutoff for statistical testing, and a hop constraint, UNAGI finds a subnetwork of
#' significant edges connecting the genes.
#' @param nodeSet A character vector of genes comprising the targets of interest.
#' @param network A network with the first column ("source") consisting of source nodes,
#' the second ("target") consisting of target nodes, and the third ("score") consisting of scores
#' @param hopConstraint The maximum number of hops to be considered between gene pairs.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' Default is FALSE.
#' @returns A unipartite subnetwork in the same format as the original networks.
#' @export
RunUNAGI <- function(nodeSet, network, hopConstraint,
                     verbose = FALSE, topX=NULL) {
  #this entire body was edited Jul 13
  if (!is.character(nodeSet) || !is.data.frame(network) || !is.numeric(hopConstraint))
    stop("Wrong input type! nodeSet must be a character vector. network must be a data frame.",
         " hopConstraint must be a scalar numeric value.")
  if (!all(c("source", "target", "score") %in% colnames(network)))
    stop("Network must have source genes in the first column, ",
         "target genes in the second column, and scores in the third column.")
  if(hopConstraint < 1){
    stop("hopConstraint must be at least 1.")
  }
  if(length(nodeSet) < 2){
    stop("Node set must contain at least 2 nodes.")
  }
  
  sigEdges <- network[, c("source", "target")]
  rownames(sigEdges) <- paste(sigEdges$source, sigEdges$target, sep = "__")
  
  ##Build subnetworks exactly like BLOBFISH 
  significantSubnetworks <- FindEdgesForHopU(geneSet = nodeSet,
                                             network = network,
                                             hopConstraint = ceiling(hopConstraint / 2),
                                             verbose = verbose, topX = topX)
  odd <- hopConstraint %% 2 != 0
  subnetwork <- FindConnectionsForAllHopCountsU(subnetworks = significantSubnetworks,
                                                verbose = verbose, odd = odd)
  return(subnetwork)
}

#analogous to BLOBFISH
#' Find the subnetwork of significant edges n / 2 hops away from each gene.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param network A network with the first column ("source") consisting of source nodes,
#' the second ("target") consisting of target nodes, and the third ("score") consisting of scores
#' @param hopConstraint The maximum number of hops to be considered for a gene.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
FindEdgesForHopU <- function(geneSet, network, hopConstraint,
                             verbose = FALSE, topX = NULL){
  # Build the significant subnetwork for each gene, up to the hop constraint.
  uniqueGeneSet <- sort(unique(geneSet))
  geneSubnetworks <- lapply(uniqueGeneSet, function(gene){
    
    # Get all significant edges for a 1-hop subnetwork.
    if(verbose == TRUE){
      message(paste("Evaluating hop 1 for gene", gene))
    }
    subnetwork1Hop <- BreadthFirstSearchU(networks = network,
                                          startingNodes = gene,
                                          nodesToExclude = c(),
                                          verbose = verbose,
                                          topX = topX)
    
    # Set the starting and excluded set for the next hop.
    excludedSubset <- gene
    startingNodes <- setdiff(unique(c(subnetwork1Hop[,1], subnetwork1Hop[,2])),
                             excludedSubset)
    topXNew <- NULL
    if(!is.null(topX)){
      topXNew <- topX * length(startingNodes)
    }
    
    # Add to the list of all subnetworks.
    allSubnetworksForGene <- list(subnetwork1Hop)
    
    
    # Loop until we reach the maximum number of hops or there are no new edges
    # to traverse.
    hop <- 2
    while(hop <= hopConstraint && length(startingNodes) > 0){
      
      # Do the next hop.
      subnetworkHops <- BreadthFirstSearchU(networks = network,
                                            startingNodes = startingNodes,
                                            nodesToExclude = excludedSubset,
                                            verbose = verbose,
                                            topX = topXNew)
      
      # Set the starting and excluded set for the next hop.
      excludedSubset <- c(excludedSubset, startingNodes)
      startingNodes <- setdiff(unique(c(subnetworkHops[,1], subnetworkHops[,2])), excludedSubset)
      
      
      if(!is.null(topX)){
        topXNew <- topX * length(startingNodes)
      }
      
      # Add to the list.
      allSubnetworksForGene[[length(allSubnetworksForGene) + 1]] <- subnetworkHops
      
      # Increment hops.
      hop <- hop + 1
    }
    return(allSubnetworksForGene)
  })
  
  # Add the names of the genes.
  names(geneSubnetworks) <- uniqueGeneSet
  return(geneSubnetworks)
}

#' Find all edges adjacent to the starting nodes, excluding the nodes
#' specified.
#' @param networks A network with the first column ("source") consisting of source nodes,
#' the second ("target") consisting of target nodes, and the third ("score") consisting of scores#' 
#' @param startingNodes The list of nodes from which to start.
#' @param nodesToExclude The list of nodes to exclude from the search.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
BreadthFirstSearchU <- function(networks, startingNodes,
                                nodesToExclude,
                                verbose = FALSE, topX = NULL){
  
  # Check that provided nodes overlap with the networks.
  if(length(setdiff(startingNodes, c(networks[,1], networks[,2]))) > 0){
    stop("ERROR: Starting nodes do not overlap with network nodes")
  }
  if(length(setdiff(nodesToExclude, c(networks[,1], networks[,2]))) > 0){
    stop("ERROR: List of nodes to exclude does not overlap with network nodes")
  }
  if(length(intersect(startingNodes, nodesToExclude)) > 0){
    stop("ERROR: Starting nodes cannot overlap with nodes to exclude")
  }
  
  
  # Identify genes and transcription factors to test, based on which of these we are
  # starting from.
  genesToTest <- setdiff(unique(c(networks[,1],networks[,2])), nodesToExclude)
  
  # Construct all edges to test based on the combination of these.
  srcGeneLongList <- rep(genesToTest, length(genesToTest))
  tgtGeneLongList <- unlist(lapply(genesToTest, function(gene){
    return(rep(gene, length(genesToTest)))
  }))
  subnetwork <- networks[
    (networks$source %in% startingNodes | networks$target %in% startingNodes) &
      !(networks$source %in% nodesToExclude | networks$target %in% nodesToExclude),
    , ]
  
  # For each edge, measure its significance.
  allEdges <- rownames(subnetwork)
  if(length(allEdges) > 0){
    
    # If topX is specified, filter again.
    significantEdges <- allEdges
    if(!is.null(topX) && length(allEdges) > topX){
      whichTopX <- order(allEdges)[1:topX]
      significantEdges <- allEdges[whichTopX]
    }
    
    # Return the edges.
    subnetwork <- networks[significantEdges, c(1:2)]
    if(verbose == TRUE){
      message(paste("Retained", length(significantEdges), "edges"))
    }
  }
  
  # Return the subnetwork.
  return(subnetwork)
}


#' A helper function for obtaining the connecting subnetwork given two spheres
#' of influence. Called from FindConnectionsForAllHopCountsU().
#' @param subnetwork1 The sphere of influence for the first gene in a pair.
#' @param subnetwork2 The sphere of influence for the second gene in a pair.
#' @param gene1 The first gene in a pair.
#' @param gene2 The second gene in a pair.
#' @param subnetworks The full subnetworks list.
#' @param hops The current hop.
#' @param odd Whether or not the maximum number of hops is an odd number.
#' @param verbose Whether or not to print detailed information about the run.
GetConnectingSubnetworkHelper <- function(subnetwork1, subnetwork2, gene1, gene2, hops,
                                          subnetworks, odd = FALSE,
                                          verbose = FALSE){
  
  # Remove edges in the gene1 neighborhood that cross gene2 and
  # edges in the gene2 neighborhood that cross gene1.
  subnetwork1 <- subnetwork1[which(subnetwork1[,1] != gene2),]
  subnetwork1 <- subnetwork1[which(subnetwork1[,2] != gene2),]
  subnetwork2 <- subnetwork2[which(subnetwork2[,1] != gene1),]
  subnetwork2 <- subnetwork2[which(subnetwork2[,2] != gene1),]
  
  # Initialize overlapping subnetwork.
  genesToRecurse1 <- c()
  genesToRecurse2 <- c()
  
  # Add edges from genes that overlap.
  overlappingGenes <- intersect(c(subnetwork1[,1], subnetwork1[, 2]), 
                                c(subnetwork2[,1], subnetwork2[, 2]))
  
  if(verbose == TRUE){
    message(paste("Hop", hops, "-", length(overlappingGenes), "overlapped between", gene1, "and", gene2))
  }
  whichSubnet1GeneCol1 <- which(subnetwork1[,1] %in% overlappingGenes)
  whichSubnet1GeneCol2 <- which(subnetwork1[,2] %in% overlappingGenes)
  whichSubnet2GeneCol1 <- which(subnetwork2[,1] %in% overlappingGenes)
  whichSubnet2GeneCol2 <- which(subnetwork2[,2] %in% overlappingGenes)
  genesToRecurse1 <- c(subnetwork1[whichSubnet1GeneCol1, 2],
                       subnetwork1[whichSubnet1GeneCol2, 1])
  genesToRecurse2 <- c(subnetwork2[whichSubnet2GeneCol1, 2],
                       subnetwork2[whichSubnet2GeneCol2, 1])
  overlappingEdges1 <- unique(c(rownames(subnetwork1)[whichSubnet1GeneCol1],
                                rownames(subnetwork1)[whichSubnet1GeneCol2]))
  overlappingEdges2 <- setdiff(unique(c(rownames(subnetwork2)[whichSubnet2GeneCol1],
                                        rownames(subnetwork2)[whichSubnet2GeneCol2])),
                               overlappingEdges1)
  
  # Add overlapping edges to subnetwork
  connectingSubnetwork <- rbind(subnetwork1[overlappingEdges1,], subnetwork2[overlappingEdges2,])
  
  # Recurse back over the number of hops.
  if(hops-1 >= 1){
    for(hop in (hops-1):1){
      subnetwork1 <- subnetworks[[gene1]][[hop]]
      subnetwork2 <- subnetworks[[gene2]][[hop]]
      
      # If starting with an odd value, go to the previous hop to recurse.
      if(odd == TRUE && hop - 1 > 0){
        subnetwork2 <- subnetworks[[gene2]][[hop-1]]
      }
      
      # Remove edges in the gene1 neighborhood that cross gene2..
      subnetwork1 <- subnetwork1[which(subnetwork1[,1] != gene2),]
      subnetwork1 <- subnetwork1[which(subnetwork1[,2] != gene2),]
      
      # If the current number of hops is even, add edges from genes connected to genes of interest.
      whichGeneConnectedToGene1Col2 <- which(subnetwork1[,2] %in% genesToRecurse1)
      whichGeneConnectedToGene1Col1 <- which(subnetwork1[,1] %in% genesToRecurse1)
      genesToRecurse1 <- unique(c(subnetwork1[whichGeneConnectedToGene1Col2, 1],
                                  subnetwork1[whichGeneConnectedToGene1Col1, 2]))
      overlappingEdges1 <- unique(c(rownames(subnetwork1)[whichGeneConnectedToGene1Col1],
                                    rownames(subnetwork1)[whichGeneConnectedToGene1Col2]))
      overlappingEdges2 <- c()
      
      # If hop - 1 > 0 or we are starting with an even number, add edges in both directions.
      # Otherwise, add edges from only gene 1.
      if(hop - 1 > 0 || odd == FALSE){
        # Remove edges in the gene1 neighborhood that cross gene1.
        subnetwork2 <- subnetwork2[which(subnetwork2[,1] != gene1),]
        subnetwork2 <- subnetwork2[which(subnetwork2[,2] != gene1),]
        
        # If the current number of hops is even, add edges from genes connected to genes of interest.
        whichGeneConnectedToGene2Col2 <- which(subnetwork2[,2] %in% genesToRecurse2)
        whichGeneConnectedToGene2Col1 <- which(subnetwork2[,1] %in% genesToRecurse2)
        genesToRecurse2 <- unique(c(subnetwork2[whichGeneConnectedToGene2Col2, 1],
                                    subnetwork2[whichGeneConnectedToGene2Col1, 2]))
        overlappingEdges2 <- unique(c(rownames(subnetwork2)[whichGeneConnectedToGene2Col1],
                                      rownames(subnetwork2)[whichGeneConnectedToGene2Col2]))
      }
      
      # Add edges.
      connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[overlappingEdges1,],
                                    subnetwork2[overlappingEdges2,])
    }
  }
  
  # Return.
  return(connectingSubnetwork)
}

#' For all hop counts up to the maximum, find subnetworks connecting each pair of
#' genes by exactly that number of hops. For instance, find each
#' containing only the significant edges meeting the hop count criteria and
#' where each network is a data frame with the following format:
#' @param subnetworks The subnetworks, generated using FindEdgesForHopU().
#' @param odd Whether or not we are looking at an odd number of hops (in which case
#' we should subtract 1 from the hop count)
#' @param verbose Whether or not to print detailed information about the run.
FindConnectionsForAllHopCountsU <- function(subnetworks, odd = FALSE, verbose = FALSE){
  
  # Find a subnetwork for all hops of distance 1.
  compositeSubnetwork <- do.call(rbind, lapply(1:(length(names(subnetworks))-1), function(i){
    genePairSpecificHopCountSubnetwork <- lapply((i+1):length(names(subnetworks)), function(j){
      
      # Get the subnetworks for the number of hops of interest.
      gene1 <- names(subnetworks)[i]
      gene2 <- names(subnetworks)[j]
      subnetwork1 <- subnetworks[[gene1]][[1]]
      oneHopSubnetworkCol1 <- subnetwork1[which(subnetwork1[,1] == gene2),]
      oneHopSubnetworkCol2 <- subnetwork1[which(subnetwork1[,2] == gene2),]
      return(rbind(oneHopSubnetworkCol1, oneHopSubnetworkCol2))
    })
    return(do.call(rbind, genePairSpecificHopCountSubnetwork))
  }))
  
  # Find a subnetwork for each hop count. If we are looking at odd numbers only and
  # we've only done one hop for each gene, that means we are only looking at hop
  # sizes of 1. Stop here.
  if(odd == FALSE || length(subnetworks[[1]]) > 1){
    hopCountSubnetworks <- lapply(1:length(subnetworks[[1]]), function(hops){
      
      # For each pair of genes, find the subnetworks for this number of hops.
      geneSpecificHopCountSubnetwork <- lapply(1:(length(names(subnetworks))-1), function(i){
        genePairSpecificHopCountSubnetwork <- lapply((i+1):length(names(subnetworks)), function(j){
          
          # Get the subnetworks for the number of hops of interest.
          gene1 <- names(subnetworks)[i]
          gene2 <- names(subnetworks)[j]
          connectingSubnetwork <- data.frame(source = NA, target = NA)[0,]
          
          # If there were no edges at this hop count for one or both genes,
          # do not evaluate.
          if(length(subnetworks[[gene1]]) >= hops && length(subnetworks[[gene2]]) >= hops){
            
            # Filter edges (even number of hops).
            if(odd == FALSE){
              subnetwork1 <- subnetworks[[gene1]][[hops]]
              subnetwork2 <- subnetworks[[gene2]][[hops]]
              connectingSubnetwork <- rbind(connectingSubnetwork,
                                            GetConnectingSubnetworkHelper(subnetwork1 = subnetwork1,
                                                                          subnetwork2 = subnetwork2,
                                                                          gene1 = gene1, gene2 = gene2,
                                                                          hops = hops, subnetworks = subnetworks,
                                                                          verbose = verbose))
              if(hops > 1){
                subnetwork2Odd <- subnetworks[[gene2]][[hops-1]]
                connectingSubnetwork <- rbind(connectingSubnetwork,
                                              GetConnectingSubnetworkHelper(subnetwork1 = subnetwork1,
                                                                            subnetwork2 = subnetwork2Odd,
                                                                            gene1 = gene1, gene2 = gene2,
                                                                            hops = hops, subnetworks = subnetworks,
                                                                            verbose = verbose))
              }
            }else if(hops > 1){
              subnetwork1 <- subnetworks[[gene1]][[hops]]
              subnetwork2 <- subnetworks[[gene2]][[hops-1]]
              connectingSubnetwork <- rbind(connectingSubnetwork,
                                            GetConnectingSubnetworkHelper(subnetwork1 = subnetwork1,
                                                                          subnetwork2 = subnetwork2,
                                                                          gene1 = gene1, gene2 = gene2,
                                                                          hops = hops, subnetworks = subnetworks,
                                                                          odd = TRUE, verbose = verbose))
              subnetwork1Even <- subnetworks[[gene1]][[hops-1]]
              connectingSubnetwork <- rbind(connectingSubnetwork,
                                            GetConnectingSubnetworkHelper(subnetwork1 = subnetwork1Even,
                                                                          subnetwork2 = subnetwork2,
                                                                          gene1 = gene1, gene2 = gene2,
                                                                          hops = hops, subnetworks = subnetworks,
                                                                          odd = TRUE, verbose = verbose))
            }
          }
          # Return the subnetwork, which should now contain all of the edges connecting the
          # gene pair at the prespecified number of hops.
          return(connectingSubnetwork)
        })
        
        # Bind together the subnetwork for each gene pair.
        connectingSubnetworkAll <- do.call(rbind, genePairSpecificHopCountSubnetwork)
        return(connectingSubnetworkAll)
      })
      
      # Bind together the subnetworks for each gene.
      return(do.call(rbind, geneSpecificHopCountSubnetwork))
    })
    
    #Fix ordering nd filtering
    compositeSubnetwork <- rbind(compositeSubnetwork, do.call(rbind, hopCountSubnetworks))
    colnames(compositeSubnetwork) <- c("source", "target")
  }
  
  # remove duplicates in real time
  edgeKeys <- paste(compositeSubnetwork$source,
                    compositeSubnetwork$target, sep = "__")
  edgeKeysReverse <- paste(compositeSubnetwork$target,
                           compositeSubnetwork$source, sep = "__")
  compositeSubnetworkDedup <- compositeSubnetwork
  
  # Remove all duplicated edges.
  i = 1
  len <- length(edgeKeys)
  while(i < len){
    whichToRemove <- setdiff(which(edgeKeys == edgeKeys[i]), i)
    if(length(whichToRemove) > 0){
      compositeSubnetworkDedup <- compositeSubnetworkDedup[-whichToRemove,]
      edgeKeys <- edgeKeys[-whichToRemove]
      edgeKeysReverse <- edgeKeysReverse[-whichToRemove]
      len <- length(edgeKeys)
    }
    i <- i + 1
  }
  
  # Remove disconnected chains.
  targetGeneCounts <- table(compositeSubnetworkDedup$target)
  disconnectedGenes <- names(targetGeneCounts)[which(targetGeneCounts == 1)]
  disconnectedGenes <- setdiff(disconnectedGenes, names(subnetworks))
  while(length(disconnectedGenes) > 0){
    disconnectedChainEdges <- c(which(compositeSubnetworkDedup$target %in% disconnectedGenes),
                                which(compositeSubnetworkDedup$source %in% disconnectedGenes))
    otherEdges <- sort(setdiff(1:nrow(compositeSubnetworkDedup), disconnectedChainEdges))
    compositeSubnetworkDedup <- compositeSubnetworkDedup[otherEdges,]
    targetGeneCounts <- table(compositeSubnetworkDedup$target)
    disconnectedGenes <- names(targetGeneCounts)[which(targetGeneCounts == 1)]
    disconnectedGenes <- setdiff(disconnectedGenes, names(subnetworks))
  }
  
  # Remove all reversed edges.
  edgeKeys <- paste(compositeSubnetworkDedup$source,
                    compositeSubnetworkDedup$target, sep = "__")
  edgeKeysReverse <- paste(compositeSubnetworkDedup$target,
                           compositeSubnetworkDedup$source, sep = "__")
  i = 1
  len <- length(edgeKeys)
  while(i < len){
    if(edgeKeysReverse[i] %in% edgeKeys){
      whichToRemove <- intersect(which(edgeKeys == edgeKeysReverse[i]),
                                 (i+1):length(edgeKeys))
      compositeSubnetworkDedup <- compositeSubnetworkDedup[-whichToRemove,]
      edgeKeys <- edgeKeys[-whichToRemove]
      edgeKeysReverse <- edgeKeysReverse[-whichToRemove]
      len <- length(edgeKeys)
    }
    i <- i + 1
  }
  
  # Set row names.
  rownames(compositeSubnetworkDedup) <- paste(compositeSubnetworkDedup[,1], compositeSubnetworkDedup[,2], sep = "__")
  return(compositeSubnetworkDedup)
}

#' Plot the networks, using different colors for transcription factors, genes of interest,
#' and additional genes.
#' @param network A data frame with the following format:
#' tf,gene
#' @param geneColorMapping Color mapping from a set of genes to a color. The
#' nodes and edges connected to them will be this color. If NULL, all genes and
#' their edges will be gray. The format is a data frame, where the first column ("gene")
#' is the name of the gene and the second ("color") is the color. All nodes must
#' be included in the color mapping.
#' @param nodeSize Size of node
#' @param edgeWidth Width of edges
#' @param vertexLabels Vector of labels to include. Labels that you do not wish to include
#' should be set to NA.
#' @param vertexLabelSize The size of label to use for the vertex, as a fraction of the default.
#' @param vertexLabelOffset Number of pixels in the offset when plotting labels.
#' Default is TRUE.
#' @export
PlotNetworkU <- function(network, nodeSize = 1,
                         edgeWidth = 0.5, vertexLabels = NA, vertexLabelSize = 0.7,
                         vertexLabelOffset = 0.5, geneColorMapping = NULL){
  
  # Convert from factor to character.
  colnames(network)[1:2] <- c("source", "target")
  network$source <- as.character(network$source)
  network$target <- as.character(network$target)
  
  # Set the node attributes.
  uniqueNodes <- unique(c(network$source, network$target))
  nodeAttrs <- data.frame(node = uniqueNodes,
                          color = rep("gray", length(uniqueNodes)),
                          size = rep(nodeSize, length(uniqueNodes)),
                          frame.width = rep(0, length(uniqueNodes)),
                          label.color = "black", label.cex = vertexLabelSize,
                          label.dist = vertexLabelOffset)
  rownames(nodeAttrs) <- uniqueNodes
  
  # Add colors.
  nodeAttrs$color <- "gray"
  if (!is.null(geneColorMapping)) {
    rownames(geneColorMapping) <- geneColorMapping[,1]
    nodeAttrs[rownames(nodeAttrs) %in% rownames(geneColorMapping), "color"] <-
      geneColorMapping[rownames(nodeAttrs)[rownames(nodeAttrs) %in% rownames(geneColorMapping)], "color"]
  }
  
  # Add gene colors.
  rownames(geneColorMapping) <- geneColorMapping[,1]
  geneColorMapping <- geneColorMapping[intersect(rownames(geneColorMapping), uniqueNodes),]
  if(!is.null(geneColorMapping)){
    nodeAttrs[rownames(geneColorMapping), "color"] <- geneColorMapping$color
  }
  # Add edge attributes.
  if(!is.null(geneColorMapping)){
    rownames(geneColorMapping) <- geneColorMapping$gene
    uniqueColors <- unique(geneColorMapping$color)
    edgeColors <- unlist(lapply(1:nrow(network), function(i){
      retcol <- "gray"
      if(geneColorMapping[network$source[i], "color"] == geneColorMapping[network$target[i], "color"]){
        retcol <- geneColorMapping[network$source[i], "color"]
      }
      return(retcol)
    }))
    network$color <- edgeColors
  }
  network$width <- edgeWidth
  # Create a graph object.
  graph <- igraph::graph_from_data_frame(network, vertices = nodeAttrs, directed = FALSE)
  V(graph)$type <- V(graph)$name %in% network$source
  
  # Plot.
  if(length(vertexLabels) == 1 && is.na(vertexLabels)){
    vertexLabels <- V(graph)$name
  }
  plot(graph, vertex.label = vertexLabels)
}
