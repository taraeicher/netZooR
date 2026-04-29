monsterAnalysis <- setClass("monsterAnalysis", slots=c("tm","nullTM","numGenes","numSamples", "logging"))
setMethod("show","monsterAnalysis",function(object){monsterPrintMonsterAnalysis(object)})

#' monsterGetTm
#'
#' acessor for the transition matrix in MONSTER object
#'
#' @param x an object of class "monsterAnalysis"
#' @export
#' @return Transition matrix
#' @examples
#' data(monsterRes)
#' tm <- monsterGetTm(monsterRes)
monsterGetTm <- function(x){
    x@tm
}

#' monsterPlotMonsterAnalysis
#'
#' plots the sum of squares of off diagonal mass (differential TF Involvement)
#'
#' @param x an object of class "monsterAnalysis"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return Plot of the dTFI for each TF against null distribution
#' @examples
#' data(yeast)
#' yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' design <- c(rep(1,25),rep(0,10),rep(NA,15))
#' #monsterRes <- monster(yeast$exp.cc, design,
#' #yeast$motif, nullPerms=10, numMaxCores=1)
#' #monsterPlotMonsterAnalysis(monsterRes)
monsterPlotMonsterAnalysis <- function(x, ...){
  monsterdTFIPlot(x,...)
}

#' monsterPrintMonsterAnalysis
#'
#' summarizes the results of a MONSTER analysis
#'
#' @param x an object of class "monster"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return Description of transition matrices in object
#' @examples
#' data(yeast)
#' yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' design <- c(rep(1,25),rep(0,10),rep(NA,15))
#' #monster(yeast$exp.cc,design,yeast$motif, nullPerms=10, numMaxCores=1)
monsterPrintMonsterAnalysis <- function(x, ...){
  if(x@logging == TRUE){
    cat("MONSTER object\n")
    cat(paste(x@numGenes, "genes\n"))
    cat(paste(x@numSamples[1],"baseline samples\n"))
    cat(paste(x@numSamples[2],"final samples\n"))
    cat(paste("Transition driven by", ncol(x@tm), "transcription factors\n"))
    cat(paste("Run with", length(x@nullTM), "randomized permutations.\n"))
  }
}

#' MOdeling Network State Transitions from Expression and Regulatory data (MONSTER)
#'
#' This function runs the MONSTER algorithm.  Biological states are characterized by distinct patterns 
#' of gene expression that reflect each phenotype's active cellular processes. 
#' Driving these phenotypes are gene regulatory networks in which transcriptions factors control 
#' when and to what degree individual genes are expressed. Phenotypic transitions, such as those that 
#' occur when disease arises from healthy tissue, are associated with changes in these  networks. 
#' MONSTER is an approach to understanding these transitions. MONSTER models phenotypic-specific 
#' regulatory networks and then estimates a "transition matrix" that converts one state to another. 
#' By examining the properties of the transition matrix, we can gain insight into regulatory 
#' changes associated with phenotypic state transition.
#' Important note: the direct regulatory network observed from gene expression is currently
#' implemented as a regular correlation as opposed to the partial correlation described 
#' in the paper.
#' There are 2 modes to run MONSTER:
#' (1) MONSTER can internally estimate a gene regulatory network, in which case 
#' it will take as input an expression data plus a motif network
#' (2) MONSTER can run on pre-computed gene regulatory networks (e.g., from PANDA),
#' in which case the `motif` argument is set to `NA`, and the `expr` argument will
#' contain the concatenated gene regulatory networks. This mode can be selected
#' by setting the `mode` argument to `regNet`.
#' Alternatively, see the `domonster` function for a quick-start way to run this mode.
#' Citation: Schlauch, Daniel, et al. "Estimating drivers of cell state transitions using gene regulatory network models." 
#' BMC systems biology 11.1 (2017): 139. https://doi.org/10.1186/s12918-017-0517-y
#' @param expr Gene Expression dataset, can be matrix or data.frame of expression values or ExpressionSet. 
#' If `mode` is set to `regNet`, MONSTER will be run in pre-computed gene regulatory network mode, in which case
#' gene regulatory networks can be passed for this argument. See also `domonster` to use this mode.
#' @param design Binary vector indicating case control partition. 1 for case and 0 for control.
#' @param motif Regulatory data.frame consisting of three columns.  For each row, a transcription factor (column 1) 
#' regulates a gene (column 2) with a defined strength (column 3), usually taken to be 0 or 1.
#' May also be NA, if MONSTER is being run on pre-computed gene regulatory networks, passed in the `expr` argument.
#' See also `domonster` to use this mode.
#' @param nullPerms number of random permutations to run (default 100).  Set to 0 to only 
#' calculate observed transition matrix. When mode is is 'buildNet' it randomly permutes the case and control expression
#' samples, if mode is 'regNet' it will randomly permute the case and control networks.
#' @param ni_method String to indicate algorithm method.  Must be one of "BERE","pearson",or "lda". Default is "BERE"
#' @param ni.coefficient.cutoff numeric to specify a p-value cutoff at the network
#' inference step.  Default is NA, indicating inclusion of all coefficients.
#' @param numMaxCores requires doParallel, foreach.  Runs MONSTER in parallel computing 
#' environment.  Set to 1 to avoid parallelization, NA will take the default parallel pool in the computer.
#' @param outputDir character vector specifying a directory or path in which 
#' which to save MONSTER results, default is NA and results are not saved.
#' @param alphaw A weight parameter between 0 and 1 specifying proportion of weight 
#' to give to indirect compared to direct evidence. The default is 0.5 to give an 
#' equal weight to direct and indirect evidence.
#' @param mode A parameter telling whether to build the regulatory networks ('buildNet') or to use provided regulatory networks
#' ('regNet'). If set to 'regNet', then the parameters motif, ni_method, ni.coefficient.cutoff, and alphaw will be set to NA. Gene regulatory
#' networks are supplied in the 'expr' variable as a TF-by-Gene matrix, by concatenating the TF-by-Gene matrices of case and control, expr has size nTFs x 2nGenes.
#' @param method Method to use in computing the transition matrix. These include "ols" (default),"kabsch","L1", and "orig" (SVD).
#' @param remove.diagonal #' Logical for returning a result containing 0s across the diagonal (default = TRUE).
#' @param nullModelType Type of null model used. If "permutation" is set as the default, MONSTER permutes expression data when in "buildNet" mode
#' and permutes network edges by gene when in "regNet" mode. "nullNetwork" generates null networks from the control network using the edge weights from
#' the control network as the mean and the edge weight variances across networks with no true signal as the variance.
#' @param nullNetworks This parameter is only used when running in "regNet" mode with the "nullNetwork" option set. This includes an H5 format file
#' containing an object named "nullNetworks"
#' or data frame with the values from networks with no signal. The first two columns must be the source and target nodes, and the remaining columns must be
#' the edge weights from each simulated network. For PANDA networks:
#' tf, gene, net1, net2, ..., netn
#' @param logging Whether or not to print logging messages for MONSTER (default is TRUE)
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom methods new
#' @return An object of class "monsterAnalysis" containing results
#' 
#' @examples
#' # Example with the network reconstruction step
#' data(yeast)
#' design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # Example with provided networks
#' \donttest{
#' pandaResult <- panda(pandaToyData$motif, pandaToyData$expression, pandaToyData$ppi)
#' case=pandaResult@regNet
#' nelemReg=dim(pandaResult@regNet)[1]*dim(pandaResult@regNet)[2]
#' nGenes=length(colnames(pandaResult@regNet))
#' control=matrix(rexp(nelemReg, rate=.1), ncol=nGenes)
#' colnames(control) = colnames(case)
#' rownames(control) = rownames(case) 
#' expr = as.data.frame(cbind(control,case))
#' design=c(rep(0,nGenes),rep(1, nGenes))
#' monsterRes <- monster(expr, design, motif=NA, nullPerms=10, numMaxCores=1, mode='regNet')
#' }
#' # alternatively, if a gene regulatory network has already been estimated,
#' # see the domonster function for quick start
monster <- function(expr, 
                    design, 
                    motif, 
                    nullPerms=100,
                    ni_method="BERE",
                    ni.coefficient.cutoff = NA,
                    numMaxCores=1, 
                    outputDir=NA, alphaw=0.5, mode='buildNet',
                    method="ols", remove.diagonal = TRUE,
                    nullModelType = "permutation",
                    nullNetworks = NA, logging = TRUE){
  if(!(ni_method %in% c("BERE","pearson", "lda"))){
    stop(paste("Supported values for ni_method are BERE, pearson, and lda.",
               ni_method, "is invalid."))
  }
  # Check for correct inputs.
  if(!nullModelType %in% c("permutation", "nullNetwork")){
    stop("Only 'permutation' and 'nullNetwork' are permitted values for nullModelType.")
  }
  if(nullModelType == "nullNetwork"){
    if(mode == "buildNet"){
      stop("Cannot run 'nullNetwork' model type in buildNet mode")
    }
    if(length(nullNetworks) == 1 && is.na(nullNetworks)){
      stop(paste("Must provide null networks when null model type is 'nullNetwork'. These can",
                 "be generated using the function GenerateNullPANDADistribution() if you are",
                 "using PANDA networks."))
    }
    if(is.data.frame(nullNetworks) && !(all(sapply(nullNetworks[1:2], is.character), sapply(nullNetworks[-(1:2)], is.numeric)))){
      stop(paste("The null network provided must have source and target nodes in the first two columns",
                 "and numeric scores in all remaining columns."))
    }else if(is.character(nullNetworks) && !file.exists(nullNetworks)){
      stop(paste("File", nullNetworks, "does not exist!"))
    }
    if(is.data.frame(nullNetworks)){
      exprNodes <- sort(unique(c(rownames(expr[,design == 0]), colnames(expr[,design == 0]))))
      nullNodes <- sort(unique(c(nullNetworks[,1], nullNetworks[,2])))
      if(length(setdiff(exprNodes, nullNodes)) > 0 || length(setdiff(nullNodes, exprNodes)) > 0){
        stop("The node lists differ between the input networks and the provided null networks.")
      }
    }
  }
  
  if(mode=='regNet'){
    motif=NA
    alphaw=NA
    ni_method=NA
    ni.coefficient.cutoff=NA
    if(length(design == 1) != length(design == 0)){
      stop('case and control have a different number of genes')
    }
  }else{
    if(is.null(motif)){
      stop("motif may not be NULL")
    }
    if(length(motif) == 1 && is.na(motif)){
      stop('Set mode to "regNet" if using as input pre-made regulatory networks.\n
           Otherwise, motif should not be NA if using buildNet mode.')
    }
  }
  # Data type checking
  expr <- monsterCheckDataType(expr)
  # Parallelize
  # Initiate cluster
  if(!is.na(numMaxCores) && numMaxCores > 1){
    # Calculate the number of cores
    numCores <- detectCores() - 4
    numCores <- min(numCores, numMaxCores)
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    if(logging == TRUE){
      print("Running null permutations in parallel")
      print(paste(numCores,"cores used"))
    }
  }
  
  iters <- nullPerms+1 # Two networks for each partition, plus observed partition
  if(logging == TRUE){
    print(paste(iters,"network transitions to be estimated"))
  }
  
  #start time
  strt  <- Sys.time()
  #loop
  if(!is.na(outputDir)){
    dir.create(file.path(outputDir))  
    dir.create(file.path(outputDir,"tms"))        
  }
  
  # Remove unassigned data 
  expr <- expr[,design%in%c(0,1)]
  design <- design[design%in%c(0,1)]
  
  # Check column order
  if(mode == 'regNet'){
    numGenes = ncol(expr)/2
    if(any(colnames(expr[,design==1]) != colnames(expr[,design==0]))){
      stop('Please provide two regulatory networks with the same gene labels and 
           the same number of genes in the same order')
    }
  }else if(mode == 'buildNet'){
    numGenes = nrow(expr)
  }
  
  # If nullNetworks is a file, we need to run the analysis one column at a time
  # to avoid memory errors. Otherwise, we can build the null networks first and
  # then run the analysis.
  if(is.character(nullNetworks)){
    transMatrices <- CalculateTransitionMatricesFromFile(file = nullNetworks, mode = mode,
                                                         expr = expr, iterations = nullPerms + 1,
                                                         design = design, logging = logging,
                                                         numMaxCores = numMaxCores,
                                                         motif = motif, 
                                                         ni_method = ni_method,
                                                         ni.coefficient.cutoff = ni.coefficient.cutoff, 
                                                         alphaw = alphaw,
                                                         remove.diagonal = remove.diagonal, 
                                                         method = method, 
                                                         outputDir = outputDir)
  }else{
    # Generate null.
    nullExprAll <- NULL
    if(nullModelType == "permutation" || mode == "buildNet"){
      nullExprAll <- GeneratePermutationNull(expr = expr, iterations = nullPerms, mode = mode,
                                             logging = logging)
    }else if(nullModelType == "nullNetwork"){
      nullExprAll <- GenerateNullFromControl(concatNet = expr, iterations = nullPerms,
                                             design = design, nullNetworks = nullNetworks,
                                             logging = logging)
    }
    
    # Calculate transition matrix for the true case/control data and for the
    # null data. Do this in a for loop if we are using one core or a foreach loop if
    # we are using multiple cores.
    transMatrices=list()
    if(numMaxCores == 1){
      for(i in seq_len(iters)){
        nullExpr <- NULL
        if(i!=1){
          nullExpr <- nullExprAll[[i-1]]
        }else{
          nullExpr <- expr
        }
        transMatrices[[i]] <- CalculateOneTransitionMatrix(i = i, mode = mode, 
                                                           nullExpr = nullExpr, 
                                                           motif = motif, design = design,
                                                           ni_method = ni_method,
                                                           ni.coefficient.cutoff = ni.coefficient.cutoff, 
                                                           alphaw = alphaw,
                                                           remove.diagonal = remove.diagonal, 
                                                           method = method, 
                                                           outputDir = outputDir,
                                                           logging = logging)
      }
    }else{
      transMatrices <- foreach(i=seq_len(iters),
                               .packages=c("netZooR","reshape2","penalized","MASS")) %dopar% {
                                 nullExpr <- NULL
                                 if(i!=1){
                                   nullExpr <- nullExprAll[[i-1]]
                                 }else{
                                   nullExpr <- expr
                                 }
                                 return(CalculateOneTransitionMatrix(i = i, mode = mode,
                                                                     nullExpr = nullExpr,
                                                                     motif = motif, design = design,
                                                                     ni_method = ni_method,
                                                                     ni.coefficient.cutoff = ni.coefficient.cutoff,
                                                                     alphaw = alphaw,
                                                                     remove.diagonal = remove.diagonal,
                                                                     method = method,
                                                                     outputDir = outputDir,
                                                                     logging = logging))
                               }
    }
    
    # Log the time.
    if(logging == TRUE){
      print(Sys.time()-strt)
    }
  }
  
  if(!is.na(numMaxCores)  && numMaxCores > 1){
    stopCluster(cl)
  }
  
  gc()
  return(
    monsterAnalysis(
      tm=transMatrices[[1]], 
      nullTM=transMatrices[-1], 
      numGenes=numGenes, 
      numSamples=c(sum(design==0), sum(design==1)),
      logging = logging))
}

#' Calculates a single transition matrix given case and control data.
#' @param i Iteration
#' @param mode "buildNet" or "regNet"
#' @param nullExpr Gene expression dataset
#' @param motif Regulatory data.frame consisting of three columns.  For each row, a transcription factor (column 1) 
#' regulates a gene (column 2) with a defined strength (column 3), usually taken to be 0 or 1.
#' May also be NA, if MONSTER is being run on pre-computed gene regulatory networks, passed in the `expr` argument.
#' See also `domonster` to use this mode.
#' @param ni_method String to indicate algorithm method.  Must be one of "BERE","pearson",or "lda". Default is "BERE"
#' @param ni.coefficient.cutoff numeric to specify a p-value cutoff at the network
#' inference step.  Default is NA, indicating inclusion of all coefficients.
#' @param alphaw A weight parameter between 0 and 1 specifying proportion of weight 
#' to give to indirect compared to direct evidence. The default is 0.5 to give an 
#' equal weight to direct and indirect evidence.
#' @param design The design matrix
#' @param remove.diagonal #' Logical for returning a result containing 0s across the diagonal (default = TRUE).
#' @param method character specifying which algorithm to use, default='ols'.
#' @param outputDir character vector specifying a directory or path in which 
#' which to save MONSTER results, default is NA and results are not saved.
#' @param logging Whether or not to print logging messages for MONSTER (default is TRUE)
#' @return A list of matrices, each one corresponding to one null model
CalculateOneTransitionMatrix <- function(i, mode, nullExpr, motif, ni_method,
                                         ni.coefficient.cutoff, alphaw, design,
                                         remove.diagonal, method, outputDir, logging){
  if(logging == TRUE){
    print(paste0("Running iteration ", i))
  }
  if(mode == 'buildNet'){
    nullExprCases <- nullExpr[,design==1]
    nullExprControls <- nullExpr[,design==0]
    
    tmpNetCases <- monsterMonsterNI(motif, nullExprCases, 
                                    method=ni_method, regularization="none",
                                    score="none", ni.coefficient.cutoff,
                                    verbose=FALSE, randomize = "none", cpp=FALSE,
                                    alphaw, logging = logging)
    tmpNetControls <- monsterMonsterNI(motif, nullExprControls, 
                                       method=ni_method, regularization="none",
                                       score="none", ni.coefficient.cutoff,
                                       verbose=FALSE, randomize = "none", cpp=FALSE,
                                       alphaw, logging = logging)
  }else if(mode == 'regNet'){
    tmpNetCases    = nullExpr[,design==1]
    tmpNetControls = nullExpr[,design==0]
  }
  transitionMatrix <- monsterTransformationMatrix(
    tmpNetControls, tmpNetCases, remove.diagonal=remove.diagonal, method=method,
    logging = logging) 
  if(logging == TRUE){
    print(paste("Finished running iteration", i))
  }
  
  if (!is.na(outputDir)){
    saveRDS(transitionMatrix,file.path(outputDir,'tms',paste0('tm_',i,'.rds')))
  }
  return(transitionMatrix)
}
#' Generates null models for MONSTER by permuting the expression levels
#' for the case and control data.
#' @param expr Gene expression dataset
#' @param iterations Number of null models to generate
#' @param mode "buildNet" or "regNet"
#' @param logging Whether or not to print logging messages for MONSTER (default is TRUE)
#' @return A list of matrices, each one corresponding to one null model
GeneratePermutationNull <- function(expr, iterations, mode, logging){
  retvals <- NULL
  if(mode == "regNet"){
    retvals <- lapply(1:iterations, function(i){
      expr[] <- expr[,sample(seq_along(colnames(expr)))]
      return(expr)
    })
  }else if(mode == "buildNet"){
    exprMat <- as.matrix(expr)
    retvals <- lapply(1:iterations, function(i){
      exprMat[] <- exprMat[sample(seq_along(c(exprMat)))]
      return(exprMat)
    })
  }
  if(logging == TRUE){
    print("Finished generating null models")
  }
  return(retvals)
}

#' Generates null model for MONSTER by re-centering null networks such that the
#' mean values are equivalent to the control network.
#' @param concatNet The concatenated case and control networks
#' @param iterations Number of null models to generate
#' @param design The design matrix
#' @param nullNetworks The null networks from which to calculate the variance
#' @param logging Whether or not to print logging messages for MONSTER (default is TRUE)
#' @return A list of matrices, each one corresponding to one null model
GenerateNullFromControl <- function(concatNet, iterations, design, nullNetworks, logging){

  # Check that the number of iterations is equal to the number of null networks.
  nullNum <- ncol(nullNetworks) - 2
  if(iterations != nullNum){
    warning(paste("Warning: You have specified", iterations, "iterations but",
                  "provided", nullNum, "null networks. MONSTER",
                  "will generate", nullNum, "iterations instead."))
    iterations <- nullNum
  }
  
  # Separate cases from controls.
  concatNet <- as.data.frame(concatNet)
  controls = concatNet[,design==0]

  # Format control matrix for comparison with null networks.
  controlsToMelt <- controls
  controlsToMelt$tf <- rownames(controls)
  controlMelt <- reshape2::melt(
    controlsToMelt,
    id.vars = "tf",
    variable.name = "gene",
    value.name = "score"
  )
  controlMelt$gene <- as.character(controlMelt$gene)
  rownames(controlMelt) <- paste(controlMelt$tf, controlMelt$gene, sep = "__")
  rownames(nullNetworks) <- paste(nullNetworks$tf, nullNetworks$gene, sep = "__")
  controlMelt <- controlMelt[rownames(nullNetworks),]

  # Find the difference in means between the control network and the null networks.
  nullNetworkScores <- nullNetworks[,3:ncol(nullNetworks)]
  nullMeans <- apply(nullNetworkScores, 1, mean)
  meanDiff <- controlMelt[,3] - nullMeans

  # Re-center the null networks.
  meanDiffMat <- as.data.frame(matrix(rep(meanDiff, (ncol(nullNetworks) - 2)), nrow = length(meanDiff)))
  controlNulls <- nullNetworks
  controlNulls[,3:ncol(controlNulls)] <- controlNulls[,3:ncol(controlNulls)] + meanDiffMat
  
  # Format as a list of networks.
  nullFormatted <- lapply(setdiff(colnames(controlNulls), c("tf", "gene")), function(c){
    meltedNull <- data.frame(tf = controlNulls$tf, gene = controlNulls$gene, value = controlNulls[,c])
    castNull <- reshape2::acast(meltedNull, tf ~ gene, value.var = "value")
    castNull <- castNull[rownames(concatNet), colnames(concatNet[,design==0])]
    concatNet[,design==1] <- castNull
    concatNetMat <- as.matrix(concatNet)
    return(concatNetMat)
  })

  # Return.
  if(logging == TRUE){
    print("Finished generating null models")
  }
  return(nullFormatted)
}

#' Calculates transition matrices using the following steps (to save memory)
#' for the real expression and each null expression matrix:
#'    1. Read in the column of the file corresponding to the null expression.
#'.   2. Center the null expression around the control values for a "noisy control".
#'.   3. Find the transition matrix between the noisy control and the cases.
#' @param expr The concatenated gene expression data / case and control networks
#' @param iterations Number of null models to generate
#' @param design The design matrix
#' @param file The file containing the null networks
#' @param mode A parameter telling whether to build the regulatory networks ('buildNet') or to use provided regulatory networks
#' ('regNet'). If set to 'regNet', then the parameters motif, ni_method, ni.coefficient.cutoff, and alphaw will be set to NA. Gene regulatory
#' networks are supplied in the 'expr' variable as a TF-by-Gene matrix, by concatenating the TF-by-Gene matrices of case and control, expr has size nTFs x 2nGenes.
#' @param logging Whether or not to print logging messages for MONSTER (default is TRUE)
#' @param numMaxCores requires doParallel, foreach.  Runs MONSTER in parallel computing 
#' environment.  Set to 1 to avoid parallelization, NA will take the default parallel pool in the computer.
#' @param motif Regulatory data.frame consisting of three columns.  For each row, a transcription factor (column 1) 
#' regulates a gene (column 2) with a defined strength (column 3), usually taken to be 0 or 1.
#' May also be NA, if MONSTER is being run on pre-computed gene regulatory networks, passed in the `expr` argument.
#' See also `domonster` to use this mode.
#' @param ni_method String to indicate algorithm method.  Must be one of "BERE","pearson",or "lda". Default is "BERE"
#' @param ni.coefficient.cutoff numeric to specify a p-value cutoff at the network
#' inference step.  Default is NA, indicating inclusion of all coefficients.
#' @param alphaw A weight parameter between 0 and 1 specifying proportion of weight 
#' to give to indirect compared to direct evidence. The default is 0.5 to give an 
#' equal weight to direct and indirect evidence.
#' @param remove.diagonal #' Logical for returning a result containing 0s across the diagonal (default = TRUE).
#' @param method character specifying which algorithm to use, default='ols'.
#' @param outputDir character vector specifying a directory or path in which 
#' which to save MONSTER results, default is NA and results are not saved.
#' @return A list of matrices, each one corresponding to one null model
CalculateTransitionMatricesFromFile <- function(file, mode,
                                                expr, iterations,
                                                design, logging, numMaxCores,
                                                motif, ni_method, ni.coefficient.cutoff, 
                                                alphaw, remove.diagonal, 
                                                method, outputDir){
  
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' is required for monster() when a file input is provided for nullNetworks.\n",
         "Install using BiocManager::install('rhdf5')", call. = FALSE)
  }
  print("starting")
  # Melt the expression data.
  control <- expr[,design == 0]
  case <- expr[,design == 1]
  
  # Open the file and extract the relevant columns.
  tf <- rhdf5::h5read(file, "matrices/tfs")
  gene <- rhdf5::h5read(file, "matrices/genes")
  sharedTF <- intersect(tf, rownames(control))
  sharedGene <- intersect(gene, colnames(control))
  controlShared <- control[sharedTF, sharedGene]
  caseShared <- case[sharedTF, sharedGene]
  tfLoc <- unlist(lapply(sharedTF, function(tfi){return(which(tf == tfi))}))
  geneLoc <- unlist(lapply(sharedGene, function(genei){return(which(gene == genei))}))

  # Extract running mean for the null network data.
  runningSum <- matrix(rep(0, nrow(controlShared) * ncol(controlShared)), nrow = nrow(controlShared))
  for(i in (seq_len(iterations-1))){
    if(logging == TRUE){
      print(paste("Calculating running mean - iteration", i))
    }
    nullMat <- rhdf5::h5read(file, paste0("matrices/", i))
    nullMatShared <- nullMat[tfLoc, geneLoc]
    runningSum <- runningSum + nullMatShared
  }
  runningMean <- runningSum / (iterations-1)

  transMatrices=list()
  if(numMaxCores == 1){
    for(i in seq_len(iterations)){

      # Create the input matrix for MONSTER containing control and case data.
      fullExprMat <- as.matrix(cbind(controlShared, caseShared))
      newDesign <- c(rep(0, ncol(controlShared)), rep(1, ncol(caseShared)))
      
      # Modify expression matrix using null data if first iteration is done.
      if(i > 1){
        
        # Extract the relevant null data column.
        nullMat <- rhdf5::h5read(file, paste0("matrices/", i-1))
        
        # Subset.
        nullMatShared <- nullMat[tfLoc, geneLoc]

        # Center around control data.
        meanDiff <- controlShared - runningMean
        controlNulls <- nullMatShared + meanDiff

        # Create the input matrix for MONSTER containing control and noisy control data (null).
        fullExprMat <- as.matrix(cbind(controlShared, controlNulls))
        newDesign <- c(rep(0, ncol(controlShared)), rep(1, ncol(controlNulls)))
        
      }
      if(logging == TRUE){
        print(paste("Calculating transition matrix - iteration", i))
      }
      transMatrices[[i]] <- CalculateOneTransitionMatrix(i = i, mode = mode, 
                                                         nullExpr = fullExprMat, 
                                                         motif = motif, 
                                                         ni_method = ni_method,
                                                         ni.coefficient.cutoff = ni.coefficient.cutoff, 
                                                         alphaw = alphaw, design = newDesign,
                                                         remove.diagonal = remove.diagonal, 
                                                         method = method, 
                                                         outputDir = outputDir,
                                                         logging = logging)
    }
  }else{
    transMatrices <- foreach(i=seq_len(iterations),
                             .packages=c("netZooR","reshape2","penalized","MASS")) %dopar% {
                               
                               # Create the input matrix for MONSTER containing control and case data.
                               fullExprMat <- as.matrix(cbind(controlShared, caseShared))
                               newDesign <- c(rep(0, ncol(controlShared)), rep(1, ncol(caseShared)))
                               
                               # Modify expression matrix using null data if first iteration is done.
                               if(i > 1){
                                 
                                 # Extract the relevant null data column.
                                 columns <- col_order[i + 1]
                                 nullDat <- rhdf5::h5read(file, paste0("nullNetworks/", columns))
                                 
                                 # Reshape so that we have a matrix.
                                 nullDF <- data.frame(tf = tf, gene = gene, score = nullDat)
                                 nullMat <- xtabs(score ~ tf + gene, data = nullDF)
                                 nullMatShared <- nullMat[sharedTF, sharedGene]
                                 
                                 # Center around control data.
                                 meanDiff <- controlShared - runningMean
                                 controlNulls <- nullMatShared + meanDiff
                                 
                                 # Create the input matrix for MONSTER containing control and noisy control data (null).
                                 fullExprMat <- as.matrix(cbind(controlShared, controlNulls))
                                 newDesign <- c(rep(0, ncol(controlShared)), rep(1, ncol(controlNulls)))
                                 
                               }
                               return(CalculateOneTransitionMatrix(i = i, mode = mode, 
                                                                   nullExpr = fullExprMat, 
                                                                   motif = motif, 
                                                                   ni_method = ni_method,
                                                                   ni.coefficient.cutoff = ni.coefficient.cutoff, 
                                                                   alphaw = alphaw, design = newDesign,
                                                                   remove.diagonal = remove.diagonal, 
                                                                   method = method, 
                                                                   outputDir = outputDir,
                                                                   logging = logging))
                             }
  }
  # Return the transition matrices.
  return(transMatrices)
}


#' Checks that data is something MONSTER can handle
#'
#' @param expr Gene Expression dataset
#' @return expr Gene Expression dataset in the proper form (may be the same as input)
#' @importFrom methods is
#' @export
#' @examples
#' expr.matrix <- matrix(rnorm(2000),ncol=20)
#' monsterCheckDataType(expr.matrix)
#' #TRUE
#' data(yeast)
#' class(yeast$exp.cc)
#' monsterCheckDataType(yeast$exp.cc)
#' #TRUE
monsterCheckDataType <- function(expr){
  if (!(is.data.frame(expr) || is.matrix(expr) || is(expr, "ExpressionSet")))
    stop("expr must be a data.frame, matrix, or ExpressionSet")
  if("ExpressionSet" %in% class(expr)){
    if (requireNamespace("Biobase", quietly = TRUE)) {
      expr <- Biobase::exprs(expr)
    } 
  }
  if(is.data.frame(expr)){
    expr <- as.matrix(expr)
  }
  return(expr)
}

globalVariables("i")

#' Bi-partite network analysis tools
#'
#' This function analyzes a bi-partite network.
#'
#' @param network.1 starting network, a genes by transcription factors data.frame with scores 
#' for the existence of edges between
#' @param network.2 final network, a genes by transcription factors data.frame with scores 
#' for the existence of edges between
#' @param by.tfs logical indicating a transcription factor based transformation.    If 
#' false, gives gene by gene transformation matrix
#' @param remove.diagonal logical for returning a result containing 0s across the diagonal
#' @param standardize logical indicating whether to standardize the rows and columns
#' @param method character specifying which algorithm to use, default='ols'
#' @param logging Whether or not to log progress (default = TRUE)
#' @return matrix object corresponding to transition matrix
#' @import MASS
#' @export
#' @examples
#' data(yeast)
#' cc.net.1 <- monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,1:20])
#' cc.net.2 <- monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,31:50])
#' monsterTransformationMatrix(cc.net.1, cc.net.2)

monsterTransformationMatrix <- function(network.1, network.2, by.tfs=TRUE, standardize=FALSE, 
                                          remove.diagonal=TRUE, method="ols", logging = TRUE){
  if(is.list(network.1)&&is.list(network.2)){
    if(by.tfs){
      net1 <- t(network.1$reg.net)
      net2 <- t(network.2$reg.net)
    } else {
      net1 <- network.1$reg.net
      net2 <- network.2$reg.net
    }
  } else if(is.matrix(network.1)&&is.matrix(network.2)){
    if(by.tfs){
      net1 <- t(network.1)
      net2 <- t(network.2)
    } else {
      net1 <- network.1
      net2 <- network.2
    }
  } else {
    stop("Networks must be lists or matrices")
  }
  
  if(!method%in%c("ols","kabsch","L1","orig")){
    stop("Invalid method.  Must be one of 'ols', 'kabsch', 'L1','orig'")
  }
  if (method == "kabsch"){
    tf.trans.matrix <- kabsch(net1,net2)
  }
  if (method == "orig"){
    svd.net1 <- svd(net1)
    tf.trans.matrix <- svd.net1$v %*% diag(1/svd.net1$d) %*% t(svd.net1$u) %*% net2
  }
  if (method == "ols"){
    net2.star <- vapply(seq_len(ncol(net1)), function(i,x,y){
      return(lm(y[,i]~x[,i])$resid)
    }, x=net1, y=net2, FUN.VALUE = numeric(dim(net1)[1]))
    tf.trans.matrix <- ginv(t(net1)%*%net1)%*%t(net1)%*%net2.star
    colnames(tf.trans.matrix) <- colnames(net1)
    rownames(tf.trans.matrix) <- colnames(net1)
    if(logging == TRUE){
      print("Using OLS method")
    }
  }
  if (method == "L1"){
    if (!requireNamespace("penalized", quietly = TRUE))
      stop("Package 'penalized' is required for method='L1'. Please install it.")
    net2.star <- vapply(seq_len(ncol(net1)), function(i,x,y){
      lm(y[,i]~x[,i])$resid
    }, x=net1, y=net2, FUN.VALUE = numeric(dim(net1)[1]))
    tf.trans.matrix <- do.call(rbind, lapply(seq_len(ncol(net1)), function(i){
      coef <- NULL
      utils::capture.output(
        l1 <- penalized::optL1(response = net2.star[,i], penalized = net1, fold=5, minlambda1=1, 
                               maxlambda1=2, model="linear", standardize=TRUE),
        z <- penalized::penalized(response = net2.star[,i], penalized = net1, 
                                  model="linear", standardize=TRUE, lambda1 = l1$lambda,
                                  trace = FALSE),
        coef <- penalized::coefficients(z, "penalized"),
        file = NULL
      )
      return(coef)
    }))
    colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
    if(logging == TRUE){
      print("Using L1 method")
    }
  }
  if (standardize){
    tf.trans.matrix <- apply(tf.trans.matrix, 1, function(x){
      x/sum(abs(x))
    })
  }
  
  if (remove.diagonal){
    diag(tf.trans.matrix) <- 0
  }
  colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
  return(tf.trans.matrix)
}

kabsch <- function(P,Q){
  
  P <- apply(P,2,function(x){
    x - mean(x)
  })
  Q <- apply(Q,2,function(x){
    x - mean(x)
  })
  covmat <- cov(P,Q)
  P.bar <- colMeans(P)
  Q.bar <- colMeans(Q)
  num.TFs <- ncol(P)        #n
  num.genes <- nrow(P)    #m
  
  #     covmat <- (t(P)%*%Q - P.bar%*%t(Q.bar)*(num.genes))
  
  svd.res <- svd(covmat-num.TFs*Q.bar%*%t(P.bar))
  
  # Note the scalar multiplier in the middle.
  # NOT A MISTAKE!
  c.k <- colSums(P %*% svd.res$v * Q %*% svd.res$u) - 
    num.genes*(P.bar%*%svd.res$v)*(Q.bar%*%svd.res$u)
  
  E <- diag(c(sign(c.k)))
  
  W <- svd.res$v %*% E %*% t(svd.res$u)
  rownames(W) <- colnames(P)
  colnames(W) <- colnames(P)
  return(W)
}

#' Transformation matrix plot
#'
#' This function plots a hierachically clustered heatmap and 
#' corresponding dendrogram of a transaction matrix
#'
#' @param monsterObj monsterAnalysis Object
#' @param method distance metric for hierarchical clustering.    
#' Default is "Pearson correlation"
#' @export
#' @rawNamespace import(stats, except= c(cov2cor,decompose,toeplitz,lowess,update,spectrum))
#' @return ggplot2 object for transition matrix heatmap
#' @examples
#' # data(yeast)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=10, numMaxCores=1)
#' data(monsterRes)
#' slot(monsterRes, "logging") <- FALSE
#' monsterHclHeatmapPlot(monsterRes)
monsterHclHeatmapPlot <- function(monsterObj, method="pearson"){
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for monsterHclHeatmapPlot. Please install it.")
  if (!requireNamespace("grid", quietly = TRUE))
    stop("Package 'grid' is required for monsterHclHeatmapPlot. Please install it.")
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Package 'reshape2' is required for monsterHclHeatmapPlot. Please install it.")
  if (!requireNamespace("ggdendro", quietly = TRUE))
    stop("Package 'ggdendro' is required for monsterHclHeatmapPlot. Please install it.")
  x <- monsterObj@tm
  if(method=="pearson"){
    dist.func <- function(y) as.dist(cor(y))
  } else {
    dist.func <- dist
  }
  x <- scale(x)
  dd.col <- as.dendrogram(hclust(dist.func(x)))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist.func(t(x))))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- x[col.ord, row.ord]
  xx_names <- attr(xx, "dimnames")
  df <- as.data.frame(xx)
  colnames(df) <- xx_names[[2]]
  df$Var1 <- xx_names[[1]]
  df$Var1 <- with(df, factor(Var1, levels=Var1, ordered=TRUE))
  mdf <- reshape2::melt(df, id.vars = "Var1")
  
  ddata_x <- ggdendro::dendro_data(dd.row)
  ddata_y <- ggdendro::dendro_data(dd.col)
  
  ### Set up a blank theme
  theme_none <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(colour=NA),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  )
  ### Set up a blank theme
  theme_heatmap <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(colour=NA),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  )
  ### Create plot components ###
  # Heatmap
  p1 <- ggplot2::ggplot(mdf, ggplot2::aes(x=variable, y=Var1)) +
    ggplot2::geom_tile(ggplot2::aes(fill=value)) + 
    ggplot2::scale_fill_gradient2() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  # Dendrogram 1
  p2 <- ggplot2::ggplot(ggdendro::segment(ddata_x)) +
    ggplot2::geom_segment(ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_none + ggplot2::theme(axis.title.x=ggplot2::element_blank())
  
  # Dendrogram 2
  p3 <- ggplot2::ggplot(ggdendro::segment(ddata_y)) +
    ggplot2::geom_segment(ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::coord_flip() + theme_none
  
  ### Draw graphic ###
  
  grid::grid.newpage()
  print(p1, vp=grid::viewport(0.80, 0.8, x=0.400, y=0.40))
  print(p2, vp=grid::viewport(0.73, 0.2, x=0.395, y=0.90))
  print(p3, vp=grid::viewport(0.20, 0.8, x=0.910, y=0.43))
}

#' Principal Components plot of transformation matrix
#'
#' This function plots the first two principal components for a 
#' transaction matrix
#'
#' @param monsterObj a monsterAnalysis object resulting from a monster analysis
#' @param title The title of the plot
#' @param clusters A vector indicating the number of clusters to compute
#' @param alpha A vector indicating the level of transparency to be plotted
#' @return ggplot2 object for transition matrix PCA
#' @export
#' @examples
#' # data(yeast)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)#' 
#' data(monsterRes)
#' slot(monsterRes, "logging") <- FALSE
#' # Color the nodes according to cluster membership
#' clusters <- kmeans(monsterGetTm(monsterRes),3)$cluster 
#' monsterTransitionPCAPlot(monsterRes, 
#' title="PCA Plot of Transition - Cell Cycle vs Stress Response", 
#' clusters=clusters)
monsterTransitionPCAPlot <-    function(monsterObj, 
                                         title="PCA Plot of Transition", 
                                         clusters=1, alpha=1){
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for monsterTransitionPCAPlot. Please install it.")
  tm.pca <- princomp(monsterObj@tm)
  odsm <- apply(monsterObj@tm,2,function(x){t(x)%*%x})
  odsm.scaled <- 2*(odsm-mean(odsm))/sd(odsm)+4
  scores.pca <- as.data.frame(tm.pca$scores)
  scores.pca <- cbind(scores.pca,'node.names'=rownames(scores.pca))
  ggplot2::ggplot(data = scores.pca, ggplot2::aes(x = Comp.1, y = Comp.2, label = node.names)) +
    ggplot2::geom_hline(yintercept = 0, colour = "gray65") +
    ggplot2::geom_vline(xintercept = 0, colour = "gray65") +
    ggplot2::geom_text(size = odsm.scaled, alpha=alpha, color=clusters) +
    ggplot2::ggtitle(title)
}

#' This function uses igraph to plot the transition matrix (directed graph) as a network.
#' The edges in the network should be read as A 'positively/negatively contributes to' the 
#' targeting of B in the target state.
#'
#' @param monsterObj monsterAnalysis Object
#' @param numEdges The number of edges to display
#' @param numTopTFs The number of TFs to display, only when rescale='significance'
#' @param rescale string to specify the order of edges. If set to 'significance', 
#' the TFs with the largest dTFI significance (smallest dTFI p-values) will be filtered first before
#' plotting the edges with the largest magnitude in the transition matrix. Otherwise
#' the filtering step will be skipped and the edges with the largest transitions will be plotted.
#' The plotted graph represents the top numEdges edges between the numTopTFs if rescale=='significance'
#' and top numEdges edges otherwise. The edge weight represents the observed transition edges standardized
#' by the null and the node size in the graph is proportional to the p-values of the dTFIs of each
#' TF. When rescale is set to 'significance', the results can be different between two MONSTER runs
#' if the number of permutations is not large enough to sample the null, that is why it is the seed should be set
#' prior to calling MONSTER to get reproducible results. If rescale is set to another value such as 'none', it will
#' produce deterministic results between two identical MONSTER runs.
#' @importFrom igraph graph.data.frame plot.igraph V E V<- E<-
#' @export
#' @return plot the transition matrix (directed graph) as a network.
#' @examples
#' # data(yeast)
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)#' 
#' data(monsterRes)
#' slot(monsterRes, "logging") <- FALSE
#' monsterTransitionNetworkPlot(monsterRes, rescale='significance')
#' monsterTransitionNetworkPlot(monsterRes, rescale='none')

monsterTransitionNetworkPlot <- function(monsterObj, numEdges=100, numTopTFs=10, rescale='significance'){
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Package 'reshape2' is required for monsterTransitionNetworkPlot. Please install it.")
  ## Calculate p-values for off-diagonals
  transitionSigmas <- function(tm.observed, tm.null){
    tm.null.mean <- apply(simplify2array(tm.null), seq_len(2), mean)
    tm.null.sd <- apply(simplify2array(tm.null), seq_len(2), sd)
    sigmas <- (tm.observed - tm.null.mean)/tm.null.sd
  }
  
  tm.sigmas <- transitionSigmas(monsterObj@tm, monsterObj@nullTM)
  diag(tm.sigmas) <- 0
  tm.sigmas.melt <- reshape2::melt(tm.sigmas)
  
  adjMat <- monsterObj@tm
  diag(adjMat) <- 0
  adjMat.melt <- reshape2::melt(adjMat)
  
  adj.combined <- merge(tm.sigmas.melt, adjMat.melt, by=c("Var1","Var2"))
  
  # adj.combined[,1] <- mappings[match(adj.combined[,1], mappings[,1]),2]
  # adj.combined[,2] <- mappings[match(adj.combined[,2], mappings[,1]),2]
  dTFI_pVals_All <- 1-2*abs(.5-monsterCalculateTmPValues(monsterObj, 
                                                             method="z-score"))
  if(rescale=='significance'){
    topTFsIncluded <- names(sort(dTFI_pVals_All)[seq_len(numTopTFs)])
    topTFIndices <- 2>(is.na(match(adj.combined[,1],topTFsIncluded)) + 
                         is.na(match(adj.combined[,2],topTFsIncluded)))
    adj.combined <- adj.combined[topTFIndices,]
  }
  
  adj.combined <- adj.combined[
    abs(adj.combined[,4])>=sort(abs(adj.combined[,4]),decreasing=TRUE)[numEdges],]
  tfNet <- graph_from_data_frame(adj.combined, directed=TRUE)
  vSize <- -log(dTFI_pVals_All)
  vSize[vSize<0] <- 0
  vSize[vSize>3] <- 3
  
  V(tfNet)$size <- vSize[V(tfNet)$name]*5
  V(tfNet)$color <- "yellow"
  E(tfNet)$width <- (abs(E(tfNet)$value.x))*15/max(abs(E(tfNet)$value.x))
  E(tfNet)$color <-ifelse(E(tfNet)$value.x>0, "blue", "red")
  
  plot.igraph(tfNet, edge.arrow.size=2, vertex.label.cex= 1.5, vertex.label.color= "black",main="")
}

#' This function plots the Off diagonal mass of an 
#' observed Transition Matrix compared to a set of null TMs
#'
#' @param monsterObj monsterAnalysis Object
#' @param rescale string indicating whether to reorder transcription
#' factors according to their statistical significance and to 
#' rescale the values observed to be standardized by the null
#' distribution ('significance'), to reorder transcription
#' factors according to the largest dTFIs ('magnitude') with the TF x axis labels proportional to their significance
#' , or finally without ordering them ('none'). When rescale is set to 'significance', 
#' the results can be different between two MONSTER runs if the number of permutations is not large enough to sample 
#' the null, that is why it is the seed should be set prior to calling MONSTER to get reproducible results. 
#' If rescale is set to another value such as 'magnitude' or 'none', it will produce deterministic results 
#' between two identical MONSTER runs.
#' @param plot.title String specifying the plot title
#' @param highlight.tfs vector specifying a set of transcription 
#' factors to highlight in the plot
#' @param nTFs number of TFs to plot in x axis. -1 takes all TFs.
#' @return ggplot2 object for transition matrix comparing observed 
#' distribution to that estimated under the null 
#' @export
#' @examples
#' # data(yeast)
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)#' 
#' data(monsterRes)
#' slot(monsterRes, "logging") <- FALSE
#' monsterdTFIPlot(monsterRes)
monsterdTFIPlot <- function(monsterObj, rescale='none', plot.title=NA, highlight.tfs=NA,
                             nTFs=-1){
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Package 'reshape2' is required for monsterdTFIPlot. Please install it.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for monsterdTFIPlot. Please install it.")
  if(is.na(plot.title)){
    plot.title <- "Differential TF Involvement"
  }
  num.iterations <- length(monsterObj@nullTM)
  # Calculate the off-diagonal squared mass for each transition matrix
  e = monsterCalculateTmStats(monsterObj)
  p.values = e$p.values
  t.values = e$t.values
  ssodm = e$ssodm
  null.ssodm.matrix = e$null.ssodm.matrix
  
  # Process the data for ggplot2
  combined.mat <- cbind(null.ssodm.matrix, ssodm)
  colnames(combined.mat) <- c(rep('Null',num.iterations),"Observed")
  
  
  if (rescale == 'significance'){
    combined.mat <- t(apply(combined.mat,1,function(x){
      (x-mean(x[-(num.iterations+1)]))/sd(x[-(num.iterations+1)])
    }))
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-t.values)]
    x.axis.size    <- 10 # pmin(15,7-log(p.values[order(p.values)]))
  } else if (rescale == 'none'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])
    x.axis.size    <- pmin(15,7-log(p.values))
  } else if (rescale == 'magnitude'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-combined.mat[, dim(combined.mat)[2]])]
    x.axis.size    <- pmin(15,7-log(p.values))
  }
  if(nTFs==-1){
    nTFs = length(x.axis.order)
  }
  null.SSODM.melt <- reshape2::melt(combined.mat)[,-1][,c(2,1)]
  null.SSODM.melt$TF<-rep(rownames(monsterObj@nullTM[[1]]),num.iterations+1)
  
  ## Plot the data
  ggplot2::ggplot(null.SSODM.melt, ggplot2::aes(x=TF, y=value))+
    ggplot2::geom_point(ggplot2::aes(color=factor(Var2), alpha = .5*as.numeric(factor(Var2))), size=2) +
    ggplot2::scale_color_manual(values = c("blue", "red")) +
    ggplot2::scale_alpha(guide = "none") +
    ggplot2::scale_x_discrete(limits = x.axis.order[seq_len(nTFs)] ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title=ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(
                                     angle = 90, hjust = 1, 
                                     size=x.axis.size,face="bold")) +
    ggplot2::ylab("dTFI") +
    ggplot2::ggtitle(plot.title)
  
}


#' Calculate statistics for a tranformation matrix
#'
#' This function powers both the p-value and t-value calculations
#' for a transformation matrix.  It calculates the off-diagonal squared mass
#' for each transition matrix, and then calculates the p-values and t-values
#' for the observed transition matrix compared to the null transition matrices.
#' It is used by the monsterCalculateTmPValues function and by the monsterdTFIPlot
#'
#' @param monsterObj monsterAnalysis Object
#' @param method one of 'z-score' or 'non-parametric'
#' @return p-values, t-values, off-diagonal squared mass for the observed transition matrix,
#' and the null off-diagonal squared mass matrix
#' @export
#' @examples
#' # data(yeast)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)
#' data(monsterRes)
#' slot(monsterRes, "logging") <- FALSE
#' monsterCalculateTmStats(monsterRes)

monsterCalculateTmStats <- function(monsterObj, method="z-score"){
  num.iterations <- length(monsterObj@nullTM)
  # Calculate the off-diagonal squared mass for each transition matrix
  null.SSODM <- lapply(monsterObj@nullTM,function(x){
    apply(x,2,function(y){t(y)%*%y})
  })
  null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
  null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))
  
  ssodm <- apply(monsterObj@tm,2,function(x){t(x)%*%x})
  
  # Get p-value (rank of observed within null ssodm)
  if(method=="non-parametric"){
    seqssodm <- seq_along(ssodm)
    names(seqssodm) <- names(ssodm)

    p.values <- vapply(seqssodm,function(i){
      1-findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
    }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)

    # tvalues are now the rank of the observed ssodm within the null ssodm matrix
    # divided by the number of iterations
    # This is the proportion of null ssodm that are less than the observed ssodm
    t.values <- vapply(seqssodm,function(i){
      findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
    }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
  } else if (method=="z-score"){
    seqssdom=seq_along(ssodm)
    names(seqssdom)=names(ssodm)

    p.values <- 1-pnorm(vapply(seqssdom,function(i){
      (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
    }, FUN.VALUE = numeric(1), USE.NAMES = TRUE))
    
    t.values <- vapply(seqssdom,function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  }, FUN.VALUE = numeric(1), USE.NAMES = TRUE)
  } else {
    stop('Undefined method')
  }
  return(list(p.values=p.values, t.values=t.values, ssodm=ssodm, null.ssodm.matrix=null.ssodm.matrix))
}

#' Calculate statistics for a transformation matrix at the matrix cell level (TF_a -> TF_b)
#'
#' This function powers both the p-value and t-value (or z-score) calculations
#' for each cell in a transformation matrix.
#'
#' @param monsterObj monsterAnalysis Object
#' @param method one of 'z-score' or 'non-parametric'
#' @return p-values, adjusted p-values, z-scores (if applicable)
#' @export
monsterCalculateTmStatsPerCell <- function(monsterObj, method="z-score"){
  
  # Check object input.
  if(!is(monsterObj, "monsterAnalysis")){
    stop("Input must be an object of class 'monsterAnalysis'.")
  }
  if(!(method %in% c("z-score", "non-parametric"))){
    stop(paste("Valid methods include 'z-score' and 'non-parametric'. Invalid method:", method))
  }
  
  # Melt each matrix.
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Package 'reshape2' is required for this method. Please install it.")
  meltedTrans <- reshape2::melt(monsterObj@tm)
  colnames(meltedTrans) <- c("Source", "Target", "Score")

  meltedNullTrans <- do.call(cbind, lapply(1:length(monsterObj@nullTM), function(i){
    
    # Melt the null matrix.
    nullMat <- monsterObj@nullTM[[i]]
    melted <- reshape2::melt(nullMat)
    colnames(melted) <- c("Source", "Target", paste0("Score", i))
    
    # Only keep the TF labels for the first null.
    retval <- as.data.frame(melted[,3])
    colnames(retval) <- paste0("Score", i)
    if(i == 1){
      retval <- melted
    }
    return(retval)
  }))

  # Initialize values to return.
  p.values <- NULL
  p.adj <- NULL
  z.scores <- NULL
  
  # Get p-value (rank of observed within null ssodm)
  if(method=="non-parametric"){
    
    # Calculate n.
    n <- ncol(meltedNullTrans) - 2
    
    # Calculate r.
    repObserved <- matrix(rep(meltedTrans$Score, n), ncol = n)
    nullTransNum <- as.matrix(meltedNullTrans[,3:ncol(meltedNullTrans)])
    isGreater <- nullTransNum > repObserved
    isGreaterBin <- isGreater
    isGreaterBin[which(isGreater == TRUE)] <- 1
    isGreaterBin[which(isGreater == FALSE)] <- 0
    r <- rowSums(isGreaterBin)

    # Calculate p-values.
    p <- (r + 1) / (n + 1)
    pDF <- data.frame(Source = meltedTrans$Source, Target = meltedTrans$Target, Score = p)
    p.values <- reshape2::acast(pDF, Source ~ Target, value.var = "Score")
    
  } else if (method=="z-score"){
    
    # Calculate z-scores.
    nullTransNum <- as.matrix(meltedNullTrans[,3:ncol(meltedNullTrans)])
    nullMeans <- rowMeans(nullTransNum)
    nullSd <- apply(nullTransNum, 1, sd)
    zDF <- data.frame(Source = meltedTrans$Source,
                      Target = meltedTrans$Target,
                      Score = (meltedTrans$Score - nullMeans) / nullSd)
    z.scores <- reshape2::acast(zDF, Source ~ Target, value.var = "Score")
    
    # Calculate p-values using a one-tailed test (we are interested in
    # transitions greater than expected).
    p.values <- pnorm(z.scores, lower.tail=FALSE)
  } else {
    stop('Undefined method')
  }
  
  # Adjust p-values using FDR.
  p.adj <- matrix(
    p.adjust(as.vector(p.values), method = "fdr"),
    nrow = nrow(p.values),
    ncol = ncol(p.values),
    dimnames = dimnames(p.values)
  )
  
  # Return.
  return(list(p.values=p.values, p.adj = p.adj, z.scores=z.scores))
}


#' Calculate p-values for a tranformation matrix
#'
#' This function calculates the significance of an observed
#' transition matrix given a set of null transition matrices
#'
#' @param monsterObj monsterAnalysis Object
#' @param method one of 'z-score' or 'non-parametric'
#' @return vector of p-values for each transcription factor
#' @export
#' @examples
#' # data(yeast)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)
#' data(monsterRes)
#' slot(monsterRes, "logging") <- FALSE
#' monsterCalculateTmPValues(monsterRes)
monsterCalculateTmPValues <- function(monsterObj, method="z-score"){

  e = monsterCalculateTmStats(monsterObj, method = method)
  p.values = e$p.values

  return(p.values)
}

globalVariables(c("Var1", "Var2","value","variable","xend","yend","y","Comp.1", "Comp.2","node.names","TF","i"))

#' Bipartite Edge Reconstruction from Expression data
#'
#' This function generates a complete bipartite network from 
#' gene expression data and sequence motif data
#' 
#' @param motif.data A motif dataset, a data.frame, matrix or exprSet containing 
#' 3 columns. Each row describes an motif associated with a transcription 
#' factor (column 1) a gene (column 2) and a score (column 3) for the motif.
#' @param expr.data An expression dataset, as a genes (rows) by samples (columns)
#' @param verbose logical to indicate printing of output for algorithm progress.
#' @param method String to indicate algorithm method.  Must be one of 
#' "BERE","pearson",or "lda". Default is "BERE".
#' Important note: the direct regulatory network observed from gene expression is currently
#' implemented as a regular correlation as opposed to the partial correlation described 
#' in the paper (please see Schlauch et al., 2017,  https://doi.org/10.1186/s12918-017-0517-y)
#' @param ni.coefficient.cutoff numeric to specify a p-value cutoff at the network
#' inference step.  Default is NA, indicating inclusion of all coefficients.
#' @param alphaw A weight parameter between 0 and 1 specifying proportion of weight 
#' to give to indirect compared to direct evidence. The default is 0.5 to give an 
#' equal weight to direct and indirect evidence.
#' @param randomize logical indicating randomization by genes, within genes or none
#' @param score String to indicate whether motif information will be 
#' readded upon completion of the algorithm
#' to give to indirect compared to direct evidence.  See documentation.
#' @param regularization String parameter indicating one of "none", "L1", "L2"
#' @param cpp logical use C++ for maximum speed, set to false if unable to run.
#' @param logging Whether or not to log progress (default = TRUE)
#' @export
#' @return matrix for inferred network between TFs and genes
#' @examples
#' data(yeast)
#' cc.net <- monsterMonsterNI(yeast$motif,yeast$exp.cc)

monsterMonsterNI <- function(motif.data, 
                              expr.data,
                              verbose=FALSE,
                              randomize="none",
                              method="BERE",
                              ni.coefficient.cutoff=NA,
                              alphaw=1.0,
                              regularization="none",
                              score="motifincluded",
                              cpp=FALSE,
                              logging = TRUE){

  if (method=="BERE"){
    result <- monsterBereFull(motif.data = motif.data, expr.data = expr.data, alpha=alphaw)
  }else{
    if(verbose)
      print('Initializing and validating')
    # Create vectors for TF names and Gene names from Motif dataset
    tf.names   <- sort(unique(motif.data[,1]))
    num.TFs    <- length(tf.names)
    if (is.null(expr.data)){
      stop("Expression data null")
    } else {
      # Use the motif data AND the expr data (if provided) for the gene list
      gene.names <- sort(intersect(motif.data[,2],rownames(expr.data)))
      num.genes  <- length(gene.names)
      
      # Filter out the expr genes without motif data
      expr.data <- expr.data[rownames(expr.data) %in% gene.names,]
      
      # Keep everything sorted alphabetically
      expr.data      <- expr.data[order(rownames(expr.data)),]
      num.conditions <- ncol(expr.data);
      if (randomize=='within.gene'){
        expr.data <- t(apply(expr.data, 1, sample))
        if(verbose)
          print("Randomizing by reordering each gene's expression")
      } else if (randomize=='by.genes'){
        rownames(expr.data) <- sample(rownames(expr.data))
        expr.data           <- expr.data[order(rownames(expr.data)),]
        if(verbose)
          print("Randomizing by reordering each gene labels")
      }
    }
    
    # Bad data checking
    if (num.genes==0){
      stop("Validating data.  No matched genes.\n
            Please ensure that gene names in expression 
            file match gene names in motif file.")
    }
    
    strt<-Sys.time()
    if(num.conditions==0) {
      stop("Number of samples = 0")
      gene.coreg <- diag(num.genes)
    } else if(num.conditions<3) {
      stop('Not enough expression conditions detected to calculate correlation.')
    } else {
      if(verbose)
        print('Verified adequate samples, calculating correlation matrix')
      if(cpp){
        # C++ implementation
        gene.coreg <- rcpp_ccorr(t(apply(expr.data, 1, function(x)(x-mean(x))/(sd(x)))))
        rownames(gene.coreg)<- rownames(expr.data)
        colnames(gene.coreg)<- rownames(expr.data)
        
      } else if(method != "pearson") {
        # Standard r correlation calculation
        gene.coreg <- cor(t(expr.data), method="pearson", use="pairwise.complete.obs")
      }
    }
    if(logging == TRUE){
      print(Sys.time()-strt)
    }
    
    if(verbose)
      print('More data cleaning')
    # Convert 3 column format to matrix format
    colnames(motif.data) <- c('TF','GENE','value')
    if (!requireNamespace("tidyr", quietly = TRUE))
      stop("Package 'tidyr' is required for this method. Please install it.")
    regulatory.network <- tidyr::spread(motif.data, GENE, value, fill=0)
    rownames(regulatory.network) <- regulatory.network[,1]
    # sort the TFs (rows), and remove redundant first column
    regulatory.network <- regulatory.network[order(rownames(regulatory.network)),-1]
    # sort the genes (columns)
    regulatory.network <- as.matrix(regulatory.network[,order(colnames(regulatory.network))])
    
    # Filter out any motifs that are not in expr dataset (if given)
    if (!is.null(expr.data)){
      regulatory.network <- regulatory.network[,colnames(regulatory.network) %in% gene.names]
    }
    
    # store initial motif network (alphabetized for rows and columns)
    #   starting.motifs <- regulatory.network
    
    if(verbose)
      print('Main calculation')
    result <- NULL
    ########################################
    if (method=="pearson"){
      tfNames = levels(motif.data$TF)
      result <- t(cor(t(expr.data),t(expr.data[rownames(expr.data)%in%tfNames,]))^2)
    } else {
      strt<-Sys.time()
      # Remove NA correlations
      gene.coreg[is.na(gene.coreg)] <- 0
      correlation.dif <- sweep(regulatory.network,1,rowSums(regulatory.network),`/`)%*%
        gene.coreg - 
        sweep(1-regulatory.network,1,rowSums(1-regulatory.network),`/`)%*%
        gene.coreg
      result <- sweep(correlation.dif, 2, apply(correlation.dif, 2, sd),'/')
      #   regulatory.network <- ifelse(res>quantile(res,1-mean(regulatory.network)),1,0)
      
      if(logging == TRUE){
        print(Sys.time()-strt)
      }
      ########################################
      if(score=="motifincluded"){
        result <- result + max(result)*regulatory.network
      }
    }
  } 
  return(result)
}

#' Bipartite Edge Reconstruction from Expression data 
#' (composite method with direct/indirect)
#'
#' This function generates a complete bipartite network from 
#' gene expression data and sequence motif data. This NI method
#' serves as a default method for inferring bipartite networks
#' in MONSTER.  Running monsterBereFull can generate these networks
#' independently from the larger MONSTER method.
#'
#' @param motif.data A motif dataset, a data.frame, matrix or exprSet 
#' containing 3 columns. Each row describes an motif associated 
#' with a transcription factor (column 1) a gene (column 2) 
#' and a score (column 3) for the motif.
#' @param expr.data An expression dataset, as a genes (rows) by 
#' samples (columns) data.frame
#' @param alpha A weight parameter specifying proportion of weight 
#' to give to indirect compared to direct evidence.  See documentation.
#' @param lambda if using penalized, the lambda parameter in the penalized logistic regression
#' @param score String to indicate whether motif information will 
#' be readded upon completion of the algorithm
#' @export
#' @return An matrix or data.frame
#' @examples
#' data(yeast)
#' monsterRes <- monsterBereFull(yeast$motif, yeast$exp.cc, alpha=.5)
monsterBereFull <- function(motif.data, 
                             expr.data, 
                             alpha=.5, 
                             lambda=10, 
                             score="motifincluded"){
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("Package 'reshape2' is required for monsterBereFull. Please install it.")
  if (!requireNamespace("penalized", quietly = TRUE))
    stop("Package 'penalized' is required for monsterBereFull. Please install it.")
  
  expr.data <- data.frame(expr.data)
  colnames(motif.data) <- c("TF", "GENE", "score")
  tfdcast <- reshape2::dcast(motif.data,TF~GENE,fill=0, value.var = "score")
  rownames(tfdcast) <- tfdcast[,1]
  tfdcast <- tfdcast[,-1]

  expr.data <- expr.data[sort(rownames(expr.data)),]
  tfdcast <- tfdcast[,sort(colnames(tfdcast)),]
  tfNames <- rownames(tfdcast)[rownames(tfdcast) %in% rownames(expr.data)]
  
  ## Filtering
  # filter out the TFs that are not in expression set
  tfdcast <- tfdcast[rownames(tfdcast)%in%tfNames,]

  # Filter out genes that aren't targetted by anything 7/28/15
  commonGenes <- intersect(colnames(tfdcast),rownames(expr.data))
  expr.data <- expr.data[commonGenes,]
  tfdcast <- tfdcast[,commonGenes]
  
  # check that IDs match
  if (prod(rownames(expr.data)==colnames(tfdcast))!=1){
    stop("ID mismatch")
  }
  ## Get direct evidence
  directCor <- t(cor(t(expr.data),t(expr.data[rownames(expr.data)%in%tfNames,]))^2)
  
  ## Get the indirect evidence    
  result <- t(apply(tfdcast, 1, function(x){
    tfTargets <- as.numeric(x)
    
    # Ordinary Logistic Reg
    # z <- glm(tfTargets ~ ., data=expr.data, family="binomial")
    
    # Penalized Logistic Reg
    expr.data[is.na(expr.data)] <- 0
    prediction <- NULL
    utils::capture.output(
      z <-penalized::penalized(response=tfTargets, penalized=expr.data, unpenalized=~0,
                               lambda2=lambda, model="logistic", standardize=TRUE),
      prediction <- penalized::predict(z, expr.data),
      file = NULL
    )
    #z <- optL1(tfTargets, expr.data, minlambda1=25, fold=5)
    return(prediction)
  }))
  
  ## Convert values to ranks

  directCor <- matrix(rank(directCor), ncol=ncol(directCor))
  result <- matrix(rank(result), ncol=ncol(result))
  consensus <- directCor*(1-alpha) + result*alpha
  rownames(consensus) <- rownames(tfdcast)
  colnames(consensus) <- rownames(expr.data)
  consensusRange <- max(consensus)- min(consensus)
  if(score=="motifincluded"){
    consensus <- as.matrix(consensus + consensusRange*tfdcast)
  }
  return(consensus)
}

globalVariables(c("expr.data","lambda","rcpp_ccorr","GENE", "TF","value"))

#' MONSTER results from example cell-cycle yeast transition
#'
#'This data contains the MONSTER result from analysis of Yeast Cell cycle, included in data(yeast).  
#'This result arbitrarily takes the first 20 gene expression samples in yeast$cc to be the baseline condition, and the final 20 samples to be the final condition.
#'
#' @docType data
#' @keywords datasets
#' @name monsterRes
#' @usage data(monsterRes)
#' @format MONSTER obj
#' #' @references Schlauch, Daniel, et al. "Estimating drivers of cell state transitions using gene regulatory network models." BMC systems biology 11.1 (2017): 1-10.
#' 
NULL

#' Toy data derived from three gene expression datasets and a mapping from transcription factors to genes.
#'
#'This data is a list containing gene expression data from three separate yeast studies along with data mapping yeast transcription factors with genes based on the presence of a sequence binding motif for each transcription factor in the vicinity of each gene. 
#' The motif data.frame, yeast$motif, describes a set of pairwise connections where a specific known sequence motif of a transcription factor was found upstream of the corresponding gene.   
#' The expression data, yeast$exp.ko, yeast$exp.cc, and yeast$exp.sr, are three gene expression datasets measured in conditions of gene knockout, cell cycle, and stress response, respectively. 
#' @docType data
#' @keywords datasets
#' @name yeast
#' @usage data(yeast)
#' @format A list containing 4 data.frames
#' @return A list of length 4
#' @references Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions. PLoS One. 2013 May 31;8(5):e64832.
#' 
NULL


#' MONSTER quick-start with pre-made regulatory networks
#'
#' This function is a wrapper to simplify usage of the \code{monster} function in
#' the case where the pair of regulatory networks to be contrasted have already
#' been estimated, either as \code{panda} objects, or represented as an
#' adjacency matrix with regulators in rows and genes in columns.
#'
#' @param exp_graph matrix or PANDA object (generated by \code{netzooR::panda}); if matrix, should be the adjacency matrix for the network with regulators in rows and genes in columns, both named. This should be the network for the experimental (case) group.
#' @param control_graph matrix or PANDA object (generated by \code{netzooR::panda}); if matrix, should be the adjacency matrix for the network with regulators in rows and genes in columns, both named. This should be the network for the control (reference) group.
#' @param nullPerms numeric; defaults to 1000. Number of null permutations to perform. See \code{monster} for more details.
#' @param numMaxCores numeric; defaults to 3. Maximum number of cores to use; will be the minimum of this number and the actual available cores. See \code{monster} for more details.
#' @param logging Whether or not to print progress (default = TRUE)
#' @param ... other arguments for \code{monster} may be passed.
#'
#' @return \code{monster} object
#' @export
#'
#' @examples
#' 
#' \donttest{
#'
#' # Generating PANDA networks for demonstration:
#' # For the purposes of this example, first partition the pandaToyData samples, then perform panda:
#' pandaResult_exp <- panda(pandaToyData$motif, pandaToyData$expression[,1:25], pandaToyData$ppi)
#' pandaResult_control <- panda(pandaToyData$motif, pandaToyData$expression[,26:50], pandaToyData$ppi)
#'
#' # function takes both panda objects and matrices, or a mixture
#' monster_res1 <- domonster(pandaResult_exp, pandaResult_control, numMaxCores = 1)
#' monster_res2 <- domonster(pandaResult_exp@regNet, pandaResult_control@regNet, numMaxCores = 1)
#' monster_res3 <- domonster(pandaResult_exp@regNet, pandaResult_control, numMaxCores = 1)
#' }
domonster <- function(exp_graph, control_graph, nullPerms = 1000, numMaxCores = 3, logging = TRUE,...){
  if('panda' %in% class(exp_graph)){
    exp_graph <- exp_graph@regNet
  }
  if('panda' %in% class(control_graph)){
    control_graph <- control_graph@regNet
  }
  
  if(!identical(colnames(control_graph), colnames(exp_graph))){
    genes_keep <- intersect(colnames(control_graph), colnames(exp_graph))
    
    exp_graph <- exp_graph[,genes_keep]
    control_graph <- control_graph[,genes_keep]
    
    if(logging == TRUE){
      cat('\nNote: column names are not identical in the input; taking the intersection of them:\n')
      cat(paste0(length(genes_keep), ' genes included \n'))
    }
  }
  if(!identical(rownames(control_graph), rownames(exp_graph))){
    reg_keep <- intersect(rownames(control_graph), rownames(exp_graph))
    
    exp_graph <- exp_graph[reg_keep,]
    control_graph <- control_graph[reg_keep,]
    
    if(logging == TRUE){
      cat('\nNote: row names are not identical in the input; taking the intersection of them:\n')
      cat(paste0(length(reg_keep), ' regulators included \n'))
    }
  }
  
  combdf <- as.data.frame(cbind(control_graph, exp_graph))
  combdes <- c(rep(0, ncol(control_graph)), rep(1, ncol(exp_graph)))
  res <- monster(expr = combdf,
                 design = combdes,
                 motif = NA,
                 mode = 'regNet',
                 nullPerms = nullPerms,
                 numMaxCores = numMaxCores,
                 logging = logging,
                 ...)
  return(res)
}
