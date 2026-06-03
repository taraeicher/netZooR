#' Run SEAHORSE in R
#' Author(s): Enakshi Saha
#' Description:
#'               SEAHORSE computes gene-gene coexpression matrix,
#'               associations between a set of given phenotypes and each gene,
#'               performs gene set enrichment analysis (GSEA) for each phenotype,
#'               using the measures of association.
#'               GSEA is performed through R package "fgsea".
#'              
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be "dichotomous", "nominal", or "continuous" 
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param compute_gene_cor : Whether or not to compute the gene-gene correlation matrix. Default is TRUE.
#' @param compute_phenotype_cor : Whether or not to compute the phenotype-phenotype association matrix. Default is TRUE.
#' @param assoc_method : The method used to infer associations between phenotypes and genes. Default
#' is "pearson". Other options are "spearman", "kendall", and "linear". The "pearson", "kendall", and "spearman" options
#' compute correlations between each phenotype and gene independently for numeric phenotypes and
#' compute an ANOVA for categorical phenotypes. linear" computes a linear regression
#' model of the form gene ~ phenotype1 + phenotype2 + ... + phenotypeN and returns the p-values
#' for each coefficient.
#' @param transform_expression : Whether or not to transform gene expression data into logCPM.
#' This parameter is only used for the "linear" association method. Use if your data are raw RNA-seq
#' counts. Default is TRUE.
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' Outputs:
#' @return results    : a list containing three objects
#'         results$coexpression: a gene x gene Pearson correlation matrix.
#'         results$phenotype_association : a list containing a vector for each phenotype
#'         results$GSEA: a list containing a matrix of GSEA results for each phenotype
#'
#' @examples
#'
#' expression_data = data.frame(matrix(rexp(200, rate=.1), ncol=10, nrow = 20))
#' rownames(expression_data) = paste("gene", 1:20, sep = "")
#' colnames(expression_data) = paste("sample", 1:10, sep = "")
#' 
#' phenotype_data = data.frame(matrix(0, ncol=2, nrow = 10))
#' colnames(phenotype_data) = c("sex", "height")
#' rownames(phenotype_data) = colnames(expression_data)
#' phenotype_data$sex = c(rep("male", nrow(phenotype_data)/2), rep("female", nrow(phenotype_data)/2))
#' phenotype_data$height = 65 + sample.int(10, nrow(phenotype_data), replace = TRUE)
#' 
#' phenotype_dictionary = c("dichotomous", "continuous")
#' 
#' pathways = list()
#' pathways$pathway1 = sample(rownames(expression_data), 5)
#' pathways$pathway2 = sample(rownames(expression_data), 3)
#' pathways$pathway3 = sample(rownames(expression_data), 7)
#'
#' # Run seahorse
#' results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways)
#'  
#' @export
seahorse <- function(expression, phenotype, phenotype_dictionary, pathways, compute_gene_cor = TRUE,
                     compute_phenotype_cor = TRUE,
                     assoc_method = "pearson", transform_expression = TRUE, pval_adj_method = "none"){
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Package 'stats' is required but not installed.")
  }
  set.seed(0)
  
  # Return an error if pval_adj_method is not an accepted method.
  if(!(pval_adj_method %in% stats::p.adjust.methods)){
    stop(paste(pval_adj_method, "is not a valid method for stats::p.adjust()."))
  }
  
  # Check that association method is valid.
  if(!assoc_method %in% c("pearson", "spearman", "kendall", "linear")){
    stop(paste(assoc_method, "is an invalid association method!"))
  }else if (assoc_method == "linear" && !requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for linear regression. Install it with BiocManager::install('limma')")
  }
  results = list()
  
  # Ensure that the column names are valid for the phenotypic data.
  colnames(phenotype) <- make.names(colnames(phenotype))
  
  # Compute coexpression of genes
  results$coexpression <- NA
  if(compute_gene_cor == TRUE){
    results$coexpression = cor(t(expression), use="pairwise.complete.obs")
  }
  
  # Compute coexpression of phenotypes
  results$phenocor <- NA
  if(compute_phenotype_cor == TRUE){
    results$phenocor = computePhenotypeCorrelations(phenotype = phenotype,
                                                    phenotype_dictionary = phenotype_dictionary,
                                                    method = assoc_method, pval_adj_method = pval_adj_method)
  }
  
  # Compute association of gene expression with phenotypes and run GSEA
  results$phenotype_association = list()
  results$GSEA = list()
  
  if(assoc_method %in% c("pearson", "spearman", "kendall")){
    corr <- computeCorrelations(expression = expression, phenotype = phenotype,
                                   pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                   method = assoc_method, pval_adj_method = pval_adj_method)
    results$phenotype_association <- corr$pheno
    results$GSEA <- corr$gsea
  }else{
    linReg <- computeLinearRegression(expression = expression, phenotype = phenotype,
                                       pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                       transform_expression = transform_expression,
                                      pval_adj_method = pval_adj_method)
    results$phenotype_association <- linReg$pheno
    results$GSEA <- linReg$gsea
  }
  return(results)
}

#' Computes linear regression for both numeric and categorical phenotypes. Returns p-values
#' computed from moderated t-statistics using ebayes() and uses t-statistics as input to the
#' GSEA analysis.
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param transform_expression : Whether or not to transform gene expression data into logCPM.
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' This parameter is only used for the "linear" association method. Use if your data are raw RNA-seq
#' counts. Default is TRUE.
computeLinearRegression <- function(expression, phenotype, phenotype_dictionary, pathways,
                                    transform_expression, pval_adj_method){
  
  # Initialize output.
  phenoAssoc <- list
  gsea <- list

  # Compute a design matrix for the phenotypic data.
  formula <- paste("~ ", paste(colnames(phenotype)[1:(ncol(phenotype) - 1)], 
                                   collapse = " + "), colnames(phenotype)[ncol(phenotype)],
                   sep = " + ")
  design <- model.matrix(object = stats::as.formula(formula), data = phenotype)
  
  # Check whether samples were removed and modify gene expression data accordingly.
  if(length(setdiff(colnames(expression), rownames(design)) > 0)){
    expressionSampsOld <- ncol(expression)
    expression <- expression[,rownames(design)]
    expressionSampsNew <- ncol(expression)
    message(paste("Out of", expressionSampsOld, "we retained", expressionSampsNew,
                  "samples."))
  }

  # Transform the gene expression data using voom.
  if(transform_expression == TRUE){
    expression <- limma::voom(counts = expression, design = design)
  }

  # Run the linear models.
  fit <- limma::lmFit(object = expression, design = design)

  # Compute empirical Bayes statistics.
  bayes <- limma::eBayes(fit = fit)

  # Extract p-values and t-statistics.
  t_stats <- bayes$t
  p_vals  <- bayes$p.value
  
  # Format p-values as a list for all covariates.
  p_list <- lapply(colnames(p_vals), function(c){
    p <- p_vals[,c]
    padj <- stats::p.adjust(p_vals[,c], method = pval_adj_method)
    return(data.frame(stat = p, padj = padj, row.names = rownames(p_vals)))
  })
  names(p_list) <- colnames(p_vals)

  phenoAssoc = p_list

  # Run GSEA
  gsea_list <- lapply(colnames(t_stats), function(c){
    tscore <- t_stats[,c]
    t_rank = sort(tscore, decreasing = T)
    fgseaRes <- fgsea::fgsea(pathways = pathways, stats = t_rank, minSize=15, maxSize=500)
    return(fgseaRes)
  })
  names(gsea_list) <- colnames(t_stats)
  
  gsea = gsea_list
  return(list(pheno = phenoAssoc, gsea = gsea))
}

#' Computes gene-phenotype correlation for continuous phenotypes, ANOVA for nominal
#' phenotypes, and t-test for dichotomous phenotypes
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param method : One of "pearson", "spearman", or "kendall".
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
computeCorrelations <- function(expression, phenotype, phenotype_dictionary, pathways, method,
                                pval_adj_method){
  phenoAssoc <- list()
  gsea <- list()
  for (i in 1:ncol(phenotype)){
    pheno = phenotype[,i]
    pheno_name = colnames(phenotype)[i]
    
    if (phenotype_dictionary[i] == "continuous"){
      output_seahorse = gsea_continuous(expression, pheno, pathways, method = method)
      output_seahorse_padj <- rep("NA", length(output_seahorse$cor))
    }else if (phenotype_dictionary[i] == "nominal") {
      output_seahorse = gsea_nominal(expression, pheno, pathways)
      output_seahorse_padj <- stats::p.adjust(output_seahorse$cor, method = pval_adj_method)
    }else{
    	output_seahorse = gsea_dichotomous(expression, pheno, pathways)
    	output_seahorse_padj <- stats::p.adjust(output_seahorse$cor, method = pval_adj_method)
    }
    phenoAssoc[[pheno_name]] = data.frame(stat = output_seahorse$cor,
                                          padj = output_seahorse_padj,
                                          row.names = names(output_seahorse$cor))
    gsea[[pheno_name]] = output_seahorse$GSEA
  }
  return(list(pheno = phenoAssoc, gsea = gsea))
}

#' Function to run GSEA for a continuous phenotype
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param method : One of "pearson", "spearman", or "kendall".
#' @export
gsea_continuous <- function(expression, pheno, pathways, method){
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  phenotype_vector = as.numeric(pheno)
  cor = unlist(apply(expression, MARGIN=1, function(x){cor(as.numeric(x), phenotype_vector, use="pairwise.complete.obs",
                                                           method = method)}))
  output_seahorse$cor = cor
  
  # Run GSEA
  cor_rank = sort(cor, decreasing = T)
  fgseaRes <- fgsea::fgsea(pathways, cor_rank, minSize=15, maxSize=500)
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}

#' Function to run GSEA for a nominal phenotype
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @export
gsea_nominal <- function(expression, pheno, pathways){
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  phenotype_vector = factor(as.character(pheno))
  phenotypes <- unique(phenotype_vector)
  if(length(phenotypes) <= 2){
    stop("Phenotype set to nominal but has 2 levels or fewer.")
  }
  cor = unlist(apply(expression, MARGIN=1, function(x){anova(lm(as.numeric(x)~phenotype_vector))$`Pr(>F)`[1]}))
  output_seahorse$cor = cor
  
  # Run GSEA
  fgseaRes <- fgsea::fgsea(pathways, -1 * log10(cor), minSize=15, maxSize=500, scoreType = "pos")
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}

#' Function to run GSEA for a dichotomous phenotype
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @export
gsea_dichotomous <- function(expression, pheno, pathways){
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  if (!requireNamespace("matrixTests", quietly = TRUE)) {
    stop("Package 'matrixTests' is required but not installed.")
  }
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  phenotype_vector = factor(as.character(pheno))
  phenotypes <- unique(phenotype_vector)
  if(length(phenotypes) > 2){
    stop("Phenotype set to dichotomous but has > 2 levels.")
  }
  group1 <- expression[,which(phenotype_vector == phenotypes[1])]
  group2 <- expression[,which(phenotype_vector == phenotypes[2])]
  tres <- matrixTests::row_t_welch(group1, group2)
  cor = tres$pvalue
  names(cor) <- rownames(tres)
  output_seahorse$cor = cor
  
  # Run GSEA
  fgseaRes <- fgsea::fgsea(pathways, -1 * log10(cor), minSize=15, maxSize=500, scoreType = "pos")
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}

#' Computes phenotype-phenotype correlation for continuous phenotypes, ANOVA for nominal
#' phenotypes, and t-test for dichotomous phenotypes
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param method : One of "pearson", "spearman", or "kendall".
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
computePhenotypeCorrelations <- function(phenotype, phenotype_dictionary, method,
                                pval_adj_method){
  
  # For each phenotype, compute all of its associations.
  phenoAssoc <- lapply(1:(ncol(phenotype)-1), function(i){
    # Extract the phenotype.
    pheno = phenotype[,i]
    pheno_name = colnames(phenotype)[i]
    phenoType <- phenotype_dictionary[i]
    
    # For dichotomous phenotypes, calculate Chi-squared + Cramer's V OR
    # Fisher-Freeman-Halton Test for all categorical comparisons and t-test for
    # all continuous comparisons.
    #
    # For nominal phenotypes, calculate calculate Chi-squared + Cramer's V OR
    # Fisher-Freeman-Halton Test for all categorical comparisons and ANOVA for
    # all continuous comparisons.
    #
    # For continuous phenotypes, calculate t-test for all dichotomous comparisons,
    # ANOVA for all nominal comparisons, and correlation for all continuous comparisons.
    categorical <- c("dichotomous", "nominal")
    output_seahorse <- NULL
    
    # Only compute upper triangular matrix.
    rangeAfter <- (i+1):ncol(phenotype)
    phenosAfter <- as.data.frame(phenotype[,rangeAfter])
    colnames(phenosAfter) <- colnames(phenotype)[rangeAfter]
    typesAfter <- phenotype_dictionary[rangeAfter]
    
    # Run associations.
    if(phenoType == "dichotomous"){
      
      # Do categorical phenotypes.
      whichCat <- which(typesAfter %in% categorical)
      phenoCat <- as.data.frame(phenosAfter[,whichCat])
      colnames(phenoCat) <- colnames(phenosAfter)[whichCat]
      rownames(phenoCat) <- rownames(phenosAfter)
      output_seahorse_cat <- phenotype_chisq(phenotype = pheno, 
                                                  phenotypesToCompare = phenoCat, 
                                                  phenotypeType = phenoType, 
                                                  phenotypesToCompareType = typesAfter[whichCat])
      output_seahorse_padj_cat <- stats::p.adjust(output_seahorse_cat$cor, method = pval_adj_method)
      
      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- phenotype_ttest(phenotype = pheno, 
                                            phenotypesToCompare = phenoCon, 
                                            phenotypeType = phenoType, 
                                            phenotypesToCompareType = typesAfter[whichContinuous])
      output_seahorse_padj_con <- stats::p.adjust(output_seahorse_con$cor, method = pval_adj_method)
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_cat$cor, output_seahorse_con$cor),
                               cramerV = c(output_seahorse_cat$V, output_seahorse_con$V),
                               padj = c(output_seahorse_padj_cat, output_seahorse_padj_con),
                               row.names = c(rownames(output_seahorse_cat), rownames(output_seahorse_con)))
      
    }else if(phenoType == "nominal"){
      
      # Do categorical phenotypes.
      whichCat <- which(typesAfter %in% categorical)
      phenoCat <- as.data.frame(phenosAfter[,whichCat])
      colnames(phenoCat) <- colnames(phenosAfter)[whichCat]
      rownames(phenoCat) <- rownames(phenosAfter)
      output_seahorse_cat <- phenotype_chisq(phenotype = pheno, 
                                                  phenotypesToCompare = phenoCat, 
                                                  phenotypeType = phenoType, 
                                                  phenotypesToCompareType = typesAfter[whichCat])
      output_seahorse_padj_cat <- stats::p.adjust(output_seahorse_cat$cor, method = pval_adj_method)
      
      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- phenotype_anova(phenotype = pheno, 
                                            phenotypesToCompare = phenoCon, 
                                            phenotypeType = phenoType, 
                                            phenotypesToCompareType = typesAfter[whichContinuous])
      output_seahorse_padj_con <- stats::p.adjust(output_seahorse_con$cor, method = pval_adj_method)
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_cat$cor, output_seahorse_con$cor),
                                   cramerV = c(output_seahorse_cat$V, output_seahorse_con$V),
                                   padj = c(output_seahorse_padj_cat, output_seahorse_padj_con),
                                   row.names = c(rownames(output_seahorse_cat), rownames(output_seahorse_con)))
      
    }else{

      # Do dichotomous phenotypes.
      whichDichotomous <- which(typesAfter == "dichotomous")
      phenoDich <- as.data.frame(phenosAfter[,whichDichotomous])
      colnames(phenoDich) <- colnames(phenosAfter)[whichDichotomous]
      rownames(phenoDich) <- rownames(phenosAfter)
      output_seahorse_dich <- phenotype_ttest(phenotype = pheno, 
                                                  phenotypesToCompare = phenoDich, 
                                                  phenotypeType = phenoType, 
                                                  phenotypesToCompareType = typesAfter[whichDichotomous])
      output_seahorse_padj_cat <- stats::p.adjust(output_seahorse_dich$cor, method = pval_adj_method)
      
      # Do nominal phenotypes.
      whichNominal <- which(typesAfter == "nominal")
      phenoNom <- as.data.frame(phenosAfter[,whichNominal])
      colnames(phenoNom) <- colnames(phenosAfter)[whichNominal]
      rownames(phenoNom) <- rownames(phenosAfter)
      output_seahorse_nom <- phenotype_anova(phenotype = pheno, 
                                            phenotypesToCompare = phenoNom, 
                                            phenotypeType = phenoType, 
                                            phenotypesToCompareType = typesAfter[whichNominal])
      output_seahorse_padj_nom <- stats::p.adjust(output_seahorse_nom$cor, method = pval_adj_method)

      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con = phenotype_cor(phenotype = pheno, 
                                            phenotypesToCompare = phenoCon, 
                                          method = method)
      output_seahorse_padj_con <- rep(NA, nrow(output_seahorse_con))
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_dich$cor, output_seahorse_nom$cor, output_seahorse_con$cor),
                                   cramerV = c(output_seahorse_dich$V, output_seahorse_nom$V, output_seahorse_con$V),
                                   padj = c(output_seahorse_padj_cat, output_seahorse_padj_nom, output_seahorse_padj_con),
                                   row.names = c(rownames(output_seahorse_dich), rownames(output_seahorse_nom), rownames(output_seahorse_con)))
      
    }
    # Add the NA values.
    varsNotAdded <- setdiff(colnames(phenotype), row.names(output_seahorse))
    additionalVars <- data.frame(stat = rep(NA, length(varsNotAdded)),
                                 cramerV = rep(NA, length(varsNotAdded)),
                                 padj = rep(NA, length(varsNotAdded)),
                                 row.names = varsNotAdded)
    output_seahorse <- rbind(output_seahorse, additionalVars)
    output_seahorse <- output_seahorse[colnames(phenotype),]
    return(output_seahorse)
  })
  names(phenoAssoc) <- colnames(phenotype)[1:(ncol(phenotype)-1)]
  
  # Compile each vector list into a matrix. Add an extra NA row at the end.
  statMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$stat)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))
  vMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$cramerV)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))
  adjMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$padj)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))

  # Return all matrices.
  phenoAssocList <- list(stat = statMat, cramerV = vMat, padj = adjMat)
  return(phenoAssocList)
}

#' Function to run associations between categorical (dichotomous or nominal) phenotypes.
#' If at least 5 samples exist in each group, run a Chi-square test and compute Cramer's V.
#' Otherwise, run a Fisher-Freeman-Halton Test.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is dichotomous or nominal.
#' @param phenotypesToCompareType A vector of types (dichotomous, nominal) corresponding to all
#' phenotypes against which to compare.
#' @returns A data frame with two vectors: "cor", which lists the Chi-square or Fisher-Freeman-Halton Test p-values, 
#' and "V" which lists the Cramer's V statistics (will be NA if Fisher-Freeman-Halton Test is computed)
phenotype_chisq <- function(phenotype, phenotypesToCompare, phenotypeType, phenotypesToCompareType){
  return(data.frame(cor = rep(NA, length(phenotypesToCompareType)),
                    V = rep(NA, length(phenotypesToCompareType)),
                    row.names = colnames(phenotypesToCompare)))
}

#' Function to compute ANOVA statistics between nominal and continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is nominal or continuous.
#' @param phenotypesToCompareType A vector of types (nominal, continuous) corresponding to all
#' phenotypes against which to compare.
#' @returns A data frame with two vectors: "cor", which lists the ANOVA p-values, 
#' and "V" which is set to NA.
phenotype_anova <- function(phenotype, phenotypesToCompare, phenotypeType, phenotypesToCompareType){
  return(data.frame(cor = rep(NA, length(phenotypesToCompareType)),
                    V = rep(NA, length(phenotypesToCompareType)),
                    row.names = colnames(phenotypesToCompare)))
}

#' Function to compute Welch's t-test statistics between dichotomous and continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is dichotomous or continuous.
#' @param phenotypesToCompareType A vector of types (dichotomous, continuous) corresponding to all
#' phenotypes against which to compare.
#' @returns A data frame with two vectors: "cor", which lists the t-test p-values, 
#' and "V" which is set to NA.
phenotype_ttest <- function(phenotype, phenotypesToCompare, phenotypeType, phenotypesToCompareType){
  return(data.frame(cor = rep(NA, length(phenotypesToCompareType)),
                    V = rep(NA, length(phenotypesToCompareType)),
                    row.names = colnames(phenotypesToCompare)))
}

#' Function to compute correlation between continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' phenotypes against which to compare.
#' @param method : One of "pearson", "spearman", or "kendall".
#' @returns A data frame with two vectors: "cor", which lists the t-test p-values, 
#' and "V" which is set to NA.
phenotype_cor <- function(phenotype, phenotypesToCompare, method){
  return(data.frame(cor = rep(NA, ncol(phenotypesToCompare)),
                    V = rep(NA, ncol(phenotypesToCompare))))
}

#' Function to run phenotype correlations for a pair of dichotomous or nominal phenotypes.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @export
gsea_dichotomous <- function(expression, pheno, pathways){
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  if (!requireNamespace("matrixTests", quietly = TRUE)) {
    stop("Package 'matrixTests' is required but not installed.")
  }
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  phenotype_vector = factor(as.character(pheno))
  phenotypes <- unique(phenotype_vector)
  if(length(phenotypes) > 2){
    stop("Phenotype set to dichotomous but has > 2 levels.")
  }
  group1 <- expression[,which(phenotype_vector == phenotypes[1])]
  group2 <- expression[,which(phenotype_vector == phenotypes[2])]
  tres <- matrixTests::row_t_welch(group1, group2)
  cor = tres$pvalue
  names(cor) <- rownames(tres)
  output_seahorse$cor = cor
  
  # Run GSEA
  fgseaRes <- fgsea::fgsea(pathways, -1 * log10(cor), minSize=15, maxSize=500, scoreType = "pos")
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}
