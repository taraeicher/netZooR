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
#'                               Types can be either "numeric" or "categorical" 
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param compute_cor : Whether or not to compute the correlation matrix. Default is TRUE.
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
#' phenotype_dictionary = c("categorical", "numeric")
#' 
#' pathways = list()
#' pathways$pathway1 = sample(rownames(expression_data), 5)
#' pathways$pathway2 = sample(rownames(expression_data), 3)
#' pathways$pathway1 = sample(rownames(expression_data), 7)
#'
#' # Run seahorse
#' results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways)
#'  
#' @export
seahorse <- function(expression, phenotype, phenotype_dictionary, pathways, compute_cor = TRUE,
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
  
  # Compute coexpression of genes
  results$coexpression <- NA
  if(compute_cor == TRUE){
    results$coexpression = cor(t(expression), use="pairwise.complete.obs")
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

#' Computes correlation for numeric phenotypes and ANOVA for categorical
#' phenotypes.
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
    
    if (phenotype_dictionary[i] == "numeric"){
      output_seahorse = gsea_numeric(expression, pheno, pathways, method = method)
      output_seahorse_padj <- rep("NA", length(output_seahorse$cor))
    }else {
      output_seahorse = gsea_categorical(expression, pheno, pathways)
      output_seahorse_padj <- stats::p.adjust(output_seahorse$cor, method = pval_adj_method)
    }
    phenoAssoc[[pheno_name]] = data.frame(stat = output_seahorse$cor,
                                          padj = output_seahorse_padj,
                                          row.names = names(output_seahorse$cor))
    gsea[[pheno_name]] = output_seahorse$GSEA
  }
  return(list(pheno = phenoAssoc, gsea = gsea))
}

#' Function to run GSEA for a numeric phenotype
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
gsea_numeric <- function(expression, pheno, pathways, method){
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

#' Function to run GSEA for a categorical phenotype
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
gsea_categorical <- function(expression, pheno, pathways){
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  phenotype_vector = factor(as.character(pheno))
  cor = unlist(apply(expression, MARGIN=1, function(x){anova(lm(as.numeric(x)~phenotype_vector))$`Pr(>F)`[1]}))
  output_seahorse$cor = cor
  
  # Run GSEA
  cor_rank = sort(cor, decreasing = T)
  fgseaRes <- fgsea::fgsea(pathways, cor_rank, minSize=15, maxSize=500, scoreType = "pos")
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}