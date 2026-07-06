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
#'                     We assume the data are adequately transformed for use with limma()
#'                     (i.e. log-scaled and/or transformed using VOOM or DeSeq2).
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be "dichotomous", "nominal", or "continuous" 
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param compute_gene_cor : Whether or not to compute the gene-gene correlation matrix. Default is TRUE.
#' @param compute_phenotype_cor : Whether or not to compute the phenotype-phenotype association matrix. Default is TRUE.
#' @param compute_gene_phenotype_cor: Whether or not to compute the gene-phenotype association matrix. Default is TRUE.
#' @param assoc_method : The method used to infer associations between phenotypes and genes. Default
#' is "spearman". Other options are "pearson", "kendall", and "linear". The "pearson", "kendall", and "spearman" options
#' compute correlations between each phenotype and gene independently for numeric phenotypes and
#' compute an ANOVA for categorical phenotypes. linear" computes a linear regression
#' model of the form gene ~ phenotype1 + phenotype2 + ... + phenotypeN and returns the p-values
#' for each coefficient.
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' @param usage_report_file Location where the usage report should be stored. Must be RDS format.
#' @param verbose Whether or not to print statements notifying user of run status.
#' Outputs:
#' @return results    : a list containing three objects
#'         results$coexpression: a gene x gene correlation matrix.
#'         results$phenotype_association : a list containing a vector for each phenotype
#'         results$GSEA: a list containing a matrix of GSEA results for each phenotype
#'         results$phenocor: an upper-triangular matrix of phenotype-phenotype associations
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
                     compute_phenotype_cor = TRUE, compute_gene_phenotype_cor = TRUE,
                     assoc_method = "spearman", pval_adj_method = "none",
                     usage_report_file = NULL, verbose = FALSE){
  
  # Check packages.
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Package 'stats' is required but not installed.")
  }
  if (!requireNamespace("matrixTests", quietly = TRUE)) {
    stop("Package 'matrixTests' is required but not installed.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE) && !is.null(usage_report_file)) {
    stop(paste("Package 'peakRAM' is required to monitor memory and time usage.",
               "Install it or set usage_report_file = NULL to turn off monitoring."))
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required but not installed.")
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
  
  # Check that file is valid.
  if(!assoc_method %in% c("pearson", "spearman", "kendall", "linear")){
    stop(paste(assoc_method, "is an invalid association method!"))
  }else if (assoc_method == "linear" && !requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for linear regression. Install it with BiocManager::install('limma')")
  }
  results = list()
  
  # Ensure that the column names are valid for the phenotypic data.
  colnames(phenotype) <- make.names(colnames(phenotype))
  
  # Check that phenotype dictionary is correct.
  for(i in 1:length(phenotype_dictionary)){
    dictType <- phenotype_dictionary[i]
    pheno <- phenotype[,i]
    phenoName <- colnames(phenotype)[i]
    
    # Check that nominal have more than 2 levels.
    levels <- unique(pheno[which(!is.na(pheno))])
    if(dictType == "nominal" && length(levels) <= 2){
      stop(paste("Phenotype", phenoName, "set to nominal but has 2 levels or fewer."))
    }
    # Check that dichotomous don't have more than 2 levels.
    if(dictType == "dichotomous" && length(levels) > 2){
      stop(paste("Phenotype", phenoName, "set to dichotomous but has more than 2 levels."))
    }
    # Check that continuous can be converted to a numeric vector.
    if(dictType == "continuous"){
      phenoChar <- as.character(pheno)
      phenoNum <- suppressWarnings(as.numeric(phenoChar))
      if(length(setdiff(which(is.na(phenoNum)), which(is.na(phenoChar)))) > 0){
        stop(paste("Phenotype", phenoName, "set to continuous but cannot be converted to numeric."))
      }
    }
  }
  
  # Check that usage report file can be created.
  if(!is.null(usage_report_file)){
    tryCatch({
      saveRDS(list(1,2,3), usage_report_file)
      unlink(usage_report_file)
    }, error = function(cond){
      stop(paste(usage_report_file, "could not be created."))
    })
    
  }

  # Compute coexpression of genes
  usageGeneGene <- NA
  results$coexpression <- NA
  if(is.null(usage_report_file)){
    if(compute_gene_cor == TRUE){
      if(verbose == TRUE){
        print("Running gene co-expression")
      }
      results$coexpression = cor(t(expression), use="pairwise.complete.obs")
      if(verbose == TRUE){
        print("Gene co-expression complete")
      }
    }
  }else{
    if(compute_gene_cor == TRUE){
      usageGeneGene <- peakRAM::peakRAM({
        if(verbose == TRUE){
          print("Running gene co-expression")
        }
        results$coexpression = cor(t(expression), use="pairwise.complete.obs")
        if(verbose == TRUE){
          print("Gene co-expression complete")
        }
      })
    }
  }
  
  # Compute coexpression of phenotypes
  usagePhenPhen <- NA
  results$phenocor <- NA
  if(is.null(usage_report_file)){
    if(compute_phenotype_cor == TRUE){
      if(verbose == TRUE){
        print("Running phenotype associations")
      }
      results$phenocor = computePhenotypeCorrelations(phenotype = phenotype,
                                                      phenotype_dictionary = phenotype_dictionary,
                                                      method = assoc_method, pval_adj_method = pval_adj_method)
      if(verbose == TRUE){
        print("Phenotype associations complete")
      }
    }
  }else{
    if(compute_phenotype_cor == TRUE){
      usagePhenPhen <- peakRAM::peakRAM({
        if(verbose == TRUE){
          print("Running phenotype associations")
        }
        results$phenocor = computePhenotypeCorrelations(phenotype = phenotype,
                                                        phenotype_dictionary = phenotype_dictionary,
                                                        method = assoc_method, pval_adj_method = pval_adj_method)
        if(verbose == TRUE){
          print("Phenotype associations complete")
        }
      })
    } 
  }
  
  # Compute association of gene expression with phenotypes and run GSEA
  usagePhenGene <- NA
  results$phenotype_association = list()
  results$GSEA = list()
  if(is.null(usage_report_file)){
    if(compute_gene_phenotype_cor == TRUE){
      if(verbose == TRUE){
        print("Running gene-phenotype associations")
      }
      if(assoc_method %in% c("pearson", "spearman", "kendall")){
        corr <- computeCorrelations(expression = expression, phenotype = phenotype,
                                    pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                    method = assoc_method, pval_adj_method = pval_adj_method)
        results$phenotype_association <- corr$pheno
        results$GSEA <- corr$gsea
      }else{
        linReg <- computeLinearRegression(expression = expression, phenotype = phenotype,
                                          pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                          pval_adj_method = pval_adj_method)
        results$phenotype_association <- linReg$pheno
        results$GSEA <- linReg$gsea
      }
      if(verbose == TRUE){
        print("Gene-phenotype associations complete")
      }
    }
  }else{
    if(compute_gene_phenotype_cor == TRUE){
      usagePhenGene <- peakRAM::peakRAM({
        if(verbose == TRUE){
          print("Running gene-phenotype associations")
        }
        if(assoc_method %in% c("pearson", "spearman", "kendall")){
          corr <- computeCorrelations(expression = expression, phenotype = phenotype,
                                      pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                      method = assoc_method, pval_adj_method = pval_adj_method)
          results$phenotype_association <- corr$pheno
          results$GSEA <- corr$gsea
        }else{
          linReg <- computeLinearRegression(expression = expression, phenotype = phenotype,
                                            pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                            pval_adj_method = pval_adj_method)
          results$phenotype_association <- linReg$pheno
          results$GSEA <- linReg$gsea
        }
        if(verbose == TRUE){
          print("Gene-phenotype associations complete")
        }
      })
    }
  }
  
  # Save the results to a usage file.
  if(!is.null(usage_report_file)){
    saveRDS(list(GeneToGeneCor = usageGeneGene, PhenToPhenCor = usagePhenPhen,
                 PhenToGeneCor = usagePhenGene), usage_report_file)
  }
  
  # Return the results.
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
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' This parameter is only used for the "linear" association method. Use if your data are raw RNA-seq
#' counts. Default is TRUE.
computeLinearRegression <- function(expression, phenotype, phenotype_dictionary, pathways,
                                    pval_adj_method){
  
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

  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  # Run correlations.
  phenotype_vector = as.numeric(pheno)
  cors <- cor(t(expression), phenotype_vector, use = "pairwise.complete.obs", method = method)
  cors <- as.numeric(cors)
  names(cors) <- rownames(expression)
  output_seahorse$cor = cors
  
  # Run GSEA
  cor_rank = sort(output_seahorse$cor, decreasing = TRUE)
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
#'                   
#' @export
gsea_nominal <- function(expression, pheno, pathways){

  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  # Run the linear models.
  phenotype_vector = factor(as.character(pheno))
  # Remove NA values.
  hasVal <- which(!is.na(phenotype_vector))
  phenotype_vector <- phenotype_vector[hasVal]
  expression <- expression[,hasVal]
  design <- model.matrix(~ phenotype_vector)
  fit <- limma::lmFit(object = expression, design = design)
  
  # Compute empirical Bayes statistics.
  # If there are no residual degrees of freedom, catch it and do not return a value.
  tryCatch({
    bayes <- limma::eBayes(fit = fit)
    
    coef_to_test <- setdiff(colnames(design), "(Intercept)")
    anova_res <- limma::topTable(bayes, coef = coef_to_test, number = Inf, sort.by = "none")
    
    # Format p-values as a list for all covariates.
    p <- anova_res[,"P.Value"]
    output_seahorse$cor <- p
    names(output_seahorse$cor) <- rownames(anova_res)
    
    # Run GSEA
    statSorted <- sort(-1 * log10(output_seahorse$cor), decreasing = TRUE)
    fgseaRes <- fgsea::fgsea(pathways, statSorted, minSize=15, maxSize=500, scoreType = "pos")
    output_seahorse$GSEA = fgseaRes
  }, error = function(cond){
    warning("Could not compute empirical Bayes statistics for this phenotypic variable. Returning empty list.")
  })
  
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
  
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  # Run the linear models.
  phenotype_vector = factor(as.character(pheno))

  # Remove NA values.
  hasVal <- which(!is.na(phenotype_vector))
  phenotype_vector <- phenotype_vector[hasVal]
  expression <- expression[,hasVal]
  
  # Check that we still have multiple values. If not, return NA for this covariate.
  output_seahorse$cor <- NA
  output_seahorse$GSEA <- NA
  if(length(unique(phenotype_vector)) > 1){
    design <- model.matrix(~ phenotype_vector)
    fit <- limma::lmFit(object = expression, design = design)
    
    # Compute empirical Bayes statistics.
    tryCatch({
      bayes <- limma::eBayes(fit = fit)

      coef_to_test <- setdiff(colnames(design), "(Intercept)")
      t_res <- limma::topTable(bayes, coef = coef_to_test, number = Inf, sort.by = "none")
      
      # Format p-values as a list for all covariates.
      p <- t_res[,"P.Value"]
      t <- t_res[,"t"]
      output_seahorse$cor <- p
      names(output_seahorse$cor) <- rownames(t_res)
      names(t) <- rownames(t_res)
      
      # Run GSEA
      fgseaRes <- NA
      tryCatch({
        statSorted <- sort(t, decreasing = TRUE)
        fgseaRes <- fgsea::fgsea(pathways, statSorted, minSize=15, maxSize=500)
      }, error = function(cond){
        print(cond)
      })
      
      output_seahorse$GSEA = fgseaRes
    }, error = function(cond){
      warning("Could not compute empirical Bayes statistics for this phenotypic variable. Returning empty list.")
    })
  }else{
    warning("Phenotype has only one level. Returning NA for all gene associations.")
  }
  
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
    
    # Initialize an empty data frame.
    emptyDF <- data.frame(cor = c(),
                          V = c(),
                          testType = c(),
                          row.names = c())
    emptyPadj <- c()
    
    # Run associations.
    if(phenoType == "dichotomous"){
      
      # Do categorical phenotypes.
      whichCat <- which(typesAfter %in% categorical)
      phenoCat <- as.data.frame(phenosAfter[,whichCat])
      colnames(phenoCat) <- colnames(phenosAfter)[whichCat]
      rownames(phenoCat) <- rownames(phenosAfter)
      output_seahorse_cat <- emptyDF
      output_seahorse_padj_cat <- emptyPadj
      if(ncol(phenoCat) > 0){
        output_seahorse_cat <- phenotype_chisq(phenotype = pheno, 
                                               phenotypesToCompare = phenoCat, 
                                               phenotypeType = phenoType)
        output_seahorse_padj_cat <- rep(NA, nrow(output_seahorse_cat))
        output_seahorse_padj_cat[which(output_seahorse_cat$testType == "FFH")] <- 
          stats::p.adjust(output_seahorse_cat[which(output_seahorse_cat$testType == "FFH"), "cor"], method = pval_adj_method)
        output_seahorse_padj_cat[which(output_seahorse_cat$testType == "Chi-square")] <- 
          stats::p.adjust(output_seahorse_cat[which(output_seahorse_cat$testType == "Chi-square"), "cor"], method = pval_adj_method)
      }
      
      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- emptyDF
      output_seahorse_padj_con <- emptyPadj
      if(ncol(phenoCon) > 0){
        output_seahorse_con <- phenotype_ttest(phenotype = pheno, 
                                              phenotypesToCompare = phenoCon, 
                                              phenotypeType = phenoType)
        output_seahorse_padj_con <- stats::p.adjust(output_seahorse_con$cor, method = pval_adj_method)
      }
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_cat$cor, output_seahorse_con$cor),
                               cramerV = c(output_seahorse_cat$V, output_seahorse_con$V),
                               padj = c(output_seahorse_padj_cat, output_seahorse_padj_con),
                               testType = c(output_seahorse_cat$testType, output_seahorse_con$testType),
                               row.names = c(rownames(output_seahorse_cat), rownames(output_seahorse_con)))
      
    }else if(phenoType == "nominal"){
      
      # Do categorical phenotypes.
      whichCat <- which(typesAfter %in% categorical)
      phenoCat <- as.data.frame(phenosAfter[,whichCat])
      colnames(phenoCat) <- colnames(phenosAfter)[whichCat]
      rownames(phenoCat) <- rownames(phenosAfter)
      output_seahorse_cat <- emptyDF
      output_seahorse_padj_cat <- emptyPadj
      if(ncol(phenoCat) > 0){
        output_seahorse_cat <- phenotype_chisq(phenotype = pheno, 
                                                    phenotypesToCompare = phenoCat, 
                                                    phenotypeType = phenoType)
        output_seahorse_padj_cat <- rep(NA, nrow(output_seahorse_cat))
        output_seahorse_padj_cat[which(output_seahorse_cat$testType == "FFH")] <- 
          stats::p.adjust(output_seahorse_cat[which(output_seahorse_cat$testType == "FFH"), "cor"], method = pval_adj_method)
        output_seahorse_padj_cat[which(output_seahorse_cat$testType == "Chi-square")] <- 
          stats::p.adjust(output_seahorse_cat[which(output_seahorse_cat$testType == "Chi-square"), "cor"], method = pval_adj_method)
      }
      
      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- emptyDF
      output_seahorse_padj_con <- emptyPadj
      if(ncol(phenoCon) > 0){
        output_seahorse_con <- phenotype_anova(phenotype = pheno, 
                                              phenotypesToCompare = phenoCon, 
                                              phenotypeType = phenoType)
        output_seahorse_padj_con <- stats::p.adjust(output_seahorse_con$cor, method = pval_adj_method)
      }
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_cat$cor, output_seahorse_con$cor),
                                   cramerV = c(output_seahorse_cat$V, output_seahorse_con$V),
                                   padj = c(output_seahorse_padj_cat, output_seahorse_padj_con),
                                   testType = c(output_seahorse_cat$testType, output_seahorse_con$testType),
                                   row.names = c(rownames(output_seahorse_cat), rownames(output_seahorse_con)))
      
    }else{

      # Do dichotomous phenotypes.
      whichDichotomous <- which(typesAfter == "dichotomous")
      phenoDich <- as.data.frame(phenosAfter[,whichDichotomous])
      colnames(phenoDich) <- colnames(phenosAfter)[whichDichotomous]
      rownames(phenoDich) <- rownames(phenosAfter)
      output_seahorse_dich <- emptyDF
      output_seahorse_padj_dich <- emptyPadj
      if(ncol(phenoDich) > 0){
        output_seahorse_dich <- phenotype_ttest(phenotype = pheno, 
                                                    phenotypesToCompare = phenoDich, 
                                                    phenotypeType = phenoType)
        output_seahorse_padj_dich <- stats::p.adjust(output_seahorse_dich$cor, method = pval_adj_method)
      }
      
      # Do nominal phenotypes.
      whichNominal <- which(typesAfter == "nominal")
      phenoNom <- as.data.frame(phenosAfter[,whichNominal])
      colnames(phenoNom) <- colnames(phenosAfter)[whichNominal]
      rownames(phenoNom) <- rownames(phenosAfter)
      output_seahorse_nom <- emptyDF
      output_seahorse_padj_nom <- emptyPadj
      if(ncol(phenoNom) > 0){
        output_seahorse_nom <- phenotype_anova(phenotype = pheno, 
                                              phenotypesToCompare = phenoNom, 
                                              phenotypeType = phenoType)
        output_seahorse_padj_nom <- stats::p.adjust(output_seahorse_nom$cor, method = pval_adj_method)
      }

      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- emptyDF
      output_seahorse_padj_con <- emptyPadj
      if(ncol(phenoCon) > 0){
        output_seahorse_con = phenotype_cor(phenotype = pheno, 
                                              phenotypesToCompare = phenoCon, 
                                            method = method)
        output_seahorse_padj_con <- rep(NA, nrow(output_seahorse_con))
      }
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_dich$cor, output_seahorse_nom$cor, output_seahorse_con$cor),
                                   cramerV = c(output_seahorse_dich$V, output_seahorse_nom$V, output_seahorse_con$V),
                                   padj = c(output_seahorse_padj_dich, output_seahorse_padj_nom, output_seahorse_padj_con),
                                   testType = c(output_seahorse_dich$testType, output_seahorse_nom$testType, output_seahorse_con$testType),
                                   row.names = c(rownames(output_seahorse_dich), rownames(output_seahorse_nom), rownames(output_seahorse_con)))
      
    }
    # Add the NA values.
    varsNotAdded <- setdiff(colnames(phenotype), row.names(output_seahorse))
    additionalVars <- data.frame(stat = rep(NA, length(varsNotAdded)),
                                 cramerV = rep(NA, length(varsNotAdded)),
                                 padj = rep(NA, length(varsNotAdded)),
                                 testType = rep(NA, length(varsNotAdded)),
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
  typeMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$testType)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))

  # Return all matrices.
  phenoAssocList <- list(stat = statMat, cramerV = vMat, padj = adjMat, testType = typeMat)
  return(phenoAssocList)
}

#' Function to run associations between categorical (dichotomous or nominal) phenotypes.
#' If at least 5 samples exist in each group, run a Chi-square test and compute Cramer's V.
#' Otherwise, run a Fisher-Freeman-Halton Test.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is dichotomous or nominal.
#' phenotypes against which to compare.
#' @returns A data frame with two vectors: "cor", which lists the Chi-square or Fisher-Freeman-Halton Test p-values, 
#' and "V" which lists the Cramer's V statistics (will be NA if Fisher-Freeman-Halton Test is computed)
phenotype_chisq <- function(phenotype, phenotypesToCompare, phenotypeType){
  
  # Generate the tables for all phenotype pairs.
  allPhenoPairTables <- lapply(phenotypesToCompare, function(phen){
    grpCounts <- table(phenotype, phen)
  })
  names(allPhenoPairTables) <- colnames(phenotypesToCompare)
  
  # Do a chi-square test everywhere else.
  # Run the comparisons one at a time because Chi-square and FFH cannot be scaled.
  chisqRes <- lapply(allPhenoPairTables, function(chisqTable) {
    
    # Check that we have 2 or more non-zero marginals. Only run the test if we do.
    rowMarginals <- rowSums(chisqTable)
    colMarginals <- colSums(chisqTable) 
    nonzeroRowMarginalCount <- length(which(rowMarginals > 0))
    nonzeroColMarginalCount <- length(which(colMarginals > 0))
    pval <- NA
    V <- NA
    type <- "Chi-square"
    if(nonzeroRowMarginalCount >= 2 && nonzeroColMarginalCount >= 2){
      # Run the test. If a warning is thrown, switch to FFH.
      chisq <- withCallingHandlers(
        chisq.test(chisqTable),
        warning = function(w) {
          if (grepl("Chi-squared approximation may be incorrect", w$message)) {
            # Do not change only within the warning scope but also within the parent scope.
            type <<- "FFH"
            invokeRestart("muffleWarning")
          }
        }
      )
      
      # Get the p-value and Cramer's V.
      pval <- chisq$p.value
      V <- sqrt(chisq$statistic / (sum(chisqTable) * min(nrow(chisqTable) - 1, 
                                                         ncol(chisqTable) - 1)))
    }else{
      warning("Less than 2 nonzero marginals in the contingency table - Chi-square will return NA")
    }
    return(list(pval = pval, v = V, testType = type))
  })
  chisqP <- unlist(lapply(chisqRes, function(res){
    return(res$pval)
  }))
  cramerV <- unlist(lapply(chisqRes, function(res){
    return(res$v)
  }))
  testType <- unlist(lapply(chisqRes, function(res){
    return(res$testType)
  }))
  corResNotLessFull <- data.frame(cor = chisqP,
                              V = cramerV,
                              testType = testType,
                              row.names = colnames(phenotypesToCompare))
  
  # Run FFS where warning was issued.
  corResNotLess <- corResNotLessFull[which(corResNotLessFull$testType == "Chi-square"),]
  switchedTables <- allPhenoPairTables[which(corResNotLessFull$testType == "FFH")]
  corResSwitched <- phenotype_ffs(switchedTables)

  # Bind together the results from the FFH and Chi-Square tests.
  corRes <- rbind(corResNotLess, corResSwitched)
  corRes <- corRes[colnames(phenotypesToCompare),]

  return(corRes)
}

#' Function to run associations between categorical (dichotomous or nominal) phenotypes using a Fisher-Freeman-Halton Test.
#' @param tables The list of tables to test, named by the phenotypes against which to compare.
#' @returns A data frame with two vectors: "cor", which lists Fisher-Freeman-Halton Test p-values, 
#' and "V" which lists the Cramer's V statistics (NA)
phenotype_ffs <- function(tables){
  
  # Run the comparisons one at a time because FFH cannot be scaled.
  pvals <- unlist(lapply(tables, function(fishTable) {
    fishP <- fisher.test(fishTable, simulate.p.value = TRUE, B = 10000)$p.value
    return(fishP)
  }))

  # Return the results.
  return(data.frame(cor = pvals,
                    V = rep(NA, length(tables)),
                    testType = rep("FFH", length(tables)),
                    row.names = names(tables)))
}

#' Function to compute ANOVA statistics between nominal and continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is nominal or continuous.
#' @returns A data frame with two vectors: "cor", which lists the ANOVA p-values, 
#' and "V" which is set to NA.
phenotype_anova <- function(phenotype, phenotypesToCompare, phenotypeType){
  
  cor <- rep(NA, ncol(phenotypesToCompare))
  
  # Case 1 - the phenotype is nominal and the phenotypes to compare are numeric.
  # Case 2 - the phenotype is numeric and the phenotypes to compare are nominal.
  if(phenotypeType == "nominal"){
    phenotype_vector = factor(as.character(phenotype))
    cor <- unlist(apply(phenotypesToCompare, MARGIN=2, function(x){
      results <- NA
      tryCatch({
        results <- anova(lm(as.numeric(as.character(x))~phenotype_vector))$`Pr(>F)`[1]
      }, error = function(cond){
        warning("In this phenotype pair, all continuous values are missing for all but one phenotype level - NA result will be returned")
      })
      return(results)
    }))
  }else{
    cor <- unlist(lapply(1:ncol(phenotypesToCompare), function(i){
      results <- NA
      phenotype_vector <- factor(as.character(phenotypesToCompare[,i]))
      tryCatch({
        results <- anova(lm(as.numeric(as.character(phenotype))~phenotype_vector))$`Pr(>F)`[1]
      }, error = function(cond){
        warning("In this phenotype pair, all continuous values are missing for all but one phenotype level - NA result will be returned")
      })
      return(results)
    }))
  }
    
  return(data.frame(cor = cor,
                    V = rep(NA, ncol(phenotypesToCompare)),
                    testType = rep("ANOVA", ncol(phenotypesToCompare)),
                    row.names = colnames(phenotypesToCompare)))
}

#' Function to compute Welch's t-test statistics between dichotomous and continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is dichotomous or continuous.
#' @returns A data frame with two vectors: "cor", which lists the t-test p-values, 
#' and "V" which is set to NA.
phenotype_ttest <- function(phenotype, phenotypesToCompare, phenotypeType){
  
  # Case 1 - the phenotype is dichotomous and the phenotypes to compare are numeric.
  # We can vectorize this and use matrixTests.
  # Case 2 - the phenotype is numeric and the phenotypes to compare are dichotomous.
  # We cannot vectorize this and use t.test instead, which defaults to Welch's t-test.
  if(phenotypeType == "dichotomous"){
    
    # Split groups.
    phenotype_vector = factor(as.character(phenotype))
    levels <- unique(phenotype_vector)
    group1 <- phenotypesToCompare[which(phenotype_vector == levels[1]),]
    group2 <- phenotypesToCompare[which(phenotype_vector == levels[2]),]
    
    # Check that thresholds are met.
    whichLevel1 <- length(which(phenotype_vector == levels[1]))
    whichLevel2 <- length(which(phenotype_vector == levels[2]))
    cor <- NA
    
    if(length(dim(group1)) >= 2 || length(dim(group2)) >= 2){
      meetsThreshold <- colSums(!is.na(group1)) >= 2 & colSums(!is.na(group2)) >= 2
      
      # Initialize result to NA.
      cor <- rep(NA, ncol(group1))
      names(cor) <- colnames(group1)
      
      # Only compute correlations if both levels are represented in the phenotype vector.
      if(whichLevel1 > 0 && whichLevel2 > 0 && length(which(meetsThreshold == TRUE)) > 0){
        # Compute t-test where thresholds are met.
        tresValid <- matrixTests::col_t_welch(
          group1[, meetsThreshold, drop = FALSE],
          group2[, meetsThreshold, drop = FALSE]
        )
        
        # Compute remaining t-tests.
        cor[meetsThreshold] <- tresValid$pvalue
      }
    }else{
      meetsThreshold <- sum(!is.na(group1)) >= 2 & sum(!is.na(group2)) >= 2
      
      # Initialize result to NA.
      names(cor) <- names(group1)
      
      # Only compute correlations if both levels are represented in the phenotype vector.
      if(whichLevel1 > 0 && whichLevel2 > 0 && length(which(meetsThreshold == TRUE)) > 0){
        tryCatch({
          # Compute t-test where thresholds are met.
          tresValid <- t.test(group1[meetsThreshold, drop = FALSE], 
                              group2[meetsThreshold, drop = FALSE])
          
          # Compute remaining t-tests.
          cor[meetsThreshold] <- tresValid$p.value
        }, error = function(cond){
          warning("Could not compute t-test (it is possible that variance is too low). Returning NA.")
        })
      }
    }
    
  }else{
    cor <- unlist(lapply(1:ncol(phenotypesToCompare), function(i){
      phenotype_vector = factor(as.character(phenotypesToCompare[,i]))
      levels <- unique(phenotype_vector)
      group1 <- phenotype[which(phenotype_vector == levels[1])]
      group2 <- phenotype[which(phenotype_vector == levels[2])]
      stat <- NA
      if(length(which(!is.na(group1))) > 2 && length(which(!is.na(group2))) > 2){
        tryCatch({
          tres <- t.test(group1, group2)
          stat = tres$p.value
        }, error = function(cond){
          warning("Could not compute t-test (it is possible that variance is too low). Returning NA.")
        })
      }
      names(stat) <- colnames(phenotypesToCompare)[i]
      return(stat)
    }))
  }
  if(length(which(is.na(cor))) > 0){
    warning("Some phenotypes did not have sufficient sample sizes to perform a t-test - NAs will be returned")
  }
  
  # Return the data frame.
  return(data.frame(cor = cor,
                    V = rep(NA, ncol(phenotypesToCompare)),
                    testType = rep("T-Test", ncol(phenotypesToCompare)),
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
  corRes = cor(phenotype, phenotypesToCompare, use="pairwise.complete.obs", method = method)
  return(data.frame(cor = unname(c(corRes)),
                    testType = rep("Cor", ncol(phenotypesToCompare)),
                    V = rep(NA, ncol(phenotypesToCompare)),
                    row.names = colnames(phenotypesToCompare)))
}

#' Generates the following files which are needed for use with the SEAHORSE UI:
#' - all_gsea_results.tsv.gz: A list of pathway enrichment scores with respect to all
#' phenotypes and all files, formatted as pathway	pval	padj	log2err	ES	NES	size	ranks	leadingEdge	varname	tissue
#' - data_dictionary.tsv.gz: A list of phenotype names, descriptions, and variable types formatted as
#' VARNAME VARDESC VARMETA TYPE
#' - geneexpression2geneexpression.tsv.gz: Gene-gene association statistics, formatted as Gene A	Gene B	Tissue	Correlation
#' - geneexpression_data.tsv.gz: Gene expression data, formatted as ENSG    SAMPID  GENE_EXPRESSION
#' - human_ensembl2symbol_map.tsv.gz: A mapping from Ensembl IDs to HGNC symbols and Entrez IDs,
#' formatted as ALIAS	ENSEMBL	SYMBOL	ENTREZID
#' - metadata.tsv.gz: The phenotype data, formatted as SAMPID	tissue	VARNAME	VALUE
#' - metadata2expression.tsv.gz: Gene-phenotype association statistics, formatted as 
#' VARNAME GENE    tissue  TEST    TESTSTAT        TESTPVALUE
#' - metadata2metadata.tsv.gz: Phenotype-phenotype association statistics, formatted as 
#' VARNAME1        VARNAME2        tissue  TEST    TESTSTAT        TESTPVALUE
#' @param input_directory : A directory containing all inputs to SEAHORSE, stored as RDS files,
#' where the input values contain slots for expression, phenotype, and phenotype dictionary. 
#' @param result_directory : A directory containing all SEAHORSE results from the expression data,
#' phenotype data, phenotype dictionary, and pathways provided, stored as RDS files. Files must follow
#' the same naming convention as in the result directory.
#' @param output_directory: The directory where you wish to store the files generated by this function.
#' @param data_dictionary: The data dictionary.
#' @param pathways The pathways used to run SEAHORSE
#' @return NULL
#' @export
seahorseFormatForUI <- function(input_directory, result_directory, output_directory,
                                data_dictionary, pathways){
  
  # Check packages.
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required but not installed.")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package 'org.Hs.eg.db' is required but not installed.")
  }
  
  # Create output directory if it doesn't already exist.
  if(!dir.exists(output_directory)){
    dir.create(output_directory)
  }
  
  # Check format of data dictionary.
  if(is(data_dictionary, "data.frame")){
    if(!identical(colnames(data_dictionary), c("VARNAME", "VARDESC", "VARMETA", "TYPE"))){
      stop("Incorrect format for data dictionary")
    }
  }else{
    stop("Incorrect format for data dictionary")
  }
  
  # Loop through all files.
  resultfiles <- list.files(result_directory)
  inputfiles <- list.files(input_directory)
  tissue_names_result <- formatTissueNames(resultfiles)
  tissue_names_input <- formatTissueNames(inputfiles)
  
  # Check that names match.
  if(length(tissue_names_result) == length(tissue_names_input)){
    if(!identical(tissue_names_result, tissue_names_input)){
      stop("Result and input file names do not match")
    }
  }else{
    stop("Number of result and input files do not match")
  }
  
  # If the output directory does not exist, create it.
  if(!dir.exists(output_directory)){
    dir.create(output_directory)
  }
  
  # Write results.
  message("Preparing to write expression...")
  writeExpression(inFiles = paste(input_directory, inputfiles, sep = "/"),
                  file = paste(output_directory, "geneexpression_data", sep = "/"))
  message("Expression file done.")
  message("Preparing to write gene mapping...")
  writeMapping(inFiles = paste(input_directory, inputfiles, sep = "/"),
               file = paste(output_directory, "human_ensembl2symbol_map", sep = "/"))
  message("Mapping file done.")
  message("Preparing to write data dictionary...")
  writeToGz(data = data_dictionary, file = paste(output_directory, "data_dictionary", sep = "/"))
  message("Data dictionary done.")
  message("Preparing to write GSEA results...")
  writePathwayEnrichment(inFiles = paste(result_directory, resultfiles, sep = "/"),
                         file = paste(output_directory, "all_gsea_results", sep = "/"),
                         pathways = pathways)
  message("GSEA result file done.")
  message("Preparing to write phenotype data...")
  writePhenotype(inFiles = paste(input_directory, inputfiles, sep = "/"),
                         file = paste(output_directory, "metadata", sep = "/"))
  message("Phenotype data done.")
  message("Preparing to write phenotype-gene associations...")
  writePhenotypeGeneAssociations(inFiles = paste(input_directory, inputfiles, sep = "/"),
                  resultFiles = paste(result_directory, resultfiles, sep = "/"),
                 file = paste(output_directory, "metadata2expression", sep = "/"))
  message("Phenotype-gene associations done.")
  message("Preparing to write phenotype associations...")
  writePhenotypePhenotypeAssociations(inFiles = paste(result_directory, resultfiles, sep = "/"),
                                 file = paste(output_directory, "metadata2metadata", sep = "/"))
  message("Phenotype associations done")
  message("Preparing to write gene associations...")
  writeGeneGene(inFiles = paste(result_directory, resultfiles, sep = "/"),
                                      file = paste(output_directory, "geneexpression2geneexpression", sep = "/"))
  message("Gene associations done.")
  
}

#' Generates the pathway enrichment scores across all tissues and writes them.
#' Format is pathway	pval	padj	log2err	ES	NES	size	ranks	leadingEdge	varname	tissue
#' @param inFiles The SEAHORSE files
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @param pathways The pathways used to run SEAHORSE
#' @return NULL
writePathwayEnrichment <- function(inFiles, file, pathways){
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))
  
  # Write the pathway enrichment results.
  pathwayEnrichment <- tryCatch({for(i in 1:length(inFiles)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])
    
    enrichmentVar <- data.frame(pathway = character(), pval = numeric(), 
                                padj = numeric(), log2err = numeric(),
                                ES = numeric(), NES = numeric(),
                                size = numeric(), ranks = character(),
                                leadingEdge = character(), varname = character(),
                                tissue = character())
    
    if(length(results$GSEA) > 0){
      # Loop through variable names.
      enrichmentVar <- do.call(rbind, lapply(names(results$GSEA), function(var){
        
        # Get full ranking of genes for this var based on whether correlation or p-value was used.
        # P-values are log-transformed and absolute values are used for correlations.
        # The largest values should then receive the lowest (best) ranking.
        transformedStats <- results$phenotype_association[[var]]$stat
        if(!all(is.null(transformedStats))){
          transformedStats[which(!is.na(results$phenotype_association[[var]]$padj))] <- -log10(transformedStats[which(!is.na(results$phenotype_association[[var]]$padj))])
          transformedStats[which(is.na(results$phenotype_association[[var]]$padj))] <- abs(transformedStats[which(is.na(results$phenotype_association[[var]]$padj))])
          pvalRanking <- order(-transformedStats)
          names(pvalRanking) <- rownames(results$phenotype_association[[var]])
        }
        
        # Build basic pathway data.
        returnDat <- data.frame(pathway = character(), pval = numeric(), 
                                padj = numeric(), log2err = numeric(),
                                ES = numeric(), NES = numeric(),
                                size = numeric(), ranks = character(),
                                leadingEdge = character(), varname = character(),
                                tissue = character())
        if(!all(is.na(results$GSEA[[var]]))){
          returnDat <- results$GSEA[[var]][,c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")]
          # Get ranking of pathway genes.
          varResults <- as.data.frame(results$GSEA[[var]])
          rankVec <- unlist(lapply(1:nrow(varResults), function(pw){
            # Get pathway genes.
            pathwayGenes <- pathways[[varResults[pw, "pathway"]]]
            # Get adjusted p-value rankings.
            rankForLeadingEdge <- pvalRanking[pathwayGenes]
            return(paste0("{", paste(rankForLeadingEdge, collapse = ","), "}"))
          }))
          returnDat$ranks <- rankVec
          # Format leading edges as string.
          leadingEdgeVec <- unlist(lapply(results$GSEA[[var]][,"leadingEdge"][[1]], function(edge){
            return(paste0("{", paste(edge, collapse = ","), "}"))
          }))
          returnDat$leadingEdge <- leadingEdgeVec
          
          # Add variable name.
          returnDat$varname <- rep(var, nrow(returnDat))
          
          # Add tissue.
          splitPath <- strsplit(inFiles[i], "/")[[1]]
          localFile <- splitPath[length(splitPath)]
          withoutFileExt <- strsplit(localFile, ".RDS")[[1]][1]
          returnDat$tissue <- rep(withoutFileExt, nrow(returnDat))
        }
        return(returnDat)
      }))
    }
    
    # Write the chunk.
    write.table(enrichmentVar, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
    
  }}, error = function(cond){
    print(cond)
    stop(paste("Could not format pathway enrichment"))})
}

#' Formats the tissue names by removing the .RDS extension from the input files.
#' @param names The list of tissue names.
#' @return A list of names.
formatTissueNames <- function(names){
  return(unlist(lapply(names, function(name){
    nameSplit <- strsplit(name, "/")[[1]]
    localName <- nameSplit[length(nameSplit)]
    result <- strsplit(localName, ".RDS")[[1]][1]
    return(result)
  })))
}

#' Generates the gene correlation results in the format expected by the UI.
#' Format is Gene A	Gene B	Tissue	Correlation
#' @param inFiles The SEAHORSE files
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @return NULL 
writeGeneGene <- function(inFiles, file){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)
  
  # Write the gene correlation results.
  geneCor <- for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])
    
    # Remove duplicate pairs.
    geneCor <- results$coexpression
    flatTable <- data.frame(A = character(), B = character(), 
                            Tissue = character(), Correlation = numeric())
    colnames(flatTable) <- c("Gene A", "Gene B", "Tissue", "Correlation")
    if(length(which(!is.na(geneCor) == TRUE)) > 0){
      geneCor[lower.tri(geneCor, diag = TRUE)] <- NA
      
      flatTable <- data.frame(
        A = rep(colnames(geneCor), times = nrow(geneCor)),
        B = rep(rownames(geneCor), each = ncol(geneCor)),
        Correlation = as.vector(geneCor)
      )
      
      # Remove NA values.
      flatTable <- na.omit(flatTable)
      
      # Add names of variables.
      colnames(flatTable) <- c("Gene A", "Gene B", "Correlation")
      
      # Add tissue.
      flatTable$Tissue <- rep(tissueNames[i], nrow(flatTable))
      
      # Rearrange columns.
      flatTable <- flatTable[,c("Gene A", "Gene B", "Tissue", "Correlation")]
    }
    
    # Write the chunk.
    write.table(flatTable, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
}

#' Generates the expression data in the format expected by the UI. Writes a
#' data frame formatted as ENSG    SAMPID  GENE_EXPRESSION as .tsv.gz.
#' @param inFiles The files to input to SEAHORSE, formatted as RDS.
#' @param file The name of the file to write. .tsv.gz. will be appended.
#' @return NULL
writeExpression <- function(inFiles, file){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)

  # Write the pathway enrichment results.
  expression <- tryCatch({for(i in 1:length(inFiles)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])

    # Flatten it.
    exprMat <- as.matrix(results$expression)
    flatTable <- data.frame(
      ENSG = rep(rownames(exprMat), times = ncol(exprMat)),
      SAMPID = rep(colnames(exprMat), each = nrow(exprMat)),
      GENE_EXPRESSION = as.vector(exprMat),
      tissue = tissueNames[i]
    )
    
    write.table(flatTable, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
  }, error = function(cond){
    print(cond)
    stop("Could not format expression")})
}

#' Writes a data frame as .csv.gz.
#' @param data The data frame
#' @param file The name of the file to write. .tsv.gz. will be appended.
#' @return NULL
writeToGz <- function(data, file){
  gz <- gzfile(paste0(file, ".tsv.gz"), "w")
  write.table(data, file = gz, sep = "\t", row.names = FALSE, quote = FALSE)
  close(gz)
}

#' Generates the gene mapping data in the format expected by the UI. Writes as
#' .tsv.gz formatted as ALIAS	ENSEMBL	SYMBOL	ENTREZID.
#' @param inFiles The files to input to SEAHORSE, formatted as RDS.
#' @param file The name of the file to write. .tsv.gz. will be appended.
#' @return NULL
writeMapping <- function(inFiles, file){
  
  # Get all unique genes.
  ensembl <- tryCatch({unique(unlist(lapply(1:length(inFiles), function(i){
    
    # Return genes from each file.
    results = readRDS(inFiles[i])
    return(rownames(results$expression))
  })))}, error = function(cond){
    stop("Incorrect format for input file")})

  # Map genes.
  mappedSymbols <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db,
    keys = ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first" # Returns one value for each ENSEMBL ID.
  )
  mappedEntrez <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db,
    keys = ensembl,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first" # Returns one value for each ENSEMBL ID.
  )
  mappedGenes <- data.frame(ENSEMBL = ensembl, SYMBOL = mappedSymbols, ENTREZID = mappedEntrez)
  mappedGenes$ALIAS <- mappedGenes$SYMBOL
  mappedGenes <- mappedGenes[,c("ALIAS", "ENSEMBL", "SYMBOL", "ENTREZID")]

  # Write the file.
  writeToGz(data = mappedGenes, file = file)
}

#' Generates the phenotype data in the format expected by the UI.
#' Format is SAMPID	tissue	VARNAME	VALUE
#' @param inFiles The files to input to SEAHORSE, formatted as RDS.
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @return NULL
writePhenotype <- function(inFiles, file){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)

  # Write the pathway enrichment results.
  pheno <- tryCatch({for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])

    # Convert to a table.
    phenoMat <- as.matrix(results$phenotype)
    flatTable <- data.frame(
      VARNAME = rep(colnames(phenoMat), each = nrow(phenoMat)),
      SAMPID = rep(rownames(phenoMat), times = ncol(phenoMat)),
      VALUE = as.vector(phenoMat)
    )

    # Add tissue
    flatTable$tissue <- rep(tissueNames[i], nrow(flatTable))

    # Rearrange columns.
    flatTable <- flatTable[,c("SAMPID", "tissue", "VARNAME", "VALUE")]
    
    write.table(flatTable, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))

  }}, error = function(cond){
    stop("Incorrect format for input file")})
}

#' Generates the phenotype data in the format expected by the UI.
#' Format is VARNAME GENE    tissue  TEST    TESTSTAT        TESTPVALUE
#' @param inFiles The SEAHORSE inputs
#' @param resultFiles The SEAHORSE results
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @return NULL
writePhenotypeGeneAssociations <- function(inFiles, resultFiles, file){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)
  
  # Write the pathway enrichment results.
  phenotypeGene <- for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results <- readRDS(resultFiles[i])
    inputs <- readRDS(inFiles[i])
    
    phenotypeToGeneDf <- data.frame(VARNAME = character(), GENE = character(), 
                                    tissue = character(), TEST = character(),
                                    TESTSTAT = numeric(), TESTPVALUE = numeric())

    if(length(results$phenotype_association) > 0){
        phenotypeToGeneList <- lapply(names(results$phenotype_association), function(var){
          outDf <- data.frame(VARNAME = character(), GENE = character(), 
                              tissue = character(), TEST = character(),
                              TESTSTAT = numeric(), TESTPVALUE = numeric())
          if(var %in% colnames(inputs$phenotype)){
            outDf <- results$phenotype_association[[var]]
            if(length(outDf) > 1){
              colnames(outDf) <- c("TESTSTAT", "TESTPVALUE")
              outDf$GENE <- rownames(outDf)
              outDf$VARNAME <- rep(var, nrow(outDf))
              outDf$tissue <- tissueNames[i]
              outDf$TEST <- "Correlation"
              dictVal <- inputs$dict[which(colnames(inputs$phenotype) == var)]
              if(dictVal == "dichotomous"){
                outDf$TEST <- "LIMMA Moderated t-test"
              }else if(dictVal == "nominal"){
                outDf$TEST <- "ANOVA"
              }
              outDf <- outDf[,c("VARNAME", "GENE", "tissue", "TEST", "TESTSTAT", "TESTPVALUE")]
            }
          }
          
          return(outDf)
        })
        phenotypeToGeneDf <- do.call(rbind, phenotypeToGeneList)
    }
    write.table(phenotypeToGeneDf, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
}

#' Generates the phenotype results in the format expected by the UI.
#' Format is VARNAME1        VARNAME2        tissue  TEST    TESTSTAT        TESTPVALUE
#' @param inFiles The SEAHORSE results
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @return NULL
writePhenotypePhenotypeAssociations <- function(inFiles, file){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)
  
  # Write the pathway enrichment results.
  phenotype <- for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])
    
    phenotypeToPhenotypeDf <- data.frame(VARNAME1 = character(), VARNAME2 = character(), tissue = character(),
                                         TEST = character(), TESTSTAT = numeric(), TESTPVALUE = numeric())
    
    if(all(!is.na(results$phenocor)) == TRUE){
      phenotypeToPhenotypeList <- lapply(1:(nrow(results$phenocor$stat)-1), function(j){
        VARNAME2 <- colnames(results$phenocor$stat)[(j+1):ncol(results$phenocor$stat)]
        TEST <- results$phenocor$testType[j,(j+1):ncol(results$phenocor$stat)]
        TESTSTAT <- results$phenocor$stat[j,(j+1):ncol(results$phenocor$stat)]
        TESTPVALUE <- results$phenocor$padj[j,(j+1):ncol(results$phenocor$padj)]
        VARNAME1 <- rep(rownames(results$phenocor$stat)[j], length((j+1):ncol(results$phenocor$padj)))
        tissue <- tissueNames[i]
        return(data.frame(VARNAME1 = VARNAME1, VARNAME2 = VARNAME2, tissue = tissue,
                          TEST = TEST, TESTSTAT = TESTSTAT, TESTPVALUE = TESTPVALUE))
      })
      phenotypeToPhenotypeDf <- do.call(rbind, phenotypeToPhenotypeList)
    }
    
    write.table(phenotypeToPhenotypeDf, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
}