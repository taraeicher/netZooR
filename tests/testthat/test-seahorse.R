context("test SEAHORSE result")

test_that("seahorse function works", {
  skip_if_not_installed("fgsea")
  set.seed(42)
  # Simulate expression data
  expression_data = data.frame(matrix(rexp(1000, rate=.1), ncol=10, nrow = 100))
  rownames(expression_data) = paste("gene", 1:100, sep = "")
  colnames(expression_data) = paste("sample", 1:10, sep = "")
  
  # Simulate phenotypic data
  phenotype_data = data.frame(matrix(0, ncol=2, nrow = 10))
  colnames(phenotype_data) = c("sex", "height")
  rownames(phenotype_data) = colnames(expression_data)
  phenotype_data$sex = c(rep("male", nrow(phenotype_data)/2), rep("female", nrow(phenotype_data)/2))
  phenotype_data$height = 65 + sample.int(10, nrow(phenotype_data), replace = T)
  
  phenotype_dictionary = c("categorical", "numeric")
  
  # Create toy pathways
  pathways = list()
  pathways$pathway1 = sample(rownames(expression_data), 50)
  pathways$pathway2 = sample(rownames(expression_data), 30)
  pathways$pathway3 = sample(rownames(expression_data), 70)
  
  # Check that seahorse returns an error if the p-value adjustment method
  # is invalid.
  expect_error(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                        pval_adj_method = "myCoolAdjustmentMethod"),
               "myCoolAdjustmentMethod is not a valid method for stats::p.adjust().")
  
  # Run seahorse
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways)

  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height") %in% names(results$GSEA)))
  expect_true(all(!is.na(results$coexpression)))
  
  # Run seahorse without correlation matrix
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_cor = FALSE)
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height") %in% names(results$GSEA)))
  expect_equal(results$phenotype_association$sex$stat, results$phenotype_association$sex$padj)
  expect_equal(results$phenotype_association$height$padj, rep("NA", length(results$phenotype_association$height$padj)))
  expect_true(is.na(results$coexpression))
  
  # Run seahorse with bonferroni correction
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_cor = FALSE, pval_adj_method = "bonferroni")
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sex$stat, method = "bonferroni"), 
               results$phenotype_association$sex$padj)
  expect_equal(results$phenotype_association$height$padj, rep("NA", length(results$phenotype_association$height$padj)))
  expect_true(is.na(results$coexpression))
  
  # Run seahorse with fdr correction
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_cor = FALSE, pval_adj_method = "fdr")
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sex$stat, method = "fdr"), 
               results$phenotype_association$sex$padj)
  expect_equal(results$phenotype_association$height$padj, rep("NA", length(results$phenotype_association$height$padj)))
  expect_true(is.na(results$coexpression))
  
  
  # Run SEAHORSE with linear regression.
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_cor = FALSE, assoc_method = "linear")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height") %in% names(results$GSEA)))
  expect_equal(results$phenotype_association$sexmale$stat, results$phenotype_association$sexmale$padj)
  expect_equal(results$phenotype_association$height$stat, results$phenotype_association$height$padj)
  expect_true(is.na(results$coexpression))
  
  # Run SEAHORSE with linear regression and Bonferroni adjustment
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_cor = FALSE, assoc_method = "linear",
                      pval_adj_method = "bonferroni")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sexmale$stat, method = "bonferroni"),
               results$phenotype_association$sexmale$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$height$stat, method = "bonferroni"),
               results$phenotype_association$height$padj)
  expect_true(is.na(results$coexpression))
  
  # Run SEAHORSE with linear regression and FDR adjustment
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_cor = FALSE, assoc_method = "linear", pval_adj_method = "fdr")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sexmale$stat, method = "fdr"),
               results$phenotype_association$sexmale$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$height$stat, method = "fdr"),
               results$phenotype_association$height$padj)
  expect_true(is.na(results$coexpression))
  
  
  # Run SEAHORSE with linear regression and the correlation matrix.
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_cor = TRUE, assoc_method = "linear")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height") %in% names(results$GSEA)))
  expect_true(all(!is.na(results$coexpression)))
})