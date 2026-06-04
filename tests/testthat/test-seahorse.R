context("test SEAHORSE result")

test_that("seahorse function works", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("matrixTests")
  set.seed(42)
  # Simulate expression data
  expression_data = data.frame(matrix(rexp(1000, rate=.1), ncol=10, nrow = 100))
  rownames(expression_data) = paste("gene", 1:100, sep = "")
  colnames(expression_data) = paste("sample", 1:10, sep = "")
  
  # Simulate phenotypic data
  phenotype_data = data.frame(matrix(0, ncol=3, nrow = 10))
  colnames(phenotype_data) = c("sex", "height", "group")
  rownames(phenotype_data) = colnames(expression_data)
  phenotype_data$sex = c(rep("male", nrow(phenotype_data)/2), rep("female", nrow(phenotype_data)/2))
  phenotype_data$height = 65 + sample.int(10, nrow(phenotype_data), replace = T)
  phenotype_data$group = c(rep("1", 3), rep("2", 4), rep("3", 3))
  
  phenotype_dictionary = c("dichotomous", "continuous", "nominal")
  
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

  # Check that seahorse returns an error if the phenotype dictionary entry is invalid.
  phenotype_bad1 <- c("nominal", "continuous", "nominal")
  expect_error(seahorse(expression_data, phenotype_data, phenotype_bad1, pathways,
                        pval_adj_method = "bonferroni"),
               "Phenotype sex set to nominal but has 2 levels or fewer.")
  phenotype_bad2 <- c("dichotomous", "continuous", "dichotomous")
  expect_error(seahorse(expression_data, phenotype_data, phenotype_bad2, pathways,
                        pval_adj_method = "bonferroni"),
               "Phenotype group set to dichotomous but has more than 2 levels.")
  phenotype_bad3 <- c("continuous", "continuous", "nominal")
  expect_error(seahorse(expression_data, phenotype_data, phenotype_bad3, pathways,
                        pval_adj_method = "bonferroni"),
               "Phenotype sex set to continuous but cannot be converted to numeric.")
  
  # Run seahorse
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways)

  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "phenocor", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height", "group") %in% names(results$GSEA)))
  expect_true(all(!is.na(results$coexpression)))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Run seahorse to test both FFH and Chi-square tests and ANOVA, t-test, and correlation.
  # We expect that "smoke" will trigger FFH in all except sex because true counts < 5,
  # "sex" will trigger FFH in "group" because true counts < 5,
  # "group" will trigger FFH in "rare group" because true counts < 5,
  # and "group" will trigger FFH in "rare group" because EXPECTED (not true) counts < 5.
  # and "grade" will trigger FFH in "rare group" because true counts < 5.
  # All other comparisons should trigger Chi-square tests.
  phenotype_data_2 <- do.call(rbind, rep(list(phenotype_data), 10))
  phenotype_data_2$smoke = c(rep("No", 90), rep("Yes", 10))
  phenotype_data_2$grade = c(rep("A", 45), rep("B", 45), rep("C", 10))
  phenotype_data_2$rare_group <- rep("common", 100)
  phenotype_data_2$rare_group[c(1:12, 14, 18, 19)] <- "rare"
  phenotype_data_2$WBC <- runif(n = 100, min = 4000, max = 11000)
  phenotype_data_2 <- phenotype_data_2[,c("smoke", "sex", "group", "height", "grade", "rare_group", "WBC")]
  phenotype_dictionary_2 <- c("dichotomous", "dichotomous", "nominal", "continuous", "nominal", "dichotomous", "continuous")
  expression_data_2 = data.frame(matrix(rexp(1000, rate=.1), ncol=100, nrow = 100))
  rownames(expression_data_2) = paste("gene", 1:100, sep = "")
  colnames(expression_data_2) = paste("sample", 1:100, sep = "")
  results <- seahorse(expression_data_2, phenotype_data_2, phenotype_dictionary_2, pathways)
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "phenocor", "GSEA") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("smoke", "sex", "height", "group", "grade", "rare_group", "WBC") %in% names(results$GSEA)))
  expect_true(all(!is.na(results$coexpression)))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("smoke", "sex", "height", "group", "grade", "rare_group", "WBC") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("smoke", "sex", "height", "group", "grade", "rare_group", "WBC") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("smoke", "sex", "height", "group", "grade", "rare_group", "WBC") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("smoke", "sex", "height", "group", "grade", "rare_group", "WBC") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("smoke", "sex", "height", "group", "grade", "rare_group", "WBC") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("smoke", "sex", "height", "group", "grade", "rare_group", "WBC") %in% colnames(results$phenocor$padj)))
  # Ensure that statistics are calculated as expected (chi-square and FFH)
  expect_true(all(!is.na(results$phenocor$stat["smoke", c("sex", "group", "grade", "rare_group")])))
  expect_true(all(!is.na(results$phenocor$stat["sex", c("group", "grade", "rare_group")])))
  expect_true(all(!is.na(results$phenocor$stat["group", c("grade", "rare_group")])))
  expect_true(all(!is.na(results$phenocor$stat["grade", "rare_group"])))
  expect_true(all(!is.na(results$phenocor$padj["smoke", c("sex", "group", "grade", "rare_group")])))
  expect_true(all(!is.na(results$phenocor$padj["sex", c("group", "grade", "rare_group")])))
  expect_true(all(!is.na(results$phenocor$padj["group", c("grade", "rare_group")])))
  expect_true(all(!is.na(results$phenocor$padj["grade", "rare_group"])))
  # FFH specifically
  expect_true(all(is.na(results$phenocor$cramerV["smoke", c("group", "grade", "rare_group")])))
  expect_true(all(is.na(results$phenocor$cramerV["sex", "group"])))
  expect_true(all(is.na(results$phenocor$cramerV["group", c("grade", "rare_group")])))
  expect_true(all(is.na(results$phenocor$cramerV["grade", "rare_group"])))
  expect_all_equal(results$phenocor$testType["smoke", c("group", "grade", "rare_group")], "FFH")
  expect_all_equal(results$phenocor$testType["sex", "group"], "FFH")
  expect_all_equal(results$phenocor$testType["group", c("grade", "rare_group")], "FFH")
  expect_all_equal(results$phenocor$testType["grade", "rare_group"], "FFH")
  # Chi-square specifically
  expect_true(all(!is.na(results$phenocor$cramerV["smoke", "sex"])))
  expect_true(all(!is.na(results$phenocor$cramerV["sex", c("grade", "rare_group")])))
  expect_all_equal(results$phenocor$testType["smoke", "sex"], "Chi-square")
  expect_all_equal(results$phenocor$testType["sex", c("grade", "rare_group")], "Chi-square")
  # Ensure that statistics are calculated as expected (ANOVA)
  expect_true(all(!is.na(results$phenocor$stat["group", "height"])))
  expect_true(all(!is.na(results$phenocor$stat["height", "grade"])))
  expect_true(all(is.na(results$phenocor$cramerV["group", "height"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "grade"])))
  expect_all_equal(results$phenocor$testType["group", "height"], "ANOVA")
  expect_all_equal(results$phenocor$testType["height", "grade"], "ANOVA")
  # Ensure that statistics are calculated as expected (t-test).
  expect_true(all(!is.na(results$phenocor$stat[c("smoke", "sex"), "height"])))
  expect_true(all(!is.na(results$phenocor$stat["height", "rare_group"])))
  expect_true(all(is.na(results$phenocor$cramerV[c("smoke", "sex"), "height"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "rare_group"])))
  expect_all_equal(results$phenocor$testType[c("smoke", "sex"), "height"], "T-Test")
  expect_all_equal(results$phenocor$testType["height", "rare_group"], "T-Test")
  # Ensure that statistics are calculated as expected (correlation)
  expect_true(all(!is.na(results$phenocor$stat["height", "WBC"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "WBC"])))
  expect_all_equal(results$phenocor$testType["height", "WBC"], "Cor")
  
  # Run seahorse without correlation matrix
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE)
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height", "group") %in% names(results$GSEA)))
  expect_equal(results$phenotype_association$sex$stat, results$phenotype_association$sex$padj)
  expect_equal(results$phenotype_association$height$padj, rep("NA", length(results$phenotype_association$height$padj)))
  expect_equal(results$phenotype_association$group$stat, results$phenotype_association$group$padj)
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Run seahorse without phenotype matrix
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_phenotype_cor = FALSE)
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height", "group") %in% names(results$GSEA)))
  expect_equal(results$phenotype_association$sex$stat, results$phenotype_association$sex$padj)
  expect_equal(results$phenotype_association$height$padj, rep("NA", length(results$phenotype_association$height$padj)))
  expect_equal(results$phenotype_association$group$stat, results$phenotype_association$group$padj)
  expect_true(is.na(results$phenocor))
  expect_true(all(!is.na(results$coexpression)))
  
  # Run seahorse with bonferroni correction
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, pval_adj_method = "bonferroni")
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height", "group") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sex$stat, method = "bonferroni"), 
               results$phenotype_association$sex$padj)
  expect_equal(results$phenotype_association$height$padj, rep("NA", length(results$phenotype_association$height$padj)))
  expect_equal(stats::p.adjust(results$phenotype_association$group$stat, method = "bonferroni"), 
               results$phenotype_association$group$padj)
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Run seahorse with fdr correction
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, pval_adj_method = "fdr")
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("sex", "height", "group") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sex$stat, method = "fdr"), 
               results$phenotype_association$sex$padj)
  expect_equal(results$phenotype_association$height$padj, rep("NA", length(results$phenotype_association$height$padj)))
  expect_equal(stats::p.adjust(results$phenotype_association$group$stat, method = "fdr"), 
               results$phenotype_association$group$padj)
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  
  # Run SEAHORSE with linear regression.
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height", "group2", "group3") %in% names(results$GSEA)))
  expect_equal(results$phenotype_association$sexmale$stat, results$phenotype_association$sexmale$padj)
  expect_equal(results$phenotype_association$height$stat, results$phenotype_association$height$padj)
  expect_equal(results$phenotype_association$group3$stat, results$phenotype_association$group3$padj)
  expect_equal(results$phenotype_association$group2$stat, results$phenotype_association$group2$padj)
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Run SEAHORSE with linear regression and Bonferroni adjustment
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear",
                      pval_adj_method = "bonferroni")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height", "group3", "group2") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sexmale$stat, method = "bonferroni"),
               results$phenotype_association$sexmale$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$height$stat, method = "bonferroni"),
               results$phenotype_association$height$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$group3$stat, method = "bonferroni"),
               results$phenotype_association$group3$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$group2$stat, method = "bonferroni"),
               results$phenotype_association$group2$padj)
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Run SEAHORSE with linear regression and FDR adjustment
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear", pval_adj_method = "fdr")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height", "group2", "group3") %in% names(results$GSEA)))
  expect_equal(stats::p.adjust(results$phenotype_association$sexmale$stat, method = "fdr"),
               results$phenotype_association$sexmale$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$height$stat, method = "fdr"),
               results$phenotype_association$height$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$group3$stat, method = "fdr"),
               results$phenotype_association$group3$padj)
  expect_equal(stats::p.adjust(results$phenotype_association$group2$stat, method = "fdr"),
               results$phenotype_association$group2$padj)
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Check that SEAHORSE runs with linear regression and a malformed column name.
  phenotype_data_mal <- phenotype_data
  colnames(phenotype_data_mal) <- c("?sex", "height in inches", "123group")
  results <- seahorse(expression_data, phenotype_data_mal, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "X.sexmale", "height.in.inches", "X123group3", "X123group2") %in% names(results$GSEA)))
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% colnames(results$phenocor$padj)))
  
  # Check that SEAHORSE runs with correlation and a malformed column name.
  results <- seahorse(expression_data, phenotype_data_mal, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "pearson")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% names(results$GSEA)))
  expect_true(is.na(results$coexpression))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("X.sex", "height.in.inches", "X123group") %in% colnames(results$phenocor$padj)))
  
  # Check that SEAHORSE runs and prints an appropriate message when NA values are
  # included.
  phenotype_data_na <- phenotype_data
  phenotype_data_na[5,"height"] <- NA
  expect_message(seahorse(expression_data, phenotype_data_na, phenotype_dictionary, pathways,
                          compute_gene_cor = FALSE, assoc_method = "linear"),
                 "Out of 10 we retained 9 samples.")
  results <- seahorse(expression_data, phenotype_data_na, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height", "group3", "group2") %in% names(results$GSEA)))
  expect_true(all(is.na(results$coexpression)))
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Run SEAHORSE with linear regression and the correlation matrix.
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = TRUE, compute_phenotype_cor = FALSE, assoc_method = "linear")
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height", "group3", "group2") %in% names(results$GSEA)))
  expect_true(all(!is.na(results$coexpression)))
})