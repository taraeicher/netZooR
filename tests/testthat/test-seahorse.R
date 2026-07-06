context("test SEAHORSE result")

test_that("seahorse function works", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("matrixTests")
  skip_if_not_installed("stats")
  skip_if_not_installed("limma")
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways))

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
  # We expect that "smoke" will trigger FFH in all except sex because EXPECTED counts < 5,
  # "group" will trigger FFH in "rare group" because EXPECTED counts < 5,
  # and "group" will trigger FFH in "rare group" because EXPECTED counts < 5.
  # and "grade" will trigger FFH in "rare group" because EXPECTED counts < 5.
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
  results <- suppressWarnings(seahorse(expression_data_2, phenotype_data_2, phenotype_dictionary_2, pathways))
  
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
  # FFH specifically
  expect_true(all(is.na(results$phenocor$cramerV["smoke", c("group", "grade", "rare_group")])))
  expect_true(all(is.na(results$phenocor$cramerV["group", c("grade", "rare_group")])))
  expect_true(all(is.na(results$phenocor$cramerV["grade", "rare_group"])))
  expect_equal(results$phenocor$stat["smoke", c("group", "grade", "rare_group")],
               results$phenocor$padj["smoke", c("group", "grade", "rare_group")])
  expect_equal(results$phenocor$stat["sex", "group"],
               results$phenocor$padj["sex", "group"])
  expect_equal(results$phenocor$stat["group", c("grade", "rare_group")],
               results$phenocor$padj["group", c("grade", "rare_group")])
  expect_equal(results$phenocor$stat["grade", "rare_group"],
               results$phenocor$padj["grade", "rare_group"])
  expect_all_equal(results$phenocor$testType["smoke", c("group", "grade", "rare_group")], "FFH")
  expect_all_equal(results$phenocor$testType["group", c("grade", "rare_group")], "FFH")
  expect_all_equal(results$phenocor$testType["grade", "rare_group"], "FFH")
  # Chi-square specifically
  expect_true(all(!is.na(results$phenocor$cramerV["smoke", "sex"])))
  expect_true(all(!is.na(results$phenocor$cramerV["sex", "group"])))
  expect_true(all(!is.na(results$phenocor$cramerV["sex", c("grade", "rare_group")])))
  expect_equal(results$phenocor$stat["smoke", "sex"],
               results$phenocor$padj["smoke", "sex"])
  expect_equal(results$phenocor$stat["sex", c("grade", "rare_group")],
               results$phenocor$padj["sex", c("grade", "rare_group")])
  expect_all_equal(results$phenocor$testType["smoke", "sex"], "Chi-square")
  expect_all_equal(results$phenocor$testType["sex", "group"], "Chi-square")
  expect_all_equal(results$phenocor$testType["sex", c("grade", "rare_group")], "Chi-square")
  # Ensure that statistics are calculated as expected (ANOVA)
  expect_true(all(!is.na(results$phenocor$stat["group", "height"])))
  expect_true(all(!is.na(results$phenocor$stat["height", "grade"])))
  expect_true(all(is.na(results$phenocor$cramerV["group", "height"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "grade"])))
  expect_equal(results$phenocor$stat["group", "height"],
               results$phenocor$padj["group", "height"])
  expect_equal(results$phenocor$stat["height", "grade"],
               results$phenocor$padj["height", "grade"])
  expect_all_equal(results$phenocor$testType["group", "height"], "ANOVA")
  expect_all_equal(results$phenocor$testType["height", "grade"], "ANOVA")
  # Ensure that statistics are calculated as expected (t-test).
  expect_true(all(!is.na(results$phenocor$stat[c("smoke", "sex"), "height"])))
  expect_true(all(!is.na(results$phenocor$stat["height", "rare_group"])))
  expect_true(all(is.na(results$phenocor$cramerV[c("smoke", "sex"), "height"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "rare_group"])))
  expect_equal(results$phenocor$stat[c("smoke", "sex"), "height"],
               results$phenocor$padj[c("smoke", "sex"), "height"])
  expect_equal(results$phenocor$stat["height", "rare_group"],
               results$phenocor$padj["height", "rare_group"])
  expect_all_equal(results$phenocor$testType[c("smoke", "sex"), "height"], "T-Test")
  expect_all_equal(results$phenocor$testType["height", "rare_group"], "T-Test")
  # Ensure that statistics are calculated as expected (correlation)
  expect_true(all(!is.na(results$phenocor$stat["height", "WBC"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "WBC"])))
  expect_true(all(is.na(results$phenocor$padj["height", "WBC"])))
  expect_all_equal(results$phenocor$testType["height", "WBC"], "Cor")
  
  # Verify the phenotype data with FDR adjustment.
  results <- suppressWarnings(seahorse(expression_data_2, phenotype_data_2, phenotype_dictionary_2, pathways,
                      pval_adj_method = "fdr"))
  
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
  # FFH specifically
  expect_true(all(is.na(results$phenocor$cramerV["smoke", c("group", "grade", "rare_group")])))
  expect_true(all(is.na(results$phenocor$cramerV["group", c("grade", "rare_group")])))
  expect_true(all(is.na(results$phenocor$cramerV["grade", "rare_group"])))
  expect_equal(p.adjust(results$phenocor$stat["smoke", c("group", "grade", "rare_group")], method = "fdr"),
               results$phenocor$padj["smoke", c("group", "grade", "rare_group")])
  expect_equal(p.adjust(results$phenocor$stat["sex", "group"], method = "fdr"),
               results$phenocor$padj["sex", "group"])
  expect_equal(p.adjust(results$phenocor$stat["group", c("grade", "rare_group")], method = "fdr"),
               results$phenocor$padj["group", c("grade", "rare_group")])
  expect_equal(p.adjust(results$phenocor$stat["grade", "rare_group"], method = "fdr"),
               results$phenocor$padj["grade", "rare_group"])
  expect_all_equal(results$phenocor$testType["smoke", c("group", "grade", "rare_group")], "FFH")
  expect_all_equal(results$phenocor$testType["group", c("grade", "rare_group")], "FFH")
  expect_all_equal(results$phenocor$testType["grade", "rare_group"], "FFH")
  # Chi-square specifically
  expect_true(all(!is.na(results$phenocor$cramerV["smoke", "sex"])))
  expect_all_equal(results$phenocor$testType["sex", "group"], "Chi-square")
  expect_true(all(!is.na(results$phenocor$cramerV["sex", "group"])))
  expect_true(all(!is.na(results$phenocor$cramerV["sex", c("grade", "rare_group")])))
  expect_equal(p.adjust(results$phenocor$stat["smoke", "sex"], method = "fdr"),
               results$phenocor$padj["smoke", "sex"])
  expect_equal(p.adjust(results$phenocor$stat["sex", c("group", "grade", "rare_group")], method = "fdr"),
               results$phenocor$padj["sex", c("group", "grade", "rare_group")])
  expect_all_equal(results$phenocor$testType["smoke", "sex"], "Chi-square")
  expect_all_equal(results$phenocor$testType["sex", c("grade", "rare_group")], "Chi-square")
  # Ensure that statistics are calculated as expected (ANOVA)
  expect_true(all(!is.na(results$phenocor$stat["group", "height"])))
  expect_true(all(!is.na(results$phenocor$stat["height", "grade"])))
  expect_true(all(is.na(results$phenocor$cramerV["group", "height"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "grade"])))
  expect_equal(p.adjust(results$phenocor$stat["group", c("height", "WBC")], method = "fdr"),
               results$phenocor$padj["group", c("height", "WBC")])
  expect_equal(p.adjust(results$phenocor$stat["height", "grade"], method = "fdr"),
               results$phenocor$padj["height", "grade"])
  expect_all_equal(results$phenocor$testType["group", "height"], "ANOVA")
  expect_all_equal(results$phenocor$testType["height", "grade"], "ANOVA")
  # Ensure that statistics are calculated as expected (t-test).
  expect_true(all(!is.na(results$phenocor$stat[c("smoke", "sex"), "height"])))
  expect_true(all(!is.na(results$phenocor$stat["height", "rare_group"])))
  expect_true(all(is.na(results$phenocor$cramerV[c("smoke", "sex"), "height"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "rare_group"])))
  expect_equal(p.adjust(results$phenocor$stat["smoke", c("height", "WBC")], method = "fdr"),
               results$phenocor$padj["smoke", c("height", "WBC")])
  expect_equal(p.adjust(results$phenocor$stat["sex", c("height", "WBC")], method = "fdr"),
               results$phenocor$padj["sex", c("height", "WBC")])
  expect_equal(p.adjust(results$phenocor$stat["rare_group", "WBC"], method = "fdr"),
               results$phenocor$padj["rare_group", "WBC"])
  expect_all_equal(results$phenocor$testType[c("smoke", "sex"), "height"], "T-Test")
  expect_all_equal(results$phenocor$testType["height", "rare_group"], "T-Test")
  # Ensure that statistics are calculated as expected (correlation)
  expect_true(all(!is.na(results$phenocor$stat["height", "WBC"])))
  expect_true(all(is.na(results$phenocor$cramerV["height", "WBC"])))
  expect_true(all(is.na(results$phenocor$padj["height", "WBC"])))
  expect_all_equal(results$phenocor$testType["height", "WBC"], "Cor")
  
  # Run seahorse without correlation matrix
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE))
  
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_phenotype_cor = FALSE))
  
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
  
  # Run seahorse without gene-phenotype matrix
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_phenotype_cor = FALSE))
  
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  expect_true(all(!is.na(results$phenocor)))
  expect_true(all(!is.na(results$coexpression)))
  expect_equal(length(results$phenotype_association), 0)
  expect_equal(length(results$GSEA), 0)
  expect_true(all(c("stat", "cramerV", "padj", "testType") %in% names(results$phenocor)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$stat)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$cramerV)))
  expect_true(all(c("sex", "height", "group") %in% rownames(results$phenocor$padj)))
  expect_true(all(c("sex", "height", "group") %in% colnames(results$phenocor$padj)))
  
  # Run seahorse with bonferroni correction
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, pval_adj_method = "bonferroni"))
  
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, pval_adj_method = "fdr"))
  
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear"))
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear",
                      pval_adj_method = "bonferroni"))
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear", pval_adj_method = "fdr"))
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data_mal, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear"))
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data_mal, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "pearson"))
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
  expect_message(suppressWarnings(seahorse(expression_data, phenotype_data_na, phenotype_dictionary, pathways,
                          compute_gene_cor = FALSE, assoc_method = "linear")),
                 "Out of 10 we retained 9 samples.")
  results <- suppressWarnings(seahorse(expression_data, phenotype_data_na, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, assoc_method = "linear"))
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
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = TRUE, compute_phenotype_cor = FALSE, assoc_method = "linear"))
  # Verify structure
  expect_type(results, "list")
  expect_true(length(results) > 0)
  # Check that results contain expected top-level keys
  expect_true(all(c("coexpression", "phenotype_association", "GSEA", "phenocor") %in% names(results)))
  # Check that phenotype names appear in sub-lists
  expect_true(all(c("(Intercept)", "sexmale", "height", "group3", "group2") %in% names(results$GSEA)))
  expect_true(all(!is.na(results$coexpression)))
  
  # Run SEAHORSE and save usage file (all).
  skip_if_not_installed("peakRAM")
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                      usage_report_file = "usage_file.RDS"))
  usage <- readRDS("usage_file.RDS")
  expect_true(all(c("GeneToGeneCor", "PhenToPhenCor", "PhenToGeneCor") %in% names(usage)))
  expect_true(length(colnames(usage$GeneToGeneCor)) == 4)
  expect_true(length(colnames(usage$PhenToPhenCor)) == 4)
  expect_true(length(colnames(usage$PhenToGeneCor)) == 4)
  unlink("usage_file.RDS")
  
  # Input invalid usage file.
  expect_error(suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                      usage_report_file = "/some_nonexistent_dir/usage_file.RDS")),
               "/some_nonexistent_dir/usage_file.RDS could not be created.")
  
  # Save usage file (all but gene-gene)
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                      usage_report_file = "usage_file.RDS"))
  usage <- readRDS("usage_file.RDS")
  expect_true(all(c("GeneToGeneCor", "PhenToPhenCor", "PhenToGeneCor") %in% names(usage)))
  expect_true(is.na(usage$GeneToGeneCor))
  expect_true(length(colnames(usage$PhenToPhenCor)) == 4)
  expect_true(length(colnames(usage$PhenToGeneCor)) == 4)
  unlink("usage_file.RDS")
  
  # Save usage file (all but phenotype-phenotype)
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = TRUE, compute_phenotype_cor = FALSE, assoc_method = "linear",
                      usage_report_file = "usage_file.RDS"))
  usage <- readRDS("usage_file.RDS")
  expect_true(all(c("GeneToGeneCor", "PhenToPhenCor", "PhenToGeneCor") %in% names(usage)))
  expect_true(is.na(usage$PhenToPhenCor))
  expect_true(length(colnames(usage$GeneToGeneCor)) == 4)
  expect_true(length(colnames(usage$PhenToGeneCor)) == 4)
  unlink("usage_file.RDS")
  
  # Save usage file (only phenotype-genotype)
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = FALSE, compute_phenotype_cor = FALSE, assoc_method = "linear",
                      usage_report_file = "usage_file.RDS"))
  usage <- readRDS("usage_file.RDS")
  expect_true(all(c("GeneToGeneCor", "PhenToPhenCor", "PhenToGeneCor") %in% names(usage)))
  expect_true(is.na(usage$PhenToPhenCor))
  expect_true(is.na(usage$GeneToGeneCor))
  expect_true(length(colnames(usage$PhenToGeneCor)) == 4)
  unlink("usage_file.RDS")
  
  # Test that statements print when verbose = TRUE.
  expect_output(suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                      compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                      verbose = TRUE)),
                "Running gene co-expression")
  expect_output(suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                         compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                         verbose = TRUE)),
                "Gene co-expression complete")
  expect_output(suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                         compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                         verbose = TRUE)),
                "Running phenotype associations")
  expect_output(suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                         compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                         verbose = TRUE)),
                "Phenotype associations complete")
  expect_output(suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                         compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                         verbose = TRUE)),
                "Running gene-phenotype associations")
  expect_output(suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways,
                         compute_gene_cor = TRUE, compute_phenotype_cor = TRUE, assoc_method = "linear",
                         verbose = TRUE)),
                "Gene-phenotype associations complete")
  
  # Check that NA values are returned when the phenotype has less than 2 levels.
  phenotype_data_1_lev <- phenotype_data
  phenotype_data_1_lev$sex <- rep("male", nrow(phenotype_data_1_lev))
  result <- suppressWarnings(seahorse(expression_data, phenotype_data_1_lev, phenotype_dictionary, pathways,
           compute_gene_cor = FALSE, compute_phenotype_cor = FALSE))
  expect_equal(result$phenotype_association$sex$stat, NA)
  expect_equal(result$GSEA$sex, NA)
  
  # Check that GSEA returns an empty list if it cannot run, e.g. if there is no pathway overlap.
  uselessPathways <- list(pathway1 = paste(rep(c("Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Neptune"), 2),
                                           c(rep(1, 8), rep(2, 8))),
                          pathway2 = paste(rep(c("Nina", "Pinta", "Santa Maria"), 6),
                                           c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3))),
                          pathway3 = c("Neutron Star", "Supernova", "Black Hole", "Red Dwarf", "White Dwarf", "Red Giant",
                                       "Blue Giant", "Nebula", "Oort Cloud", "Brown Dwarf", "Comet", "Asteroid", "Dark Matter",
                                       "Dark Energy", "Binary Star System", "Black Dwarf"))
  result <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, uselessPathways,
                     compute_gene_cor = FALSE, compute_phenotype_cor = FALSE))
  expect_equal(nrow(result$GSEA$sex), 0)
  expect_equal(nrow(result$GSEA$height), 0)
  expect_equal(nrow(result$GSEA$group), 0)
  
  # Check that a phenotype pair result will be set to NA if we have 2 or more zero marginals.
  phenotype_data_2_zeros <- phenotype_data_2
  phenotype_data_2_zeros[which(phenotype_data_2$grade == "A")[1:10], "sex"] <- "female"
  phenotype_data_2_zeros[which(phenotype_data_2$grade == "A")[2:20], "sex"] <- "male"
  phenotype_data_2_zeros[which(phenotype_data_2$grade == "A")[21:45], "sex"] <- NA
  phenotype_data_2_zeros[which(phenotype_data_2$grade == "B"), "sex"] <- NA
  phenotype_data_2_zeros[which(phenotype_data_2$grade == "C"), "sex"] <- NA
  result <- suppressWarnings(seahorse(expression_data_2, phenotype_data_2_zeros, phenotype_dictionary_2, pathways,
                     compute_gene_cor = FALSE, compute_phenotype_cor = TRUE))
  expect_true(is.na(result$phenocor$stat["sex", "grade"]))
  
  # Check that a phenotype pair result will be set to NA if we have all continuous data missing
  # for all but one level (nominal).
  phenotype_data_2_na <- phenotype_data_2
  phenotype_data_2_na[which(phenotype_data_2$grade == "A"), "WBC"] <- NA
  phenotype_data_2_na[which(phenotype_data_2$grade == "B"), "WBC"] <- NA
  result <- suppressWarnings(seahorse(expression_data_2, phenotype_data_2_na, phenotype_dictionary_2, pathways,
                     compute_gene_cor = FALSE, compute_phenotype_cor = TRUE))
  expect_true(is.na(result$phenocor$stat["grade", "WBC"]))
  
  # Check that a phenotype pair result will be set to NA if we have all continuous data missing
  # for one level (dichotomous).
  phenotype_data_2_na <- phenotype_data_2
  phenotype_data_2_na[which(phenotype_data_2$sex == "female"), "WBC"] <- NA
  result <- suppressWarnings(seahorse(expression_data_2, phenotype_data_2_na, phenotype_dictionary_2, pathways,
                     compute_gene_cor = FALSE, compute_phenotype_cor = TRUE))
  expect_true(is.na(result$phenocor$stat["sex", "WBC"]))
  
  # Check that a gene-phenotype association will be set to NA if there is no variance in either
  # phenotype level.
  phenotype_data_smaller <- phenotype_data_2[,c("sex", "WBC")]
  phenotype_dictionary_smaller <- c("dichotomous", "continuous")
  phenotype_data_smaller[which(phenotype_data_smaller$sex == "female"), "WBC"] <- 1
  phenotype_data_smaller[which(phenotype_data_smaller$sex == "male"), "WBC"] <- 1
  result <- suppressWarnings(seahorse(expression_data_2, phenotype_data_smaller, phenotype_dictionary_smaller, pathways,
                     compute_gene_cor = FALSE, compute_phenotype_cor = TRUE))
  expect_true(is.na(result$phenocor$stat["sex", "WBC"]))
})
test_that(".tsv.gz output works", {
  # Check packages.
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  # Simulate expression data
  expression_data = data.frame(matrix(rexp(1000, rate=.1), ncol=10, nrow = 100))
  rownames(expression_data) = sprintf("ENSG%011d", 1:100)
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
  
  # Make input data and write it.
  dir.create("~/tmpInputDir")
  input <- list(expression = expression_data, phenotype = phenotype_data, dict = phenotype_dictionary)
  saveRDS(input, "~/tmpInputDir/tissue1.RDS")
  saveRDS(input, "~/tmpInputDir/tissue2.RDS")
  saveRDS(input, "~/tmpInputDir/tissue3.RDS")
  saveRDS(c(1, 2, 3), "~/tmpInputDir/badInput.RDS")
  
  # Get toy SEAHORSE results and write them.
  dir.create("~/tmpResultDir")
  results <- suppressWarnings(seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways))
  saveRDS(results, "~/tmpResultDir/tissue1.RDS")
  saveRDS(results, "~/tmpResultDir/tissue2.RDS")
  saveRDS(results, "~/tmpResultDir/tissue3.RDS")
  saveRDS(c(1, 2, 3), "~/tmpResultDir/badInput.RDS")
  
  # Make the data dictionary.
  dataDict <- data.frame(VARNAME = colnames(phenotype_data),
                         VARDESC = c("Sex (male or female)",
                                     "Height in inches",
                                     "Group membership"),
                         VARMETA = rep("phenotype", 3),
                         TYPE = c("string", "decimal", "integer, encoded value"))
  dataDictCols <- c("VARNAME", "VARDESC", "VARMETA", "TYPE")
  
  # Create the output directory.
  outputDir <- "~/tmpOutputDir/"
  dir.create(outputDir)
  
  # Check that an error is thrown if the SEAHORSE input object is formatted incorrectly.
  expect_error(suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                   result_directory = "~/tmpResultDir",
                                   output_directory = outputDir,
                                   data_dictionary = dataDict,
                                   pathways = pathways)),
               "Could not format expression")
  unlink("~/tmpInputDir/badInput.RDS")
  unlink("~/tmpResultDir/badInput.RDS")
  
  # Check that an error is thrown if the SEAHORSE result object is formatted incorrectly.
  saveRDS(input, "~/tmpInputDir/badInput.RDS")
  saveRDS(c(1, 2, 3), "~/tmpResultDir/badInput.RDS")
  expect_error(suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                   result_directory = "~/tmpResultDir",
                                   output_directory = outputDir,
                                   data_dictionary = dataDict,
                                   pathways = pathways)),
               "Could not format pathway enrichment")
  unlink("~/tmpInputDir/badInput.RDS")
  unlink("~/tmpResultDir/badInput.RDS")
  
  # Check that an error is thrown if the data dictionary is formatted incorrectly.
  badDataDict <- dataDict[,c(1:2)]
  expect_error(suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                   result_directory = "~/tmpResultDir",
                                   output_directory = outputDir,
                                   data_dictionary = badDataDict,
                                   pathways = pathways)),
               "Incorrect format for data dictionary")

  # Check that an error is thrown if the number of input and result files do not match.
  saveRDS(results, "~/tmpInputDir/tissue4.RDS")
  expect_error(suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                   result_directory = "~/tmpResultDir",
                                   output_directory = outputDir,
                                   data_dictionary = dataDict,
                                   pathways = pathways)),
               "Number of result and input files do not match")
  unlink("~/tmpInputDir/tissue4.RDS")
  
  # Check that an error is thrown if the input and result file names do not match.
  saveRDS(results, "~/tmpInputDir/badName.RDS")
  unlink("~/tmpInputDir/tissue2.RDS")
  expect_error(suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                   result_directory = "~/tmpResultDir",
                                   output_directory = outputDir,
                                   data_dictionary = dataDict,
                                   pathways = pathways)),
               "Result and input file names do not match")
  saveRDS(input, "~/tmpInputDir/tissue2.RDS")
  unlink("~/tmpInputDir/badName.RDS")
  
  # Check that the files exist and are the correct format when we read them.
  resultFormat <- suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                result_directory = "~/tmpResultDir",
                                output_directory = outputDir,
                                data_dictionary = dataDict,
                                pathways = pathways))
  readFile <- function(fname){
    con <- gzfile(fname, "rt")
    data <- read.delim(con, sep = "\t", header = TRUE)
    return(data)
  }
  # Map genes.
  mappedSymbols <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(expression_data),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first" # Returns one value for each ENSEMBL ID.
  )
  mappedEntrez <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(expression_data),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first" # Returns one value for each ENSEMBL ID.
  )
  mappingResults <- data.frame(ENSEMBL = rownames(expression_data), SYMBOL = mappedSymbols, ENTREZID = mappedEntrez)
  mappingResults$ALIAS <- mappingResults$SYMBOL
  mappingResults <- mappingResults[,c("ALIAS", "ENSEMBL", "SYMBOL", "ENTREZID")]
  mappingResults$ENTREZID <- as.numeric(mappingResults$ENTREZID)
  rownames(mappingResults) <- NULL
  
  expect_equal(colnames(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               c("pathway", "pval","padj",	"log2err",	"ES",	"NES",	"size",	"ranks",	"leadingEdge",	"varname",	"tissue"))
  expect_equal(colnames(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               c("VARNAME","VARDESC",	"VARMETA",	"TYPE"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               c("Gene.A","Gene.B",	"Tissue",	"Correlation"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               c("ENSG","SAMPID",	"GENE_EXPRESSION", "tissue"))
  expect_equal(readFile(paste0(outputDir, "human_ensembl2symbol_map.tsv.gz")),
               mappingResults)
  expect_equal(colnames(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               c("SAMPID","tissue",	"VARNAME", "VALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               c("VARNAME","GENE",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               c("VARNAME1","VARNAME2",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(nrow(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               3 * ncol(phenotype_data) * length(pathways))
  expect_equal(nrow(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               nrow(dataDict))
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               (nrow(expression_data) * (nrow(expression_data) - 1)) / 2 * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               nrow(expression_data) * ncol(expression_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               nrow(phenotype_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               nrow(expression_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               (ncol(phenotype_data) * (ncol(phenotype_data) - 1)) / 2 * 3)
  
  # Check that the output directory is created if we didn't create it.
  outputDir = "~/newOutputDir/"
  resultFormat <- suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                result_directory = "~/tmpResultDir",
                                output_directory = outputDir,
                                data_dictionary = dataDict,
                                pathways = pathways))
  
  expect_equal(colnames(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               c("pathway", "pval","padj",	"log2err",	"ES",	"NES",	"size",	"ranks",	"leadingEdge",	"varname",	"tissue"))
  expect_equal(colnames(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               c("VARNAME","VARDESC",	"VARMETA",	"TYPE"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               c("Gene.A","Gene.B",	"Tissue",	"Correlation"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               c("ENSG","SAMPID",	"GENE_EXPRESSION", "tissue"))
  expect_equal(readFile(paste0(outputDir, "human_ensembl2symbol_map.tsv.gz")),
               mappingResults)
  expect_equal(colnames(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               c("SAMPID","tissue",	"VARNAME", "VALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               c("VARNAME","GENE",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               c("VARNAME1","VARNAME2",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(nrow(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               3 * ncol(phenotype_data) * length(pathways))
  expect_equal(nrow(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               nrow(dataDict))
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               (nrow(expression_data) * (nrow(expression_data) - 1)) / 2 * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               nrow(expression_data) * ncol(expression_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               nrow(phenotype_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               nrow(expression_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               (ncol(phenotype_data) * (ncol(phenotype_data) - 1)) / 2 * 3)
  unlink("~/newOutputDir", recursive = TRUE)
  
  # Check that function still works if we did not infer gene-gene associations.
  outputDir <- "~/tmpOutputDir/"
  resultsNoGene <- results
  resultsNoGene$coexpression <- NA
  saveRDS(resultsNoGene, "~/tmpResultDir/tissue1.RDS")
  saveRDS(resultsNoGene, "~/tmpResultDir/tissue2.RDS")
  saveRDS(resultsNoGene, "~/tmpResultDir/tissue3.RDS")

  resultFormat <- suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                result_directory = "~/tmpResultDir",
                                output_directory = outputDir,
                                data_dictionary = dataDict,
                                pathways = pathways))
  expect_equal(colnames(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               c("pathway", "pval","padj",	"log2err",	"ES",	"NES",	"size",	"ranks",	"leadingEdge",	"varname",	"tissue"))
  expect_equal(colnames(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               c("VARNAME","VARDESC",	"VARMETA",	"TYPE"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               c("Gene.A","Gene.B",	"Tissue",	"Correlation"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               c("ENSG","SAMPID",	"GENE_EXPRESSION", "tissue"))
  expect_equal(readFile(paste0(outputDir, "human_ensembl2symbol_map.tsv.gz")),
               mappingResults)
  expect_equal(colnames(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               c("SAMPID","tissue",	"VARNAME", "VALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               c("VARNAME","GENE",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               c("VARNAME1","VARNAME2",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(nrow(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               3 * ncol(phenotype_data) * length(pathways))
  expect_equal(nrow(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               nrow(dataDict))
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               0)
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               nrow(expression_data) * ncol(expression_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               nrow(phenotype_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               nrow(expression_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               (ncol(phenotype_data) * (ncol(phenotype_data) - 1)) / 2 * 3)
  
  # Check that function still works if we did not infer gene-phenotype associations.
  resultsNoGenePhen <- results
  resultsNoGenePhen$phenotype_association <- list()
  resultsNoGenePhen$GSEA <- list()
  saveRDS(resultsNoGenePhen, "~/tmpResultDir/tissue1.RDS")
  saveRDS(resultsNoGenePhen, "~/tmpResultDir/tissue2.RDS")
  saveRDS(resultsNoGenePhen, "~/tmpResultDir/tissue3.RDS")
  
  resultFormat <- suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                      result_directory = "~/tmpResultDir",
                                      output_directory = outputDir,
                                      data_dictionary = dataDict,
                                      pathways = pathways))
  expect_equal(colnames(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               c("pathway", "pval","padj",	"log2err",	"ES",	"NES",	"size",	"ranks",	"leadingEdge",	"varname",	"tissue"))
  expect_equal(colnames(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               c("VARNAME","VARDESC",	"VARMETA",	"TYPE"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               c("Gene.A","Gene.B",	"Tissue",	"Correlation"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               c("ENSG","SAMPID",	"GENE_EXPRESSION", "tissue"))
  expect_equal(readFile(paste0(outputDir, "human_ensembl2symbol_map.tsv.gz")),
               mappingResults)
  expect_equal(colnames(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               c("SAMPID","tissue",	"VARNAME", "VALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               c("VARNAME","GENE",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               c("VARNAME1","VARNAME2",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(nrow(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               0)
  expect_equal(nrow(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               nrow(dataDict))
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               (nrow(expression_data) * (nrow(expression_data) - 1)) / 2 * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               nrow(expression_data) * ncol(expression_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               nrow(phenotype_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               0)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               (ncol(phenotype_data) * (ncol(phenotype_data) - 1)) / 2 * 3)
  
  # Check that function still works if we did not infer phenotype-phenotype associations.
  resultsNoPhen <- results
  resultsNoPhen$phenocor <- NA
  saveRDS(resultsNoPhen, "~/tmpResultDir/tissue1.RDS")
  saveRDS(resultsNoPhen, "~/tmpResultDir/tissue2.RDS")
  saveRDS(resultsNoPhen, "~/tmpResultDir/tissue3.RDS")
  
  resultFormat <- suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                      result_directory = "~/tmpResultDir",
                                      output_directory = outputDir,
                                      data_dictionary = dataDict,
                                      pathways = pathways))
  expect_equal(colnames(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               c("pathway", "pval","padj",	"log2err",	"ES",	"NES",	"size",	"ranks",	"leadingEdge",	"varname",	"tissue"))
  expect_equal(colnames(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               c("VARNAME","VARDESC",	"VARMETA",	"TYPE"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               c("Gene.A","Gene.B",	"Tissue",	"Correlation"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               c("ENSG","SAMPID",	"GENE_EXPRESSION", "tissue"))
  expect_equal(readFile(paste0(outputDir, "human_ensembl2symbol_map.tsv.gz")),
               mappingResults)
  expect_equal(colnames(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               c("SAMPID","tissue",	"VARNAME", "VALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               c("VARNAME","GENE",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               c("VARNAME1","VARNAME2",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(nrow(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               3 * ncol(phenotype_data) * length(pathways))
  expect_equal(nrow(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               nrow(dataDict))
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               (nrow(expression_data) * (nrow(expression_data) - 1)) / 2 * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               nrow(expression_data) * ncol(expression_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               nrow(phenotype_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               nrow(expression_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               0)
  
  # Check that function still works if we did not infer anything.
  resultsNothing <- results
  resultsNothing$coexpression <- NA
  resultsNothing$phenocor <- NA
  resultsNothing$phenotype_association <- list()
  resultsNothing$GSEA <- list()
  saveRDS(resultsNothing, "~/tmpResultDir/tissue1.RDS")
  saveRDS(resultsNothing, "~/tmpResultDir/tissue2.RDS")
  saveRDS(resultsNothing, "~/tmpResultDir/tissue3.RDS")
  
  resultFormat <- suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                      result_directory = "~/tmpResultDir",
                                      output_directory = outputDir,
                                      data_dictionary = dataDict,
                                      pathways = pathways))
  expect_equal(colnames(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               c("pathway", "pval","padj",	"log2err",	"ES",	"NES",	"size",	"ranks",	"leadingEdge",	"varname",	"tissue"))
  expect_equal(colnames(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               c("VARNAME","VARDESC",	"VARMETA",	"TYPE"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               c("Gene.A","Gene.B",	"Tissue",	"Correlation"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               c("ENSG","SAMPID",	"GENE_EXPRESSION", "tissue"))
  expect_equal(readFile(paste0(outputDir, "human_ensembl2symbol_map.tsv.gz")),
               mappingResults)
  expect_equal(colnames(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               c("SAMPID","tissue",	"VARNAME", "VALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               c("VARNAME","GENE",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               c("VARNAME1","VARNAME2",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(nrow(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               0)
  expect_equal(nrow(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               nrow(dataDict))
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               0)
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               nrow(expression_data) * ncol(expression_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               nrow(phenotype_data) * ncol(phenotype_data) * 3)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               0)
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               0)
  
  # Check that function still works if there is only one tissue.
  saveRDS(results, "~/tmpResultDir/tissue1.RDS")
  unlink("~/tmpResultDir/tissue2.RDS")
  unlink("~/tmpResultDir/tissue3.RDS")
  unlink("~/tmpInputDir/tissue2.RDS")
  unlink("~/tmpInputDir/tissue3.RDS")
  
  resultFormat <- suppressWarnings(seahorseFormatForUI(input_directory = "~/tmpInputDir",
                                      result_directory = "~/tmpResultDir",
                                      output_directory = outputDir,
                                      data_dictionary = dataDict,
                                      pathways = pathways))
  expect_equal(colnames(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               c("pathway", "pval","padj",	"log2err",	"ES",	"NES",	"size",	"ranks",	"leadingEdge",	"varname",	"tissue"))
  expect_equal(colnames(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               c("VARNAME","VARDESC",	"VARMETA",	"TYPE"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               c("Gene.A","Gene.B",	"Tissue",	"Correlation"))
  expect_equal(colnames(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               c("ENSG","SAMPID",	"GENE_EXPRESSION", "tissue"))
  expect_equal(readFile(paste0(outputDir, "human_ensembl2symbol_map.tsv.gz")),
               mappingResults)
  expect_equal(colnames(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               c("SAMPID","tissue",	"VARNAME", "VALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               c("VARNAME","GENE",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(colnames(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               c("VARNAME1","VARNAME2",	"tissue", "TEST", "TESTSTAT", "TESTPVALUE"))
  expect_equal(nrow(readFile(paste0(outputDir, "all_gsea_results.tsv.gz"))),
               ncol(phenotype_data) * length(pathways))
  expect_equal(nrow(readFile(paste0(outputDir, "data_dictionary.tsv.gz"))),
               nrow(dataDict))
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression2geneexpression.tsv.gz"))),
               (nrow(expression_data) * (nrow(expression_data) - 1)) / 2)
  expect_equal(nrow(readFile(paste0(outputDir, "geneexpression_data.tsv.gz"))),
               nrow(expression_data) * ncol(expression_data))
  expect_equal(nrow(readFile(paste0(outputDir, "metadata.tsv.gz"))),
               nrow(phenotype_data) * ncol(phenotype_data))
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2expression.tsv.gz"))),
               nrow(expression_data) * ncol(phenotype_data))
  expect_equal(nrow(readFile(paste0(outputDir, "metadata2metadata.tsv.gz"))),
               (ncol(phenotype_data) * (ncol(phenotype_data) - 1)) / 2)
  
  # Remove files.
  unlink("~/tmpInputDir", recursive = TRUE)
  unlink("~/tmpResultDir", recursive = TRUE)
  unlink("~/tmpOutputDir", recursive = TRUE)
})