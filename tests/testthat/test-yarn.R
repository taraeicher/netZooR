context("test YARN helper functions")

# Helper: create a minimal ExpressionSet for testing
make_tiny_eset <- function(ngenes = 20, nsamples = 10) {
  set.seed(42)
  expr_mat <- matrix(rpois(ngenes * nsamples, lambda = 100),
                     nrow = ngenes, ncol = nsamples,
                     dimnames = list(paste0("gene", 1:ngenes),
                                    paste0("sample", 1:nsamples)))
  fdata <- data.frame(
    chromosome_name = rep(c("1", "2", "X", "Y", "MT"), length.out = ngenes),
    gene_biotype = rep(c("protein_coding", "lincRNA"), length.out = ngenes),
    row.names = paste0("gene", 1:ngenes)
  )
  pdata <- data.frame(
    tissue = rep(c("brain", "liver"), each = nsamples / 2),
    row.names = paste0("sample", 1:nsamples)
  )
  Biobase::ExpressionSet(
    assayData = expr_mat,
    phenoData = Biobase::AnnotatedDataFrame(pdata),
    featureData = Biobase::AnnotatedDataFrame(fdata)
  )
}

test_that("filterGenes() removes genes on specified chromosomes", {
  obj <- make_tiny_eset()
  filtered <- filterGenes(obj, labels = c("X", "Y", "MT"),
                          featureName = "chromosome_name")
  remaining <- Biobase::fData(filtered)[, "chromosome_name"]
  expect_true(!any(remaining %in% c("X", "Y", "MT")))
  expect_true(nrow(filtered) < nrow(obj))
})

test_that("filterGenes() keepOnly mode keeps only specified labels", {
  obj <- make_tiny_eset()
  filtered <- filterGenes(obj, labels = "protein_coding",
                          featureName = "gene_biotype", keepOnly = TRUE)
  remaining <- Biobase::fData(filtered)[, "gene_biotype"]
  expect_true(all(remaining == "protein_coding"))
})

test_that("filterMissingGenes() removes zero-sum genes", {
  obj <- make_tiny_eset()
  # Set some genes to all zeros
  Biobase::exprs(obj)[1:3, ] <- 0

  filtered <- filterMissingGenes(obj)
  expect_equal(unname(nrow(filtered)), 17L)
})

test_that("filterMissingGenes() respects threshold", {
  obj <- make_tiny_eset()
  # With a very high threshold, most genes should be filtered
  filtered <- filterMissingGenes(obj, threshold = 1e6)
  expect_true(unname(nrow(filtered)) < unname(nrow(obj)))
})

test_that("filterMissingGenes() keeps all genes when none are missing", {
  obj <- make_tiny_eset()
  filtered <- filterMissingGenes(obj, threshold = 0)
  expect_equal(unname(nrow(filtered)), unname(nrow(obj)))
})

test_that("filterSamples() removes specified samples by name", {
  obj <- make_tiny_eset()
  filtered <- filterSamples(obj, ids = c("sample1", "sample2"))
  expect_equal(unname(ncol(filtered)), 8L)
  expect_true(!any(colnames(filtered) %in% c("sample1", "sample2")))
})

test_that("filterSamples() removes by group label", {
  obj <- make_tiny_eset()
  filtered <- filterSamples(obj, ids = "brain", groups = "tissue")
  remaining <- Biobase::pData(filtered)[, "tissue"]
  expect_true(!any(remaining == "brain"))
})

test_that("filterSamples() keepOnly keeps only specified", {
  obj <- make_tiny_eset()
  filtered <- filterSamples(obj, ids = "brain", groups = "tissue",
                            keepOnly = TRUE)
  remaining <- Biobase::pData(filtered)[, "tissue"]
  expect_true(all(remaining == "brain"))
})

test_that("plotCMDS() returns coordinates without plotting", {
  obj <- make_tiny_eset()
  coords <- plotCMDS(obj, comp = 1:2, plotFlag = FALSE)
  expect_true(is.matrix(coords) || is.data.frame(coords))
  expect_equal(ncol(coords), 2)
  expect_equal(unname(nrow(coords)), unname(ncol(obj)))
})

test_that("plotCMDS() with samples=FALSE returns gene coordinates", {
  obj <- make_tiny_eset()
  coords <- plotCMDS(obj, comp = 1:2, samples = FALSE, plotFlag = FALSE, n = 10)
  expect_true(is.matrix(coords) || is.data.frame(coords))
  expect_equal(ncol(coords), 2)
  expect_equal(unname(nrow(coords)), 10)
})

test_that("plotCMDS() with plotting works", {
  obj <- make_tiny_eset()
  coords <- plotCMDS(obj, comp = 1:2, plotFlag = TRUE)
  expect_true(!is.null(coords))
  graphics.off()
})

test_that("checkMisAnnotation() returns coordinates without error", {
  obj <- make_tiny_eset()
  # Add a sex-linked phenotype for trivial checking
  result <- checkMisAnnotation(obj, phenotype = "tissue",
                               controlGenes = "all", plotFlag = FALSE)
  expect_true(!is.null(result))
})

test_that("checkMisAnnotation() with controlGenes filter works", {
  obj <- make_tiny_eset()
  result <- checkMisAnnotation(obj, phenotype = "tissue",
                               controlGenes = "X",
                               columnID = "chromosome_name",
                               plotFlag = FALSE)
  expect_true(!is.null(result))
})

test_that("checkMisAnnotation() with vector phenotype works", {
  obj <- make_tiny_eset()
  pheno_vec <- factor(rep(c("A", "B"), each = 5))
  result <- checkMisAnnotation(obj, phenotype = pheno_vec,
                               controlGenes = "all", plotFlag = FALSE)
  expect_true(!is.null(result))
})

test_that("checkTissuesToMerge() returns CMDS results", {
  obj <- make_tiny_eset()
  # Add major/minor group columns
  Biobase::pData(obj)$major <- rep(c("brain", "liver"), each = 5)
  Biobase::pData(obj)$minor <- rep(c("cortex", "hippo", "cortex", "left", "right"), 2)
  
  result <- checkTissuesToMerge(obj, majorGroups = "major",
                                minorGroups = "minor", plotFlag = FALSE)
  expect_true(!is.null(result))
  expect_type(result, "list")
  graphics.off()
})

test_that("checkTissuesToMerge() with vector groups works", {
  obj <- make_tiny_eset()
  major <- rep(c("brain", "liver"), each = 5)
  
  result <- checkTissuesToMerge(obj, majorGroups = major,
                                minorGroups = "tissue", plotFlag = FALSE)
  expect_true(!is.null(result))
  graphics.off()
})

test_that("checkTissuesToMerge() with filterFun works", {
  obj <- make_tiny_eset()
  Biobase::pData(obj)$major <- rep(c("brain", "liver"), each = 5)
  
  result <- checkTissuesToMerge(obj, majorGroups = "major",
                                minorGroups = "tissue",
                                filterFun = function(x) filterMissingGenes(x),
                                plotFlag = FALSE)
  expect_true(!is.null(result))
  graphics.off()
})

test_that("extractMatrix() returns raw counts", {
  obj <- make_tiny_eset()
  mat <- extractMatrix(obj, normalized = FALSE, log = FALSE)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), dim(Biobase::exprs(obj)))
  expect_equal(mat, Biobase::exprs(obj))
})

test_that("extractMatrix() returns log-transformed counts", {
  obj <- make_tiny_eset()
  mat <- extractMatrix(obj, normalized = FALSE, log = TRUE)
  expect_true(is.matrix(mat))
  expect_equal(mat, log2(Biobase::exprs(obj) + 1))
})

test_that("extractMatrix() with normalized=TRUE requires normalizedMatrix", {
  obj <- make_tiny_eset()
  expect_error(extractMatrix(obj, normalized = TRUE),
               "normalizedMatrix assayData missing")
})

test_that("extractMatrix() works on plain matrix input", {
  mat <- matrix(1:20, nrow = 4)
  result <- extractMatrix(mat, log = FALSE)
  expect_equal(result, mat)
  
  result_log <- extractMatrix(mat, log = TRUE)
  expect_equal(result_log, log2(mat + 1))
})

test_that("filterLowGenes() removes low-expression genes", {
  skip_if_not_installed("edgeR")
  obj <- make_tiny_eset()
  # Set some genes to very low counts
  Biobase::exprs(obj)[1:5, ] <- 0
  
  filtered <- filterLowGenes(obj, groups = "tissue", threshold = 1)
  expect_true(unname(nrow(filtered)) < unname(nrow(obj)))
})

test_that("filterLowGenes() with vector groups works", {
  skip_if_not_installed("edgeR")
  obj <- make_tiny_eset()
  Biobase::exprs(obj)[1:5, ] <- 0
  groups <- rep(c("A", "B"), each = 5)
  
  filtered <- filterLowGenes(obj, groups = groups, threshold = 1)
  expect_true(unname(nrow(filtered)) <= unname(nrow(obj)))
})

test_that("filterLowGenes() with custom minSamples works", {
  skip_if_not_installed("edgeR")
  obj <- make_tiny_eset()
  
  filtered <- filterLowGenes(obj, groups = "tissue", threshold = 1, minSamples = 2)
  expect_true(unname(nrow(filtered)) <= unname(nrow(obj)))
})

test_that("normalizeTissueAware() with qsmooth method works", {
  obj <- make_tiny_eset()
  
  result <- normalizeTissueAware(obj, groups = "tissue",
                                 normalizationMethod = "qsmooth")
  
  expect_s4_class(result, "ExpressionSet")
  expect_true("normalizedMatrix" %in% names(Biobase::assayData(result)))
  norm_mat <- Biobase::assayData(result)[["normalizedMatrix"]]
  expect_equal(dim(norm_mat), dim(Biobase::exprs(obj)))
})

test_that("normalizeTissueAware() with quantile method works", {
  skip_if_not_installed("preprocessCore")
  obj <- make_tiny_eset()
  
  result <- normalizeTissueAware(obj, groups = "tissue",
                                 normalizationMethod = "quantile")
  
  expect_s4_class(result, "ExpressionSet")
  expect_true("normalizedMatrix" %in% names(Biobase::assayData(result)))
})

test_that("normalizeTissueAware() quantile with single group works", {
  skip_if_not_installed("preprocessCore")
  obj <- make_tiny_eset()
  groups <- rep("single_group", ncol(obj))
  
  result <- normalizeTissueAware(obj, groups = groups,
                                 normalizationMethod = "quantile")
  
  expect_s4_class(result, "ExpressionSet")
  expect_true("normalizedMatrix" %in% names(Biobase::assayData(result)))
})

test_that("normalizeTissueAware() with vector groups works", {
  obj <- make_tiny_eset()
  groups <- factor(rep(c("brain", "liver"), each = 5))
  
  result <- normalizeTissueAware(obj, groups = groups,
                                 normalizationMethod = "qsmooth")
  
  expect_s4_class(result, "ExpressionSet")
})

test_that("extractMatrix() on normalized data works", {
  obj <- make_tiny_eset()
  norm_obj <- normalizeTissueAware(obj, groups = "tissue",
                                   normalizationMethod = "qsmooth")
  
  mat <- extractMatrix(norm_obj, normalized = TRUE)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), dim(Biobase::exprs(obj)))
})

test_that("qsmooth() performs quantile smoothing", {
  obj <- make_tiny_eset()
  groups <- factor(Biobase::pData(obj)$tissue)
  
  result <- qsmooth(obj, groups = groups, log = TRUE)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(Biobase::exprs(obj)))
  expect_equal(rownames(result), rownames(Biobase::exprs(obj)))
  expect_equal(colnames(result), colnames(Biobase::exprs(obj)))
})

test_that("qsmooth() with norm.factors works", {
  obj <- make_tiny_eset()
  groups <- factor(Biobase::pData(obj)$tissue)
  norm.factors <- rep(1, ncol(obj))
  
  result <- qsmooth(obj, groups = groups, norm.factors = norm.factors)
  expect_true(is.matrix(result))
})

test_that("qsmooth() with plot=TRUE works", {
  obj <- make_tiny_eset()
  groups <- factor(Biobase::pData(obj)$tissue)
  
  result <- qsmooth(obj, groups = groups, plot = TRUE)
  expect_true(is.matrix(result))
  graphics.off()
})

test_that("qsmooth() with log=FALSE works", {
  obj <- make_tiny_eset()
  groups <- factor(Biobase::pData(obj)$tissue)
  
  result <- qsmooth(obj, groups = groups, log = FALSE)
  expect_true(is.matrix(result))
})

test_that("qstats() computes quantile statistics", {
  set.seed(42)
  exprs <- matrix(rnorm(200), nrow = 20, ncol = 10)
  groups <- rep(c("A", "B"), each = 5)
  
  result <- qstats(exprs, groups, window = 0.05)
  
  expect_type(result, "list")
  expect_true("Q" %in% names(result))
  expect_true("Qref" %in% names(result))
  expect_true("Qhat" %in% names(result))
  expect_true("smoothWeights" %in% names(result))
  expect_equal(length(result$Qref), nrow(exprs))
})

test_that("plotHeatmap() works on ExpressionSet", {
  skip_if_not_installed("gplots")
  obj <- make_tiny_eset()
  
  result <- plotHeatmap(obj, n = 5, normalized = FALSE, log = TRUE, trace = "none")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5)
  graphics.off()
})

test_that("plotDensity() generates density plot", {
  skip_if_not_installed("quantro")
  obj <- make_tiny_eset()
  
  expect_error(
    plotDensity(obj, groups = "tissue"),
    NA
  )
  graphics.off()
})

test_that("plotDensity() with legendPos works", {
  skip_if_not_installed("quantro")
  obj <- make_tiny_eset()
  
  expect_error(
    plotDensity(obj, groups = "tissue", legendPos = "topleft"),
    NA
  )
  graphics.off()
})
