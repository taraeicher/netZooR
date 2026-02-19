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
  expect_equal(nrow(filtered), 17)
})

test_that("filterMissingGenes() respects threshold", {
  obj <- make_tiny_eset()
  # With a very high threshold, most genes should be filtered
  filtered <- filterMissingGenes(obj, threshold = 1e6)
  expect_true(nrow(filtered) < nrow(obj))
})

test_that("filterSamples() removes specified samples by name", {
  obj <- make_tiny_eset()
  filtered <- filterSamples(obj, ids = c("sample1", "sample2"))
  expect_equal(ncol(filtered), 8)
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
  expect_equal(nrow(coords), ncol(obj))
})

test_that("checkMisAnnotation() returns coordinates without error", {
  obj <- make_tiny_eset()
  # Add a sex-linked phenotype for trivial checking
  result <- checkMisAnnotation(obj, phenotype = "tissue",
                               controlGenes = "all", plotFlag = FALSE)
  expect_true(!is.null(result))
})
