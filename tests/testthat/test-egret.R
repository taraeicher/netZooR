test_that("EGRET function works", {
  skip_if_offline()
  skip_on_cran()
  
  # Use package bundled test data
  toy_qbic_path <- system.file("extdata", "toy_qbic.txt", package = "netZooR", mustWork = TRUE)
  toy_genotype_path <- system.file("extdata", "toy_genotype.vcf", package = "netZooR", mustWork = TRUE)
  toy_motif_path <- system.file("extdata", "toy_motif_prior.txt", package = "netZooR", mustWork = TRUE)
  toy_expr_path <- system.file("extdata", "toy_expr.txt", package = "netZooR", mustWork = TRUE)
  toy_ppi_path <- system.file("extdata", "toy_ppi_prior.txt", package = "netZooR", mustWork = TRUE)
  toy_eqtl_path <- system.file("extdata", "toy_eQTL.txt", package = "netZooR", mustWork = TRUE)
  toy_map_path <- system.file("extdata", "toy_map.txt", package = "netZooR", mustWork = TRUE)
  
  qbic <- read.table(file = toy_qbic_path, header = FALSE)
  vcf <- read.table(toy_genotype_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                    colClasses = c("character", "numeric", "character", "character", "character",
                                   "character", "character", "character", "character", "character"))
  motif <- read.table(toy_motif_path, sep = "\t", header = FALSE)
  expr <- read.table(toy_expr_path, header = FALSE, sep = "\t", row.names = 1)
  ppi <- read.table(toy_ppi_path, header = FALSE, sep = "\t")
  qtl <- read.table(toy_eqtl_path, header = FALSE)
  nameGeneMap <- read.table(toy_map_path, header = FALSE)
  tag <- "my_toy_egret_run"
  
  # Run in temp directory to avoid polluting working dir
  old_wd <- setwd(tempdir())
  on.exit(setwd(old_wd), add = TRUE)
  
  runEgret(qtl, vcf, qbic, motif, expr, ppi, nameGeneMap, tag)
  
  # Check output files were created
  expect_true(file.exists(paste0(tag, "_egret.RData")))
  expect_true(file.exists(paste0(tag, "_panda.RData")))
  
  # Load and verify output
  load(paste0(tag, "_egret.RData"))
  load(paste0(tag, "_panda.RData"))
  expect_true(is.matrix(regnetE))
  expect_true(is.matrix(regnetP))
  expect_equal(dim(regnetE), dim(regnetP))
  expect_true(all(is.finite(regnetE)))
  expect_true(all(is.finite(regnetP)))
  
  # Cleanup
  file.remove(paste0(tag, "_egret.RData"))
  file.remove(paste0(tag, "_panda.RData"))
  file.remove(paste0("priors_", tag, ".txt"))
})