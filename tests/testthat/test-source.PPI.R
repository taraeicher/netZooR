context("test source PPI")

test_that("sourcePPI works", {
  skip_if_offline()
  skip_if_not_installed("STRINGdb")
  suppressWarnings(load("./testDataset.RData"))
  tf <- data.frame(matrix(c("Rv0022c","Rv0023","Rv0042c","Rv0043c","Rv0047c","Rv0054"), 
                          nrow=6, byrow=T),stringsAsFactors=FALSE)
  # STRINGdb Version 10
  if(R.Version()$major=="3"){
    actual_PPI_V10 <- sourcePPI(tf,"10",83332)
    ppiV10$from = as.factor(ppiV10$from)
    ppiV10$to = as.factor(ppiV10$to)
    expect_equal(actual_PPI_V10, ppiV10)
  }
  # STRINGdb Version 11
  else if(R.Version()$major=="4"){
    #actual_PPI_V11 <- sourcePPI(tf,"11",83332)
    string_db=STRINGdb::STRINGdb$new(version="11", species=83332)
    # change the colname to "TF"
    colnames(tf) <- "TF"
    # map the TF to STRINGdb dataset
    TF_mapped <-  string_db$map(tf,"TF",removeUnmappedRows=FALSE)
    # collect the interactions between the TF of interest
    ppi_tmp <- string_db$get_interactions(TF_mapped$STRING_id)[,c(1,2)]
    # store the PPI by using original identifier.
    actual_PPI_V11 <- data.frame(from=TF_mapped[match(ppi_tmp$from,TF_mapped$STRING_id),1], to=TF_mapped[match(ppi_tmp$to,TF_mapped$STRING_id),1])
    # assign "score"column  to "1"
    actual_PPI_V11$score <- "1"
    #expect_equal(actual_PPI_V11, ppiV11)
    expect_equal(
      actual_PPI_V11,
      ppiV11,
      info = paste(
        paste(string_db$score_threshold),
        "\n",
        paste(string_db$version),
        "\n",
        paste(string_db$species),
        "\n",
        paste(string_db$input_directory),
        "\n",
        paste(packageVersion("STRINGdb")),
        "\n",
        paste(BiocManager::version())
      ))
  }
}) 
