context("db")

has_internet <- function() {
  ! is.null(curl::nslookup("r-project.org", error = FALSE))
}

test_that("ucsc connection works", {
  skip("run db tests manually.")
  skip_if_not_installed("curl")
  skip_if_not_installed("RMySQL")

  if (!has_internet()) skip("no internet connection")

  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  db <- db_ucsc("hg38")
  refgene <- tbl(db, "refGene")
  expect_is(refgene, "tbl_mysql")
  expect_equal(ncol(refgene), 16)
})

test_that("ensembl connection works", {
  skip("run db tests manually.")
  skip_if_not_installed("curl")
  skip_if_not_installed("RMySQL")

  if (!has_internet()) skip("no internet connection")

  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  db <- db_ensembl("spermophilus_tridecemlineatus_core_67_2")
  gene <- tbl(db, "gene")
  expect_is(gene, "tbl_mysql")
  expect_equal(ncol(gene), 18)
})
