context("bed12_to_exons")

test_that("invalid ncol causes an error", {
  x <- tibble::tribble()
  expect_error(bed12_to_exons(x))
})

test_that("BED12 is parsed correctly", {
  x <- read_bed12(valr_example("mm9.refGene.bed.gz"))
  expect_equal(nrow(bed12_to_exons(x)), 1683)
  expect_equal(ncol(bed12_to_exons(x)), 6)
})
