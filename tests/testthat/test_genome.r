library(Rbedtools)
context("genome files")

genome_path <- system.file('extdata', 'hg19.chrom.sizes.gz', package = 'Rbedtools')
genome <- read_genome(genome_path)

test_that("genomes are correctly read", {
  expect_equal(nrow(genome), 93)
  expect_equal(colnames(genome), c('chrom','size'))
  expect_is(genome, "data.frame")
})
