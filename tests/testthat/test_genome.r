context("genome file functions")

genome_path <- system.file('extdata', 'hg19.chrom.sizes.gz', package = 'Rbedtools')
genome <- read_genome(genome_path)

test_that("genomes are correctly read", {
  expect_equal(nrow(genome), 93)
  expect_equal(colnames(genome), c('chrom','size'))
  expect_is(genome, "data.frame")
})

test_that("unbounded intervals are removed", {
  
  bed_tbl <- dplyr::tibble(
   ~chrom, ~start, ~end,
   "chr1", -100,   500,
   "chr1", 100,    1e9,
   "chr1", 500,    1000
  )
  expect_equal(nrow(bound_intervals(bed_tbl, genome)), 1)
})
