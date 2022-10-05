genome_path <- system.file("extdata", "hg19.chrom.sizes.gz", package = "valr")
genome <- read_genome(genome_path)

test_that("genomes are correctly read", {
  expect_equal(nrow(genome), 25)
  expect_equal(colnames(genome), c("chrom", "size"))
  expect_true("data.frame" %in% class(genome))
})

test_that("unbounded intervals are removed", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", -100, 500,
    "chr1", 100, 1e9,
    "chr1", 500, 1000
  )
  expect_equal(nrow(bound_intervals(x, genome)), 1)
})

test_that("trim param removes dangling intervals", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 249250721
  )
  res <- bound_intervals(x, genome, trim = TRUE)
  expect_equal(res$end, 249250621)
})

test_that("duplicate chroms throw an error.", {
  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 1000,
    "chr1", 1000
  )

  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1, 10000
  )

  expect_error(bed_complement(x, genome))
})
