context("bed_jaccard")

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 10,     20,
  "chr1", 30,     40
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 15,     20
)

test_that("jaccard coeff is calculated correctly", {
  res <- bed_jaccard(x, y)
  expect_equal(res$jaccard, 0.25)
})

test_that("jaccard coeff is calc'd for large data sets", {
  genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

  x <- bed_random(genome, n = 1e5, seed = 10000)
  y <- bed_random(genome, n = 1e5, seed = 20000)

  res <- bed_jaccard(x, y)
  expect_equal(round(res$jaccard, 3), 0.016)
})

test_that("jaccard with grouped inputs are calculated", {
   genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

  x <- bed_random(genome, n = 1e5, seed = 10000)
  y <- bed_random(genome, n = 1e5, seed = 20000)

  res <- bed_jaccard(group_by(x, chrom),
                     group_by(y, chrom))

  expect_equal(nrow(res), 24)
  expect_true('chrom' %in% names(res))
})
