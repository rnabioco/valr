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
