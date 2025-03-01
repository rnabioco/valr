# fmt: skip
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 10, 20,
  "chr1", 30, 40,
  "chr1", 51, 52
)

# fmt: skip
y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 15, 25,
  "chr1", 51, 52
)

# fmt: skip
genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 500
)

test_that("fisher p.value is correct", {
  res <- bed_fisher(x, y, genome)
  expect_equal(round(res$p.value, 5), 0.00385)
})
