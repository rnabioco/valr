context("bed_absdist")

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 10000,
  "chr2", 10000,
  "chr3", 10000
)

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 75, 125
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 50, 100,
  "chr1", 100, 150
)

test_that("absdist calculation is correct", {
  res <- bed_absdist(x, y, genome)
  expect_true(res$.absdist_scaled == 0.005)
})

test_that("self absdist is 0", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 5, 15,
    "chr1", 50, 150,
    "chr2", 1000, 2000,
    "chr3", 3000, 4000
  )

  res <- bed_absdist(x, x, genome)
  expect_true(sum(res$.absdist) == 0)
})

test_that("x ivls without matching y-ivls chroms are reported with absdist = NA", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 5, 15,
    "chr1", 50, 150,
    "chr2", 1000, 2000,
    "chr3", 3000, 4000
  )

  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 25, 125,
    "chr1", 150, 250,
    "chr1", 550, 580,
    "chr2", 1, 1000,
    "chr2", 2000, 3000
  )
  res <- bed_absdist(x, y, genome)
  expect_true("chr3" %in% res$chrom)
  expect_true(any(is.na(res$.absdist)))
  expect_true(any(is.na(res$.absdist_scaled)))
})


test_that("ensure that absdist is calculated with respect to input tbls issue#108", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "A",
    "chr1", 200, 400, "B",
    "chr1", 300, 500, "C",
    "chr2", 125, 175, "C",
    "chr2", 150, 200, "A"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "C",
    "chr1", 200, 400, "B",
    "chr1", 300, 500, "A",
    "chr2", 125, 175, "C",
    "chr2", 350, 500, "A"
  )

  genome <- tibble::tribble(
    ~chrom, ~size,
    "chr1", 10000,
    "chr2", 10000
  )

  orig <- bed_absdist(x, y, genome)
  x <- arrange(x, chrom, start)
  x <- group_by(x, group, chrom)
  y <- arrange(y, chrom, start)
  y <- group_by(y, group, chrom)

  res <- bed_absdist(x, y, genome)

  expect_true(any(orig[, c(1:5)] != res[, c(1:5)]))
})
