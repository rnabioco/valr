context("bed_coverage")

x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 0, 10,
  "chr1", 15, 20,
  "chr1", 21, 25
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 3, 15
)

a <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 20, 70, 6, 25, "+",
  "chr1", 50, 100, 1, 25, "-",
  "chr1", 200, 250, 3, 25, "+",
  "chr2", 80, 130, 5, 25, "-",
  "chr2", 150, 200, 4, 25, "+",
  "chr2", 180, 230, 2, 25, "-"
)

b <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 0, 75, 19, 75, "+",
  "chr1", 3, 78, 15, 75, "+",
  "chr1", 74, 149, 7, 75, "+",
  "chr1", 77, 152, 16, 75, "-",
  "chr1", 92, 167, 6, 75, "+",
  "chr1", 116, 191, 9, 75, "+",
  "chr1", 128, 203, 12, 75, "-",
  "chr1", 132, 207, 10, 75, "-",
  "chr1", 143, 218, 20, 75, "-",
  "chr1", 163, 238, 8, 75, "-",
  "chr2", 6, 81, 17, 75, "-",
  "chr2", 39, 114, 4, 75, "+",
  "chr2", 74, 149, 11, 75, "-",
  "chr2", 77, 152, 1, 75, "+",
  "chr2", 114, 189, 3, 75, "-",
  "chr2", 127, 202, 5, 75, "-",
  "chr2", 131, 206, 13, 75, "+",
  "chr2", 137, 212, 2, 75, "-",
  "chr2", 139, 214, 18, 75, "-",
  "chr2", 163, 238, 14, 75, "+"
)

c <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 0, 50,
  "chr1", 12, 20
)

test_that("default coverage works", {
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand, ~.ints, ~.cov, ~.len, ~.frac,
    "chr1", 20, 70, 6, 25, "+", 2, 50, 50, 1.0000000,
    "chr1", 50, 100, 1, 25, "-", 5, 50, 50, 1.0000000,
    "chr1", 200, 250, 3, 25, "+", 4, 38, 50, 0.7600000,
    "chr2", 80, 130, 5, 25, "-", 6, 50, 50, 1.0000000,
    "chr2", 150, 200, 4, 25, "+", 7, 50, 50, 1.0000000,
    "chr2", 180, 230, 2, 25, "-", 6, 50, 50, 1.0000000
  )
  res <- bed_coverage(a, b)
  expect_true(all(res == pred))
})

test_that("coverage of stranded tbls can be calc", {
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand, ~.ints, ~.cov, ~.len, ~.frac,
    "chr1", 20, 70, 6, 25, "+", 2, 50, 50, 1.0000000,
    "chr1", 50, 100, 1, 25, "-", 1, 23, 50, 0.4600000,
    "chr1", 200, 250, 3, 25, "+", 0, 0, 50, 0.000000,
    "chr2", 80, 130, 5, 25, "-", 4, 50, 50, 1.0000000,
    "chr2", 150, 200, 4, 25, "+", 3, 50, 50, 1.0000000,
    "chr2", 180, 230, 2, 25, "-", 4, 34, 50, 0.6800000
  )
  res <- bed_coverage(group_by(a, strand), group_by(b, strand)) %>%
    arrange(chrom, start)
  expect_true(all(res == pred))
})

test_that(" strand_opp coverage works (strand_opp = TRUE)", {
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand, ~.ints, ~.cov, ~.len, ~.frac,
    "chr1", 20, 70, 6, 25, "+", 0, 0, 50, 0.0000000,
    "chr1", 50, 100, 1, 25, "-", 4, 50, 50, 1.0000000,
    "chr1", 200, 250, 3, 25, "+", 4, 38, 50, 0.760000,
    "chr2", 80, 130, 5, 25, "-", 2, 50, 50, 1.0000000,
    "chr2", 150, 200, 4, 25, "+", 4, 50, 50, 1.0000000,
    "chr2", 180, 230, 2, 25, "-", 2, 50, 50, 1.0000000
  )
  res <- bed_coverage(group_by(a, strand), group_by(flip_strands(b), strand)) %>%
    arrange(chrom, start)
  expect_true(all(res == pred))
})

test_that("ensure that coverage is calculated with respect to input tbls issue#108", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 200, "B",
    "chr1", 200, 400, "A",
    "chr1", 500, 600, "C",
    "chr2", 125, 175, "C",
    "chr2", 150, 200, "A",
    "chr3", 100, 300, "A"
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~group,
    "chr1", 100, 199, "A",
    "chr1", 200, 400, "B",
    "chr1", 500, 600, "A",
    "chr2", 125, 175, "C",
    "chr2", 350, 500, "A",
    "chr3", 500, 600, "A"
  )

  x <- arrange(x, chrom, start)
  x <- group_by(x, group, chrom)
  y <- arrange(y, chrom, start)
  y <- group_by(y, group, chrom)

  res <- bed_coverage(x, y)
  expect_true(res[1, ".ints"] == 0)
})

# from https://github.com/arq5x/bedtools2/blob/master/test/coverage/test-coverage.sh
test_that("Test the last record in file with no overlaps is reported", {
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 10,
    "chr1", 15, 20,
    "chr1", 21, 25
  )
  y <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 3, 15
  )
  res <- bed_coverage(x, y)
  expect_equal(res$.cov, c(7, 0, 0))
})

test_that("Test that simple chr	0	100 works", {
  z <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 0, 100
  )
  res <- bed_coverage(z, z)
  expect_equal(nrow(res), 1)
})
