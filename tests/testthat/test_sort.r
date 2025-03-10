test_that("intervals can be sorted by size", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 2000,
    "chr1", 200, 400,
    "chr2", 100, 200,
    "chr1", 400, 700
  )
  res <- bed_sort(x, by_size = TRUE)
  expect_equal(res$start, c(100, 200, 400, 1000))
})

test_that("intervals can be sorted by size and chrom", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 2000,
    "chr1", 200, 400,
    "chr2", 100, 200,
    "chr1", 400, 700
  )
  res <- bed_sort(x, by_size = TRUE, by_chrom = TRUE)
  expect_equal(res$start, c(200, 400, 1000, 100))
})

test_that("intervals can be reverse sorted by size", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 2000,
    "chr1", 200, 400,
    "chr2", 100, 200,
    "chr1", 400, 700
  )
  res <- bed_sort(x, by_size = TRUE, reverse = TRUE)
  expect_equal(res$start, c(1000, 400, 200, 100))
})

test_that("intervals can be reverse sorted by size and chrom", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 2000,
    "chr1", 200, 400,
    "chr2", 100, 200,
    "chr1", 400, 700
  )
  res <- bed_sort(x, by_size = TRUE, by_chrom = TRUE, reverse = TRUE)
  expect_equal(res$start, c(1000, 400, 200, 100))
})

test_that("intervals can be sorted by chrom", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 2000,
    "chr2", 100, 200,
    "chr1", 200, 400,
    "chr3", 400, 700
  )
  res <- bed_sort(x, by_chrom = TRUE)
  expect_equal(res$start, c(200, 1000, 100, 400))
})

test_that("intervals can be reverse sorted by start and chrom", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 2000,
    "chr2", 100, 200,
    "chr1", 200, 400,
    "chr3", 400, 700
  )
  res <- bed_sort(x, by_chrom = TRUE, reverse = TRUE)
  expect_equal(res$start, c(1000, 200, 100, 400))
})

test_that("ties in start are sorted by end", {
  # fmt: skip
  x <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 2000,
    "chr1", 1000, 400
  )

  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end,
    "chr1", 1000, 400,
    "chr1", 1000, 2000
  )

  res <- bed_sort(x)
  expect_equal(res, pred, ignore_attr = TRUE)
})

# from https://github.com/arq5x/bedtools2/blob/master/test/sort/test-sort.sh
# fmt: skip
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr7", 240L, 560L, "x", "86", "-",
  "chr7", 210L, 525L, "d", "21", "+",
  "chr7", 2100L, 2310L, "e", "32", "+",
  "chr3", 40L, 260L, "p", "41", "-",
  "chr3", 10L, 220L, "f", "12", "+",
  "chr3", 100L, 320L, "g", "96", "-",
  "chr9", 140L, 160L, "z", "05", "+",
  "chr9", 110L, 120L, "a", "81", "-",
  "chr9", 1100L, 1120L, "b", "12", "+"
)

test_that("Test default", {
  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr3", 10L, 220L, "f", "12", "+",
    "chr3", 40L, 260L, "p", "41", "-",
    "chr3", 100L, 320L, "g", "96", "-",
    "chr7", 210L, 525L, "d", "21", "+",
    "chr7", 240L, 560L, "x", "86", "-",
    "chr7", 2100L, 2310L, "e", "32", "+",
    "chr9", 110L, 120L, "a", "81", "-",
    "chr9", 140L, 160L, "z", "05", "+",
    "chr9", 1100L, 1120L, "b", "12", "+"
  )

  res <- bed_sort(x)
  expect_equal(res$start, pred$start)
})

test_that("Test by_size = TRUE", {
  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr9", 110L, 120L, "a", "81", "-",
    "chr9", 140L, 160L, "z", "05", "+",
    "chr9", 1100L, 1120L, "b", "12", "+",
    "chr3", 10L, 220L, "f", "12", "+",
    "chr7", 2100L, 2310L, "e", "32", "+",
    "chr3", 40L, 260L, "p", "41", "-",
    "chr3", 100L, 320L, "g", "96", "-",
    "chr7", 210L, 525L, "d", "21", "+",
    "chr7", 240L, 560L, "x", "86", "-"
  )

  res <- bed_sort(x, by_size = TRUE)
  expect_equal(res$start - res$end, pred$start - pred$end)
})

test_that("Test by_size = TRUE, reverse = TRUE", {
  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr7", 240L, 560L, "x", "86", "-",
    "chr7", 210L, 525L, "d", "21", "+",
    "chr3", 40L, 260L, "p", "41", "-",
    "chr3", 100L, 320L, "g", "96", "-",
    "chr3", 10L, 220L, "f", "12", "+",
    "chr7", 2100L, 2310L, "e", "32", "+",
    "chr9", 140L, 160L, "z", "05", "+",
    "chr9", 1100L, 1120L, "b", "12", "+",
    "chr9", 110L, 120L, "a", "81", "-"
  )

  res <- bed_sort(x, by_size = TRUE, reverse = TRUE)
  expect_equal(res$start - res$end, pred$start - pred$end)
})

test_that("Test by_chrom = TRUE", {
  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr3", 10L, 220L, "f", "12", "+",
    "chr3", 40L, 260L, "p", "41", "-",
    "chr3", 100L, 320L, "g", "96", "-",
    "chr7", 210L, 525L, "d", "21", "+",
    "chr7", 240L, 560L, "x", "86", "-",
    "chr7", 2100L, 2310L, "e", "32", "+",
    "chr9", 110L, 120L, "a", "81", "-",
    "chr9", 140L, 160L, "z", "05", "+",
    "chr9", 1100L, 1120L, "b", "12", "+"
  )

  res <- bed_sort(x, by_chrom = TRUE)
  expect_equal(res$start, pred$start)
})

test_that("Test by_size = TRUE, by_chrom = TRUE", {
  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr3", 10L, 220L, "f", "12", "+",
    "chr3", 40L, 260L, "p", "41", "-",
    "chr3", 100L, 320L, "g", "96", "-",
    "chr7", 2100L, 2310L, "e", "32", "+",
    "chr7", 210L, 525L, "d", "21", "+",
    "chr7", 240L, 560L, "x", "86", "-",
    "chr9", 110L, 120L, "a", "81", "-",
    "chr9", 140L, 160L, "z", "05", "+",
    "chr9", 1100L, 1120L, "b", "12", "+"
  )

  res <- bed_sort(x, by_chrom = TRUE, by_size = TRUE)
  expect_equal(res$start, pred$start)
})

test_that("Test by_size = TRUE, by_chrom = TRUE, reverse = TRUE", {
  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr3", 40L, 260L, "p", "41", "-",
    "chr3", 100L, 320L, "g", "96", "-",
    "chr3", 10L, 220L, "f", "12", "+",
    "chr7", 240L, 560L, "x", "86", "-",
    "chr7", 210L, 525L, "d", "21", "+",
    "chr7", 2100L, 2310L, "e", "32", "+",
    "chr9", 140L, 160L, "z", "05", "+",
    "chr9", 1100L, 1120L, "b", "12", "+",
    "chr9", 110L, 120L, "a", "81", "-"
  )

  res <- bed_sort(x, by_chrom = TRUE, by_size = TRUE, reverse = TRUE)
  expect_equal(res$start, pred$start)
})

test_that("Test by_size = TRUE, by_chrom = TRUE, reverse = TRUE", {
  # fmt: skip
  y <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 30L, 40L, "RegionD", 0L, "+",
    "chr1", 25L, 25L, "RegionC", 0L, "+",
    "chr1", 10L, 20L, "RegionA", 0L, "+",
    "chr1", 24L, 30L, "RegionB", 0L, "+"
  )
  # fmt: skip
  pred <- tibble::tribble(
    ~chrom, ~start, ~end, ~name, ~score, ~strand,
    "chr1", 10L, 20L, "RegionA", 0L, "+",
    "chr1", 24L, 30L, "RegionB", 0L, "+",
    "chr1", 25L, 25L, "RegionC", 0L, "+",
    "chr1", 30L, 40L, "RegionD", 0L, "+"
  )

  res <- bed_sort(y)
  expect_equal(res$start, pred$start)
})
