context("bed_subtract")

x <- tibble::tribble(
  ~ chrom, ~ start, ~ end,
  "chr1", 100, 200,
  "chr1", 250, 400,
  "chr1", 500, 600,
  "chr1", 1000, 2000
) %>%
  group_by(chrom)

y <- tibble::tribble(
  ~ chrom, ~ start, ~ end,
  "chr1", 150, 175,
  "chr1", 525, 575,
  "chr1", 1100, 1200,
  "chr1", 1400, 1600
) %>%
  group_by(chrom)

test_that("any = TRUE eliminates overlapping intervals", {
  res <- bed_subtract(x, y, any = TRUE)
  pred <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 250, 400
  )
  expect_equal(res, pred)
})

test_that("fully contained y intervals generate new intervals", {
  res <- bed_subtract(x, y)
  expect_equal(nrow(res), 8)
})

test_that("left dangling y intervals adjust x starts", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 100, 200
  ) %>%
    group_by(chrom)

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 75, 150
  ) %>%
    group_by(chrom)

  res <- bed_subtract(x, y)
  expect_equal(res$start, 150)
})

test_that("right dangling y intervals adjust x ends", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 100, 200
  ) %>%
    group_by(chrom)

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 175, 250
  ) %>%
    group_by(chrom)

  res <- bed_subtract(x, y)
  expect_equal(res$end, 175)
})

test_that("fully contained x intervals are removed", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 100, 200
  ) %>%
    group_by(chrom)

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 50, 250
  ) %>%
    group_by(chrom)

  res <- bed_subtract(x, y)
  expect_equal(nrow(res), 0)
})

test_that("subtractions from x bed_tbl with more chroms than y are captured", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 100, 200,
    "chr3", 400, 500
  )

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr3", 425, 475
  )

  res <- bed_subtract(x, y)
  expect_true("chr3" %in% res$chrom)
})

test_that("non-overlapping intervals from different chrom are not dropped", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 100, 200,
    "chr3", 400, 500
  )

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr3", 425, 475
  )

  res <- bed_subtract(x, y)
  expect_true("chr1" %in% res$chrom)
})

a <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
  "chr1", 10, 20, "a1", 1, "+",
  "chr1", 50, 70, "a2", 2, "-"
)

b <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
  "chr1", 18, 25, "b1", 1, "-",
  "chr1", 80, 90, "b2", 2, "+"
)

test_that("tbls grouped by strand are processed", {
  res <- bed_subtract(group_by(a, strand), group_by(b, strand))
  expect_equal(nrow(res), 2)
  expect_true(all(res == a))
})

test_that("longest merged y intervals are used for subtraction", {
  x <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 500, 600
  )

  y <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 510, 580,
    "chr1", 550, 575
  )

  res <- bed_subtract(x, y)
  expect_true(max(res$start) == 580)
})

test_that("all intervals are dropped in large dataset", {
  bg <- read_bedgraph(valr_example("hela.h3k4.chip.bg.gz"))
  res <- bed_subtract(bg, bg)
  expect_true(nrow(res) == 0)
})

# from https://github.com/arq5x/bedtools2/blob/master/test/subtract/test-subtract.sh
a <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
  "chr1", 10, 20, "a1", 1, "+",
  "chr1", 50, 70, "a2", 2, "-"
)

b <- tibble::tribble(
  ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
  "chr1", 18, 25, "b1", 1, "-",
  "chr1", 80, 90, "b2", 2, "+"
)

test_that("test baseline subtraction", {
  c <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 10, 18, "a1", 1, "+",
    "chr1", 50, 70, "a2", 2, "-"
  )
  res <- bed_subtract(a, b)
  expect_equivalent(res, c)
})

test_that("test any = TRUE subtraction", {
  c <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 50, 70
  )
  res <- bed_subtract(a, b, any = TRUE)
  expect_equal(res, c)
})

test_that("test with 2 DBs", {
  b2 <- tibble::tribble(
    ~ chrom, ~ start, ~ end,
    "chr1", 5, 15,
    "chr1", 55, 65
  )

  c <- tibble::tribble(
    ~ chrom, ~ start, ~ end, ~ name, ~ score, ~ strand,
    "chr1", 15, 18, "a1", 1, "+",
    "chr1", 50, 55, "a2", 2, "-",
    "chr1", 65, 70, "a2", 2, "-"
  )
  res <- bed_subtract(bed_subtract(a, b), b2)
  expect_equivalent(res, c)
})
