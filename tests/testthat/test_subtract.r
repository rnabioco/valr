context("bed_subtract")

x <- tibble::frame_data(
 ~chrom, ~start, ~end,
 "chr1", 100,    200,
 "chr1", 250,    400,
 "chr1", 500,    600,
 "chr1", 1000,   2000
) %>% group_by(chrom)

y <- tibble::frame_data(
 ~chrom, ~start, ~end,
 "chr1", 150,    175,
 "chr1", 525,    575,
 "chr1", 1100,   1200,
 "chr1", 1400,   1600
) %>% group_by(chrom)

test_that("any = TRUE eliminates overlapping intervals", {
  res <- bed_subtract(x, y, any = TRUE)
  expect_equal(nrow(res), 1)
})

test_that("fully contained y intervals generate new intervals", {
  res <- bed_subtract(x, y)
  expect_equal(nrow(res), 8)
})

test_that("left dangling y intervals adjust x starts", {
  x <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  ) %>% group_by(chrom)
  
  y <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 75,    150
  ) %>% group_by(chrom)
  
  res <- bed_subtract(x, y)
  expect_equal(res$start, 150)
})

test_that("right dangling y intervals adjust x ends", {
  x <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  ) %>% group_by(chrom)
  
  y <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 175,    250
  ) %>% group_by(chrom)
  
  res <- bed_subtract(x, y)
  expect_equal(res$end, 175)
})

test_that("fully contained x intervals are removed", {
  x <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 100,    200
  ) %>% group_by(chrom)
  
  y <- tibble::frame_data(
    ~chrom, ~start, ~end,
    "chr1", 50,    250
  ) %>% group_by(chrom)
  
  res <- bed_subtract(x, y)
  expect_equal(nrow(res), 0)
  
})

