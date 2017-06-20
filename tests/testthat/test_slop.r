context("bed_slop")

genome <- tibble::tribble(
 ~chrom, ~size,
 "chr1", 5000
)

x <- tibble::tribble(
 ~chrom, ~start, ~end, ~name, ~score, ~strand,
 "chr1", 500,    1000, ".",   ".",     "+",
 "chr1", 1000,   1500, ".",   ".",     "-"
)

test_that("left arg works", {
  dist <- 100
  out <- x %>% bed_slop(genome, left = dist)
  expect_true(all(x$start - out$start == dist))
})

test_that("right arg works", {
  dist <- 100
  out <- x %>% bed_slop(genome, right = dist)
  expect_true(all(out$end - x$end == dist))
})

test_that("both arg works", {
  dist <- 100
  out <- x %>% bed_slop(genome, both = dist)
  expect_true(all(x$start - out$start == dist))
  expect_true(all(out$end - x$end == dist))
})

test_that("both with fraction works", {
  res <- bed_slop(x, genome, both = 0.5, fraction = TRUE)
  expect_equal(res$start, c(250, 750))
  expect_equal(res$end, c(1250, 1750))
})

test_that("left / right with fraction works", {
  res <- bed_slop(x, genome, left = 0.5, fraction = TRUE)
  expect_equal(res$start, c(250, 750))
  expect_equal(res$end, c(1000, 1500))
})

test_that("left, fraction, strand works", {
  res <- bed_slop(x, genome, left = 0.5, fraction = TRUE, strand = TRUE)
  expect_equal(res$start, c(250, 1000))
  expect_equal(res$end, c(1000, 1750))
})

test_that("right, fraction, strand works", {
  res <- bed_slop(x, genome, right = 0.5, fraction = TRUE, strand = TRUE)
  expect_equal(res$start, c(500, 750))
  expect_equal(res$end, c(1250, 1500))
})

test_that("strand with left works", {
  res <- bed_slop(x, genome, left = 100, strand = TRUE)
  expect_equal(res$start, c(400, 1000))
  expect_equal(res$end, c(1000, 1600))
})
