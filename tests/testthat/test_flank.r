context("bed_flank")

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
  out <- bed_flank(x, genome, left = dist)
  expect_true(all(x$start - out$start == dist))
})

test_that("right arg works", {
  dist <- 100
  out <- bed_flank(x, genome, right = dist)
  expect_true(all(out$end - x$end == dist))
})

test_that("both arg works", {
  dist <- 100
  out <- bed_flank(x, genome, both = dist)
  expect_equal(nrow(out), 4)
})

test_that("args with left and right works", {
  dist_l <- 100
  dist_r <- 50
  out <- bed_flank(x, genome, left = dist_l, right = dist_r)
  expect_equal(nrow(out), 4)
  # test left side
  expect_true(x$start[1] - out$start[1] == 100)
  # test right side
  expect_true(out$end[3] - x$end[1] == 50)
})

test_that("all left and right intervals are reported with both arg", {
  dist <- 100
  out_left <- bed_flank(x, genome, left = dist)
  out_right <- bed_flank(x, genome, right = dist)
  out_both <- bed_flank(x, genome, both = dist)
  out_left_right <- dplyr::bind_rows(out_left, out_right) %>% arrange(chrom, start)
  expect_true(all(out_both == out_left_right))
})

test_that("fraction arg works", {
  dist <- 0.1
  out <- bed_flank(x, genome, both = dist, fraction = TRUE)
  expect_true(nrow(out) == 4)
  expect_true(all(out$end - out$start == 50))
})

test_that("strand arg with both works", {
  dist <- 100
  out <- bed_flank(x, genome, both = dist, strand = TRUE)
  out_nostrand <- bed_flank(x, genome, both = dist)
  expect_true(nrow(out) == 4)
  expect_true(all(out == out_nostrand))
})

test_that("strand arg with left works", {
  dist <- 100
  out <- bed_flank(x, genome, left = dist, strand = TRUE)
  expect_true(nrow(out) == 2)
  plus_in <- dplyr::filter(x, strand == "+")
  plus_out <- dplyr::filter(out, strand == "+")
  minus_in <- dplyr::filter(x, strand == "-")
  minus_out <- dplyr::filter(out, strand == "-")
  expect_true(plus_in$start > plus_out$start)
  expect_true(minus_in$start < minus_out$start)
})

test_that("strand arg with right works", {
  dist <- 100
  out <- bed_flank(x, genome, right = dist, strand = TRUE)
  expect_true(nrow(out) == 2)
  plus_in <- dplyr::filter(x, strand == "+")
  plus_out <- dplyr::filter(out, strand == "+")
  minus_in <- dplyr::filter(x, strand == "-")
  minus_out <- dplyr::filter(out, strand == "-")
  expect_true(plus_in$start < plus_out$start)
  expect_true(minus_in$start > minus_out$start)
})

test_that("strand arg with left and right works", {
  dist_l <- 100
  dist_r <- 50
  out <- bed_flank(x, genome, left = dist_l, right = dist_r, strand = TRUE)
  expect_true(nrow(out) == 4)
  # test left side plus strand
  expect_true(x$start[1] - out$start[1] == 100)
  # test right side plus strand
  expect_true(out$end[3] - x$end[1] == 50)
  # test left side minus strand
  expect_true(x$start[2] - out$start[2] == 50)
  # test right side minus strand
  expect_true(out$end[4] - x$end[2] == 100)
})

test_that("strand arg with both and fraction works", {
  dist <- 0.2
  out <- bed_flank(x, genome, both = dist, strand = TRUE, fraction = TRUE)
  out_nostrand <- bed_flank(x, genome, both = dist, fraction = TRUE)
  expect_true(nrow(out) == 4)
  expect_true(all(out == out_nostrand))
})

test_that("strand arg with left and fraction works", {
  dist <- 0.2
  out <- bed_flank(x, genome, left = dist, strand = TRUE, fraction = TRUE)
  expect_true(nrow(out) == 2)
  # test left side plus strand
  expect_true(x$start[1] - out$start[1] == 100)
  # test left side minus strand
  expect_true(out$end[2] - x$end[2] == 100)
})

test_that("strand arg with right and fraction works", {
  dist <- 0.2
  out <- bed_flank(x, genome, right = dist, strand = TRUE, fraction = TRUE)
  out <- dplyr::arrange(out, desc(start))
  expect_true(nrow(out) == 2)
  # test right side plus strand
  expect_true(out$end[1] - x$end[1] == 100)
  # test right side minus strand
  expect_true(x$start[2] - out$start[2] == 100)
})

test_that("strand arg with left and right and fraction works", {
  dist_l <- 0.2
  dist_r <- 0.1
  out <- bed_flank(x, genome, left = dist_l, right = dist_r, strand = TRUE, fraction = TRUE)
  expect_true(nrow(out) == 4)
  # test left side plus strand
  expect_true(x$start[1] - out$start[1] == 100)
  # test right side plus strand
  expect_true(out$end[3] - x$end[1] == 50)
  # test right side minus strand
  expect_true(x$start[2] - out$start[2] == 50)
  # test left side minus strand
  expect_true(out$end[4] - x$end[2] == 100)
})

test_that("intervals are not reported off of chromosomes", {
    dist <- 600
    out <- bed_flank(x, genome, left = dist)
    #test left side
    expect_true(nrow(out) == 1)
    expect_true(out$start[1] == 400)
    #test right side
    dist <- 3501
    out <- bed_flank(x, genome, right = dist)
    expect_true(nrow(out) == 1)
    expect_true(out$end[1] == 4501)
  })
