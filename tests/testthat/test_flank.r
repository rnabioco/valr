context('bed_flank')

genome <- tibble::tribble(
 ~chrom, ~size,
 "chr1", 5000
)

x <- tibble::tribble(
 ~chrom, ~start, ~end, ~name, ~score, ~strand,
 "chr1", 500,    1000, '.',   '.',     '+',
 "chr1", 1000,   1500, '.',   '.',     '-'
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

test_that("all left and right intervals are reported with both arg", {
  dist <- 100
  out_left <- bed_flank(x, genome, left = dist)
  out_right <- bed_flank(x, genome, right = dist)
  out_both <- bed_flank(x, genome, both = dist)
  out_left_right <- dplyr::bind_rows(out_left, out_right) %>% bed_sort()
  expect_true(all(out_both == out_left_right))
})
# test fraction
# test off chrom
