context('bed_flank')

genome <- tibble::frame_data(
 ~chrom, ~size,
 "chr1", 5000
)

x <- tibble::frame_data(
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

# test fraction
# test off chrom
